#ifndef IGA_TEMPORAL_UNSTRUCTURED_VTK_HDF_HPP
#define IGA_TEMPORAL_UNSTRUCTURED_VTK_HDF_HPP

#include "PartitionedVtkOutput.hpp"
#include "TemporalVtkHdf.hpp"

namespace iga {

// One complete grid per file. Topology and Int64 IDs are static; point and
// cell fields are losslessly compressed FP64. Moving coordinates are optional.
// This is a serial writer: replicated native solvers call it on root only.
class TemporalUnstructuredVtkHdfWriter
{
public:
	TemporalUnstructuredVtkHdfWriter(const std::filesystem::path& path,
		const VtkPartition& initial,bool moving=false,int compression=4,
		std::int64_t resume_steps=0)
		: moving_(moving),compression_(compression),points_(initial.point_ids.size()),
		  cells_(initial.cell_ids.size()),topology_hash_(TopologyHash(initial,moving)),
		  point_schema_(hdf_detail::ArraySchema(initial.point_arrays)),
		  cell_schema_(hdf_detail::ArraySchema(initial.cell_arrays))
	{
		Validate(initial,0.);
		if(points_==0||cells_==0||compression<0||compression>9||resume_steps<0)
			throw std::invalid_argument("invalid temporal unstructured VTKHDF configuration");
		if(!path.parent_path().empty())std::filesystem::create_directories(path.parent_path());
		file_=hdf_detail::RequireHandle(resume_steps
			?H5Fopen(path.string().c_str(),H5F_ACC_RDWR,H5P_DEFAULT)
			:H5Fcreate(path.string().c_str(),H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT),
			H5Fclose,"cannot open exclusive/resumed VTKHDF output: "+path.string());
		if(resume_steps)OpenExisting(initial,resume_steps);
		else Create(initial);
	}

	void Append(double time,const VtkPartition& piece)
	{
		if(closing_||failed_)throw std::runtime_error("VTKHDF writer is closed or failed");
		Validate(piece,time);
		if(time<=last_time_)throw std::invalid_argument("VTKHDF times must strictly increase");
		failed_=true;
		const auto write=[&](const std::string& path,hsize_t offset,hsize_t rows,
			int components,const auto* data){
			auto dataset=hdf_detail::RequireHandle(H5Dopen2(root_.get(),path.c_str(),
				H5P_DEFAULT),H5Dclose,"cannot open VTKHDF dataset: "+path);
			hdf_detail::WriteRows(dataset.get(),offset,rows,components,data);
		};
		if(moving_)write("Points",steps_*points_,points_,3,piece.grid.points.data());
		for(bool cell:{false,true}){
			const std::string group=cell?"CellData":"PointData";
			const auto count=cell?cells_:points_;
			const std::int64_t offset=steps_*count,zero=0;
			for(const auto& array:cell?piece.cell_arrays:piece.point_arrays){
				write(group+"/"+array.name,offset,count,array.components,array.values.data());
				write("Steps/"+group+"Offsets/"+array.name,steps_,1,1,&offset);
			}
			write("Steps/"+group+"Offsets/"+(cell?"GlobalCellIds":"GlobalPointIds"),
				steps_,1,1,&zero);
		}
		const std::int64_t zero=0,one=1,point_offset=moving_?steps_*points_:0;
		write("Steps/Values",steps_,1,1,&time);
		write("Steps/PointOffsets",steps_,1,1,&point_offset);
		for(const auto* name:{"PartOffsets","CellOffsets","ConnectivityIdOffsets"})
			write(std::string("Steps/")+name,steps_,1,1,&zero);
		write("Steps/NumberOfParts",steps_,1,1,&one);
		hdf_detail::Require(H5Fflush(file_.get(),H5F_SCOPE_GLOBAL),"cannot flush VTKHDF data");
		auto steps=hdf_detail::OpenGroup(root_.get(),"Steps");
		hdf_detail::ReplaceScalarAttribute<std::int64_t>(steps.get(),"NSteps",steps_+1);
		hdf_detail::Require(H5Fflush(file_.get(),H5F_SCOPE_GLOBAL),"cannot commit VTKHDF step");
		++steps_;last_time_=time;failed_=false;
	}

	void Close()
	{
		closing_=true;
		if(!file_)return;
		hdf_detail::Require(H5Fflush(file_.get(),H5F_SCOPE_GLOBAL),"cannot finalize VTKHDF");
		root_.CloseChecked("cannot close VTKHDF root");
		if(H5Fget_obj_count(file_.get(),H5F_OBJ_ALL|H5F_OBJ_LOCAL)!=1)
			throw std::runtime_error("VTKHDF has outstanding handles");
		file_.CloseChecked("cannot close VTKHDF file");
	}

	double LastTime() const { return last_time_; }

private:
	static std::uint64_t TopologyHash(const VtkPartition& piece,bool moving)
	{
		std::uint64_t hash=1469598103934665603ULL;
		const auto add=[&](const auto& values){
			const std::uint64_t size=values.size();
			hash=hdf_detail::HashBytes(hash,&size,sizeof(size));
			hash=hdf_detail::HashBytes(hash,values.data(),values.size()*sizeof(values[0]));
		};
		add(piece.point_ids);add(piece.cell_ids);add(piece.grid.connectivity);
		add(piece.grid.offsets);add(piece.grid.types);
		if(!moving)add(piece.grid.points);
		return hash;
	}
	void Validate(const VtkPartition& piece,double time) const
	{
		ValidateVtkPartition(piece,time);
		if(TopologyHash(piece,moving_)!=topology_hash_
			||hdf_detail::ArraySchema(piece.point_arrays)!=point_schema_
			||hdf_detail::ArraySchema(piece.cell_arrays)!=cell_schema_)
			throw std::invalid_argument("VTKHDF topology, static coordinates or schema changed");
		for(const auto* arrays:{&piece.point_arrays,&piece.cell_arrays})
			for(const auto& array:*arrays)
				if(array.name.find('/')!=std::string::npos||array.name=="."||array.name=="..")
					throw std::invalid_argument("invalid VTKHDF array name");
	}
	void Create(const VtkPartition& piece)
	{
		root_=hdf_detail::CreateGroup(file_.get(),"VTKHDF");
		const hsize_t two=2;
		auto space=hdf_detail::RequireHandle(H5Screate_simple(1,&two,nullptr),H5Sclose,"version space");
		auto attribute=hdf_detail::RequireHandle(H5Acreate2(root_.get(),"Version",H5T_NATIVE_INT64,
			space.get(),H5P_DEFAULT,H5P_DEFAULT),H5Aclose,"version attribute");
		const std::int64_t version[2]={2,1};
		hdf_detail::Require(H5Awrite(attribute.get(),H5T_NATIVE_INT64,version),"write version");
		hdf_detail::WriteStringAttribute(root_.get(),"Type","UnstructuredGrid");
		const auto fixed=[&](const std::string& name,const auto& values){
			hdf_detail::WriteFixedDataset(root_.get(),name,values.data(),{values.size()},compression_);
		};
		const std::int64_t connectivity=piece.grid.connectivity.size();
		hdf_detail::WriteFixedDataset(root_.get(),"NumberOfPoints",&points_,{1},0);
		hdf_detail::WriteFixedDataset(root_.get(),"NumberOfCells",&cells_,{1},0);
		hdf_detail::WriteFixedDataset(root_.get(),"NumberOfConnectivityIds",&connectivity,{1},0);
		fixed("Connectivity",piece.grid.connectivity);
		std::vector<std::int64_t> offsets{0};
		offsets.insert(offsets.end(),piece.grid.offsets.begin(),piece.grid.offsets.end());
		fixed("Offsets",offsets);
		fixed("Types",std::vector<std::uint8_t>(piece.grid.types.begin(),piece.grid.types.end()));
		if(moving_)hdf_detail::CreateExpandableDataset<double>(root_.get(),"Points",3,compression_);
		else hdf_detail::WriteFixedDataset(root_.get(),"Points",piece.grid.points.data(),
			{static_cast<hsize_t>(points_),3},compression_);
		auto steps=hdf_detail::CreateGroup(root_.get(),"Steps");
		hdf_detail::WriteScalarAttribute<std::int64_t>(steps.get(),"NSteps",0);
		hdf_detail::CreateExpandableDataset<double>(steps.get(),"Values",1,0);
		for(const auto* name:{"PartOffsets","NumberOfParts","PointOffsets","CellOffsets","ConnectivityIdOffsets"})
			hdf_detail::CreateExpandableDataset<std::int64_t>(steps.get(),name,1,0,
				std::string(name)=="CellOffsets"||std::string(name)=="ConnectivityIdOffsets");
		for(bool cell:{false,true}){
			const std::string group=cell?"CellData":"PointData";
			auto data=hdf_detail::CreateGroup(root_.get(),group);
			auto positions=hdf_detail::CreateGroup(steps.get(),group+"Offsets");
			const auto* id=cell?"GlobalCellIds":"GlobalPointIds";
			const auto& ids=cell?piece.cell_ids:piece.point_ids;
			hdf_detail::WriteFixedDataset(data.get(),id,ids.data(),{ids.size()},compression_);
			hdf_detail::CreateExpandableDataset<std::int64_t>(positions.get(),id,1,0);
			for(const auto& array:cell?piece.cell_arrays:piece.point_arrays){
				hdf_detail::CreateExpandableDataset<double>(data.get(),array.name,array.components,compression_);
				hdf_detail::CreateExpandableDataset<std::int64_t>(positions.get(),array.name,1,0);
			}
		}
		auto metadata=hdf_detail::CreateGroup(file_.get(),"TubularFlowIGA");
		hdf_detail::WriteScalarAttribute<std::uint64_t>(metadata.get(),"TopologyHash",topology_hash_);
		hdf_detail::WriteScalarAttribute<std::int64_t>(metadata.get(),"MovingPoints",moving_);
		hdf_detail::WriteStringAttribute(metadata.get(),"PointArraySchema",point_schema_);
		hdf_detail::WriteStringAttribute(metadata.get(),"CellArraySchema",cell_schema_);
	}
	void OpenExisting(const VtkPartition& piece,std::int64_t expected_steps)
	{
		root_=hdf_detail::OpenGroup(file_.get(),"VTKHDF");
		auto metadata=hdf_detail::OpenGroup(file_.get(),"TubularFlowIGA");
		auto steps=hdf_detail::OpenGroup(root_.get(),"Steps");
		steps_=hdf_detail::ReadScalarAttribute<std::int64_t>(steps.get(),"NSteps");
		if(steps_!=expected_steps
			||hdf_detail::ReadScalarAttribute<std::uint64_t>(metadata.get(),"TopologyHash")!=topology_hash_
			||hdf_detail::ReadScalarAttribute<std::int64_t>(metadata.get(),"MovingPoints")!=moving_
			||hdf_detail::ReadStringAttribute(metadata.get(),"PointArraySchema")!=point_schema_
			||hdf_detail::ReadStringAttribute(metadata.get(),"CellArraySchema")!=cell_schema_)
			throw std::runtime_error("VTKHDF resume identity or published step count differs");
		const auto dimensions=[&](const std::string& path,hsize_t rows,int components,bool matrix=false){
			auto dataset=hdf_detail::RequireHandle(H5Dopen2(root_.get(),path.c_str(),H5P_DEFAULT),H5Dclose,"resume dataset");
			auto space=hdf_detail::RequireHandle(H5Dget_space(dataset.get()),H5Sclose,"resume space");
			const int rank=H5Sget_simple_extent_ndims(space.get());
			hsize_t sizes[2]={0,1};
			if(rank!=(components>1||matrix?2:1)||H5Sget_simple_extent_dims(space.get(),sizes,nullptr)<0
				||sizes[0]!=rows||sizes[1]!=static_cast<hsize_t>(components))
				throw std::runtime_error("VTKHDF resume has incomplete or trailing rows: "+path);
		};
		dimensions("Points",points_*(moving_?steps_:1),3);
		for(const auto* name:{"Values","PartOffsets","NumberOfParts","PointOffsets","CellOffsets","ConnectivityIdOffsets"})
			dimensions(std::string("Steps/")+name,steps_,1,
				std::string(name)=="CellOffsets"||std::string(name)=="ConnectivityIdOffsets");
		for(bool cell:{false,true}){
			const std::string group=cell?"CellData":"PointData";
			const auto* id=cell?"GlobalCellIds":"GlobalPointIds";
			dimensions(group+"/"+id,cell?cells_:points_,1);
			dimensions("Steps/"+group+"Offsets/"+id,steps_,1);
			for(const auto& array:cell?piece.cell_arrays:piece.point_arrays){
				dimensions(group+"/"+array.name,steps_*(cell?cells_:points_),array.components);
				dimensions("Steps/"+group+"Offsets/"+array.name,steps_,1);
			}
		}
		auto values=hdf_detail::RequireHandle(H5Dopen2(steps.get(),"Values",H5P_DEFAULT),H5Dclose,"resume time");
		last_time_=hdf_detail::ReadRowValue<double>(values.get(),steps_-1);
		if(!std::isfinite(last_time_))throw std::runtime_error("nonfinite VTKHDF resume time");
	}
	bool moving_=false,closing_=false,failed_=false;
	int compression_=4;
	std::int64_t points_=0,cells_=0,steps_=0;
	std::uint64_t topology_hash_=0;
	std::string point_schema_,cell_schema_;
	double last_time_=-std::numeric_limits<double>::infinity();
	hdf_detail::Handle file_,root_;
};

} // namespace iga
#endif
