#include "ParallelBezierVisualization.hpp"
#include "PetscBezierVisualization.hpp"
#define main SerialBezierRegressionMain
#include "test_bezier_visualization.cpp"
#undef main
#include <iostream>

namespace {
void RequireParallel(bool value,const char* message)
{
	if(!value)throw std::runtime_error(message);
}
double Control(std::int32_t node,int component)
{
	if(node<0||node>=112)throw std::out_of_range("missing control node");
	const double x=(node%7)/3.,y=((node/7)%4)/3.,z=(node/28)/3.;
	return component==0?x:component==1?y:component==2?z:x+2*y+3*z;
}
void Run(MPI_Comm comm,const fs::path& root)
{
	int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
	const std::vector<iga::VtkArraySchema> fields{{"velocity",3},{"scalar",1}};
	const std::vector<iga::Element> global{MakeElement(0,0,0,true),MakeElement(1,1,0,true)};
	iga::CollectiveLocalStage(comm,"serial oracle write",[&] { if(rank==0) { fs::create_directories(root);WriteDatabase(root/"oracle.ntiga",global,112); } });
	iga::BezierVisualizationMesh serial;std::vector<iga::VtkPointArray> serial_values;
	std::map<std::string,std::size_t> serial_indices;
	iga::CollectiveLocalStage(comm,"serial oracle build",[&] {
		iga::Database database((root/"oracle.ntiga").string());serial=iga::BuildBezierVisualizationMesh(database);
		std::vector<iga::VtkPointArray> controls{{"velocity",3,{}},{"scalar",1,{}}};
		for(int node=0;node<112;++node) { for(int component=0;component<3;++component)controls[0].values.push_back(Control(node,component));controls[1].values.push_back(Control(node,3)); }
		serial_values=iga::ExtractBezierPointArrays(serial,controls);
		for(std::size_t point=0;point<serial.points.size();++point) {
			iga::BezierPointSignature signature;
			for(auto entry=serial.signature_offsets[point];entry<serial.signature_offsets[point+1];++entry)
				signature.key.push_back({serial.signature_nodes[entry],iga::detail::QuantizeCoefficient(serial.signature_coefficients[entry],1e12)});
			serial_indices.emplace(iga::EncodeBezierPointSignature(signature),point);
		}
	});
	Vec state=nullptr,saved=nullptr;
	iga::RequireCollectivePetscSuccess(comm,"Bezier test state create",VecCreateMPI(comm,PETSC_DECIDE,448,&state));
	iga::CollectiveLocalStage(comm,"Bezier test owned state seed",[&] {
		PetscInt first=0,last=0;
		RequireParallel(!VecGetOwnershipRange(state,&first,&last),"cannot query test state range");
		for(PetscInt row=first;row<last;++row)RequireParallel(!VecSetValue(state,row,Control(row/4,row%4),INSERT_VALUES),"cannot seed test state");
	});
	iga::RequireCollectivePetscSuccess(comm,"Bezier test state assembly begin",VecAssemblyBegin(state));
	iga::RequireCollectivePetscSuccess(comm,"Bezier test state assembly end",VecAssemblyEnd(state));
	iga::RequireCollectivePetscSuccess(comm,"Bezier test state snapshot create",VecDuplicate(state,&saved));
	iga::RequireCollectivePetscSuccess(comm,"Bezier test state snapshot",VecCopy(state,saved));
	for(int layout=0;layout<3;++layout) {
		std::vector<iga::Element> local;
		for(auto element:global) {
			const auto owner=layout==0?0:layout==1?static_cast<int>(element.id%ranks):ranks-1-static_cast<int>(element.id%ranks);
			if(owner==rank) { element.owner=rank;local.push_back(element); }
		}
		std::reverse(local.begin(),local.end());
		const auto reference_piece=iga::BuildParallelBezierPartition(comm,local,2,fields,Control);
		const auto extracted=iga::BuildPetscBezierPartition(state,local,2,112,fields);
		const auto& piece=extracted.piece;
		iga::CollectiveLocalStage(comm,"PETSc Bezier subset comparison",[&] {
			std::set<std::int32_t> requested;
			for(const auto& element:local)requested.insert(element.connectivity.begin(),element.connectivity.end());
			RequireParallel(extracted.requested_nodes==requested.size()&&extracted.requested_rows==4*requested.size()
				&&extracted.global_rows==448,"PETSc subset count is incorrect");
			RequireParallel(piece.grid.points==reference_piece.grid.points&&piece.grid.connectivity==reference_piece.grid.connectivity
				&&piece.point_ids==reference_piece.point_ids&&piece.cell_ids==reference_piece.cell_ids,"PETSc Bezier geometry differs");
			for(std::size_t array=0;array<piece.point_arrays.size();++array)
				RequireParallel(piece.point_arrays[array].values==reference_piece.point_arrays[array].values,"PETSc selected field differs from local oracle");
			std::cout<<"petsc_bezier rank="<<rank<<" ranks="<<ranks<<" layout="<<layout<<" selected_rows="<<extracted.requested_rows<<" global_rows=448\n";
		});
		iga::CollectiveLocalStage(comm,"parallel serial Bezier comparison",[&] {
			RequireParallel(piece.cell_ids.size()==local.size(),"local cell count changed");
			for(std::size_t point=0;point<piece.point_ids.size();++point) {
				const auto id=static_cast<std::uint64_t>(piece.point_ids[point]);
				const auto key=iga::EncodeBezierPointSignature(iga::BuildBezierPointSignature(global.at(id/64),id%64));
				const auto reference=serial_indices.at(key);
				for(unsigned q=0;q<3;++q) {
					RequireParallel(piece.grid.points[3*point+q]==serial.points[reference][q],"parallel Bezier coordinate differs from serial");
					RequireParallel(piece.point_arrays[0].values[3*point+q]==serial_values[0].values[3*reference+q],"parallel Bezier velocity differs from serial");
				}
				RequireParallel(piece.point_arrays[1].values[point]==serial_values[1].values[reference],"parallel Bezier scalar differs from serial");
			}
		});
		iga::WriteParallelVtkSnapshot(comm,root/("layout"+std::to_string(layout)),piece,.25);
	}
	for(const auto& mode:std::vector<std::string>{"missing","duplicate","coordinate","callback","jacobian"}) {
		std::vector<iga::Element> local;
		if(rank==0) {
			local=global;
			if(mode=="missing")local.pop_back();
			if(mode=="duplicate")local[1]=local[0];
			if(mode=="coordinate")local[1].bezier_points[0][0]+=.01;
			if(mode=="jacobian") { local.resize(1);for(auto& point:local[0].bezier_points)point[0]=-point[0]; }
		}
		bool rejected=false;std::string message;
		try { iga::BuildParallelBezierPartition(comm,local,mode=="jacobian"?1:2,fields,[&](auto node,int component) { if(mode=="callback")throw std::runtime_error("missing local field");return Control(node,component); }); }
		catch(const std::exception& error) { rejected=true;message=error.what(); }
		iga::CollectiveLocalStage(comm,"parallel Bezier rejection",[&] {
			RequireParallel(rejected,"invalid parallel Bezier geometry accepted");
			const auto expected=mode=="missing"?"coverage":mode=="duplicate"?"duplicate global point occurrence":mode=="coordinate"?"inconsistent shared coordinates":mode=="callback"?"missing local field":"Jacobian";
			RequireParallel(message.find(expected)!=std::string::npos,"unexpected Bezier rejection stage");
		});
	}
	std::vector<iga::Element> owned;
	if(rank==0)owned=global;
	for(const auto& mode:std::vector<std::string>{"shape","late-cap"}) {
		iga::PointIdentityLimits limits;if(mode=="late-cap")limits.max_wire_bytes=8;
		bool rejected=false;std::string message;
		try { iga::BuildPetscBezierPartition(state,owned,2,mode=="shape"&&rank==ranks-1?113:112,fields,limits); }
		catch(const std::exception& error) { rejected=true;message=error.what(); }
		iga::CollectiveLocalStage(comm,"PETSc Bezier expected rejection",[&] {
			RequireParallel(rejected,"invalid PETSc Bezier build accepted");
			RequireParallel(message.find(mode=="shape"?"state size":"wire cap")!=std::string::npos,"unexpected PETSc Bezier rejection");
		});
		const auto retry=iga::BuildPetscBezierPartition(state,owned,2,112,fields);
		iga::CollectiveLocalStage(comm,"PETSc Bezier retry",[&] { RequireParallel(retry.piece.cell_ids.size()==owned.size(),"PETSc Bezier retry failed"); });
	}
	PetscBool equal=PETSC_FALSE;
	iga::RequireCollectivePetscSuccess(comm,"Bezier test state equality",VecEqual(state,saved,&equal));
	iga::CollectiveLocalStage(comm,"Bezier test immutable source",[&] { RequireParallel(equal==PETSC_TRUE,"Bezier extraction changed source Vec"); });
	iga::RequireCollectivePetscSuccess(comm,"Bezier test state cleanup",VecDestroy(&state));
	iga::RequireCollectivePetscSuccess(comm,"Bezier test snapshot cleanup",VecDestroy(&saved));
	if(rank==0)std::cout<<"petsc_bezier ranks="<<ranks<<" subset_parity=exact rejections=2 retries=passed source_unchanged=true\n";
	if(rank==0)std::cout<<"parallel_bezier ranks="<<ranks<<" layouts=3 serial_parity=exact rejections=5 passed\n";
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int rank=0,ranks=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try {
		RequireParallel(argc==2&&ranks==3,"require three ranks and new output directory");Run(PETSC_COMM_WORLD,fs::path(argv[1])/"world");
		MPI_Comm group=MPI_COMM_NULL;MPI_Comm_split(PETSC_COMM_WORLD,rank==0?0:1,rank,&group);
		Run(group,fs::path(argv[1])/(rank==0?"self":"pair"));MPI_Comm_free(&group);
	} catch(const std::exception& error) { std::cerr<<"rank "<<rank<<": "<<error.what()<<'\n';MPI_Abort(PETSC_COMM_WORLD,1);return 1; }
	PetscFinalize();return 0;
}
