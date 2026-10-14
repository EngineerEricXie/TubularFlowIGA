// Run a multidomain graph case whose domains are a native tetrahedral ALE
// fluid and a tetrahedral solid joined by one FSI coupling. The graph is
// validated with the surface layouts of the matching interface, then the
// strong Dirichlet--Neumann coordinator advances the pair step by step.
#include "MultidomainConfig.hpp"
#include "NativeTetFsiCase.hpp"
#include "SimulationGraph.hpp"

#include <petscksp.h>

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

namespace fs = std::filesystem;

struct Options
{
	fs::path graph_case,output_directory;
};

Options ParseOptions(int argc,char** argv)
{
	Options options;
	for(int i=1;i<argc;++i){
		const std::string argument=argv[i];
		if((argument=="--graph-case"||argument=="--output-dir")&&i+1<argc){
			(argument=="--graph-case"?options.graph_case:options.output_directory)=argv[++i];
		}
		// Anything else is a PETSc option, already consumed by PetscInitialize.
	}
	if(options.graph_case.empty()||options.output_directory.empty())
		throw std::invalid_argument(
			"usage: native_tet_fsi_graph --graph-case ROOT --output-dir DIR [PETSc options]");
	return options;
}

const iga::GraphDomainDefinition& DomainById(const iga::MultidomainConfiguration& configuration,
	const std::string& id)
{
	for(const auto& domain:configuration.domains)if(domain.id==id)return domain;
	throw std::runtime_error("graph FSI coupling names an undeclared domain");
}

int InterfaceLabel(const iga::GraphDomainDefinition& domain,const std::string& interface_id)
{
	for(const auto& surface:domain.surface_interfaces)
		if(surface.first==interface_id){
			if(surface.second.size()!=1)
				throw std::runtime_error("native FSI surface '"+domain.id+"."+interface_id
					+"' must cover exactly one boundary label");
			return *surface.second.begin();
		}
	throw std::runtime_error("graph FSI coupling names an undeclared surface interface");
}

} // namespace

int main(int argc,char** argv)
{
	if(PetscInitialize(&argc,&argv,nullptr,nullptr))return 1;
	int rank=0,ranks=1,status=0;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try{
		const auto options=ParseOptions(argc,argv);
		const auto root=iga::multidomain_detail::CanonicalGraphCaseRoot(options.graph_case);
		const auto configuration=iga::ReadMultidomainConfiguration(
			(root/"simulation_config.json").string());
		if(configuration.fsi_couplings.size()!=1||!configuration.graph.Edges().empty())
			throw std::invalid_argument("native_tet_fsi_graph runs one FSI coupling without"
				" pressure/flow couplings; use iga_multidomain_flow for scalar graphs");
		const auto& coupling=configuration.fsi_couplings.front();
		const auto& fluid_domain=DomainById(configuration,coupling.fluid_domain_id);
		const auto& solid_domain=DomainById(configuration,coupling.structure_domain_id);
		const int label=InterfaceLabel(fluid_domain,coupling.fluid_interface_id);
		if(InterfaceLabel(solid_domain,coupling.structure_interface_id)!=label)
			throw std::invalid_argument("native FSI fluid and solid surfaces must share one boundary label");
		iga::NativeTetFsiCase scenario;
		scenario.fluid=iga::ReadNativeTetFsiFluidCase(root/fluid_domain.case_directory/"fluid.json");
		scenario.solid=iga::ReadNativeTetFsiSolidCase(root/solid_domain.case_directory/"solid.json");
		scenario.interface_label=label;
		scenario.dt_s=configuration.time.dt_s;
		scenario.steps=configuration.time.steps;
		iga::NativeTetFsiNames names;
		names.fluid_domain=fluid_domain.id;names.fluid_interface=coupling.fluid_interface_id;
		names.solid_domain=solid_domain.id;names.solid_interface=coupling.structure_interface_id;
		names.edge=coupling.id;
		iga::NativeTetFsiSession session(scenario,iga::ReadNativeTetFsiSolverOptions(ranks),names);
		// The full graph: both domains with their surface catalog and layout,
		// joined by the FSI edge. Construction validates the pair contract.
		std::vector<iga::DomainNode> nodes;
		const std::array<const iga::GraphDomainDefinition*,2> pair{{&fluid_domain,&solid_domain}};
		for(const auto* domain:pair){
			const bool fluid=domain==&fluid_domain;
			const auto& surface=fluid?session.FluidSurface():session.SolidSurface();
			nodes.emplace_back(domain->id,domain->kind,domain->ports,domain->species_bindings,
				std::vector<iga::DistributedSurfaceInterface>{surface},
				std::map<iga::SurfaceInterfaceRef,iga::DistributedSurfaceLayout>{{surface.id,session.Layout()}});
		}
		const iga::SimulationGraph graph(std::move(nodes),{},{},{session.Edge()});
		const auto& output=options.output_directory;
		for(int step=1;step<=scenario.steps;++step)session.RunStep(step,output,rank,ranks);
		session.WriteSeries(output);
		if(rank==0){
			std::ofstream manifest(output/"graph_binding_manifest.json");
			manifest<<"{\n  \"schema_version\": 1,\n  \"graph_case\": \""<<root.generic_string()
				<<"\",\n  \"fsi_coupling\": \""<<coupling.id
				<<"\",\n  \"fluid_domain\": \""<<fluid_domain.id
				<<"\",\n  \"fluid_kind\": \""<<iga::DomainKindName(fluid_domain.kind)
				<<"\",\n  \"structure_domain\": \""<<solid_domain.id
				<<"\",\n  \"structure_kind\": \""<<iga::DomainKindName(solid_domain.kind)
				<<"\",\n  \"interface_label\": "<<label
				<<",\n  \"interface_identity_sha256\": \""<<session.InterfaceIdentity()
				<<"\",\n  \"layout_identity_sha256\": \""<<session.Layout().layout_identity_sha256
				<<"\",\n  \"fsi_edge_identity_sha256\": \""<<iga::BuildFsiCouplingEdgeIdentitySha256(session.Edge())
				<<"\",\n  \"mpi_ranks\": "<<ranks<<",\n  \"accepted_steps\": "<<scenario.steps<<"\n}\n";
			if(!manifest)throw std::runtime_error("cannot write graph binding manifest");
			std::cout<<"native_fsi_graph_complete domains="<<graph.Domains().size()
				<<" fsi_edges="<<graph.FsiEdges().size()<<" accepted_steps="<<scenario.steps<<'\n';
		}
	}catch(const std::exception& error){
		if(rank==0)std::cerr<<"native_tet_fsi_graph: "<<error.what()<<'\n';
		status=1;
	}
	PetscFinalize();return status;
}
