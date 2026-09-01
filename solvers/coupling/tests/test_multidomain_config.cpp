#include "MultidomainConfig.hpp"
#include "BifurcationMultidomainCase.hpp"
#include "SequentialMultidomainCase.hpp"

#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void RequireRejected(const std::function<void()>& operation)
{
	bool rejected = false;
	try { operation(); }
	catch (const std::exception&) { rejected = true; }
	assert(rejected);
}

std::string ValidConfiguration()
{
	return R"json({
  "schema_version": 5,
  "time": {"dt": 0.01, "steps": 8},
  "start_domain": "upstream",
  "execution": {
    "kind": "aitken",
    "maximum_iterations": 12,
    "pressure_relative_tolerance": 1e-5,
    "pressure_reference_pa": 100,
    "flow_relative_tolerance": 1e-8,
    "relaxation_factor": 0.5,
    "minimum_relaxation": 0.1,
    "maximum_relaxation": 0.9
  },
  "domains": [
    {
      "id": "upstream", "dimension": "1d", "kind": "network_flow",
      "case": "domains/upstream", "inlet_policy": "configured_open_loop",
      "ports": [
        {"id":"terminal","locator_kind":"runtime_port","locator":"outlet:2",
         "provides":["area","flow_rate","mean_pressure"],"requires":["mean_pressure"]},
        {"id":"root_state","locator_kind":"runtime_port","locator":"root",
         "provides":["area","flow_rate","mean_pressure"],"requires":[]}
      ]
    },
    {
      "id": "roi", "dimension": "3d", "kind": "body_fitted_iga_flow",
      "case": "domains/roi", "database": "domains/roi/roi.ntiga",
      "ports": [
        {"id":"inlet","locator_kind":"boundary_label","locator":"1",
         "native_to_outward_sign":-1,
         "provides":["area","flow_rate","mean_pressure"],"requires":["flow_rate"]},
        {"id":"outlet","locator_kind":"boundary_label","locator":"2",
         "provides":["area","flow_rate","mean_pressure"],"requires":["mean_pressure"]},
        {"id":"wall","locator_kind":"boundary_label","locator":"0",
         "provides":["flow_rate"],"requires":[]}
      ]
    },
    {
      "id": "downstream", "dimension": "1d", "kind": "one_d_flow",
      "case": "domains/downstream", "inlet_policy": "coupled_root",
      "ports": [
        {"id":"root","locator_kind":"runtime_port","locator":"root",
         "provides":["area","flow_rate","mean_pressure"],"requires":["flow_rate"]},
        {"id":"terminal_state","locator_kind":"runtime_port","locator":"outlet:2",
         "provides":["area","flow_rate","mean_pressure"],"requires":[]}
      ]
    }
  ],
  "couplings": [
    {"id":"upstream_to_roi","a":{"domain":"upstream","port":"terminal"},
     "b":{"domain":"roi","port":"inlet"},"mode":"pressure_flow",
     "initial_pressure_pa":10},
    {"id":"roi_to_downstream","a":{"domain":"roi","port":"outlet"},
     "b":{"domain":"downstream","port":"root"},"mode":"pressure_flow",
     "initial_pressure_pa":5}
  ]
})json";
}

std::string Replace(std::string text, const std::string& from, const std::string& to)
{
	const auto position = text.find(from);
	if (position == std::string::npos) throw std::runtime_error("test replacement source missing");
	text.replace(position, from.size(), to);
	return text;
}

std::string InsertBeforeLast(std::string text, const std::string& marker,
	const std::string& insertion)
{
	const auto position = text.rfind(marker);
	if (position == std::string::npos) throw std::runtime_error("test insertion marker missing");
	text.insert(position, insertion);
	return text;
}

} // namespace

int main()
{
	const auto configuration = iga::ParseMultidomainConfiguration(ValidConfiguration());
	assert(configuration.schema_version == 5);
	assert(configuration.time.dt_s == 0.01 && configuration.time.steps == 8);
	assert(configuration.start_domain_id == "upstream");
	assert(configuration.execution.kind == iga::GraphExecutionKind::Aitken);
	assert(configuration.execution.maximum_iterations == 12);
	assert(configuration.domains.size() == 3);
	assert(configuration.domains[0].one_d_inlet_policy
		== iga::OneDInletPolicy::ConfiguredOpenLoop);
	assert(configuration.domains[1].database == "domains/roi/roi.ntiga");
	assert(configuration.domains[1].ports[0].orientation.native_to_outward_sign == -1);
	assert(configuration.domains[2].one_d_inlet_policy
		== iga::OneDInletPolicy::CoupledRoot);
	assert(configuration.graph.Domains().size() == 3);
	assert(configuration.graph.Edges().size() == 2);
	assert(configuration.initial_pressure_pa.at("upstream_to_roi") == 10.0);
	const auto plan = iga::MakeSequentialPlan(configuration.graph,
		configuration.start_domain_id);
	assert(plan.domain_ids
		== std::vector<std::string>({"upstream", "roi", "downstream"}));
	assert(plan.edge_ids
		== std::vector<std::string>({"upstream_to_roi", "roi_to_downstream"}));
	assert(iga::MakeSequentialPressureFlowPlan(configuration.graph,
		configuration.start_domain_id).domain_ids == plan.domain_ids);
	assert(configuration.graph.Port({"roi", "wall"}).requires.empty());
	const auto resolved = iga::ResolveSequentialOneDThreeDOneD(configuration);
	assert(resolved.upstream_domain_id == "upstream");
	assert(resolved.three_d_domain_id == "roi");
	assert(resolved.downstream_domain_id == "downstream");
	assert(resolved.upstream_terminal == iga::PortRef({"upstream", "terminal"}));
	assert(resolved.downstream_terminal_observation
		== iga::PortRef({"downstream", "terminal_state"}));
	const auto controls = iga::PressureFlowControlsFor(configuration.execution);
	assert(controls.method == iga::PressureFlowIterationMethod::Aitken);
	assert(controls.maximum_iterations == 12);
	RequireRejected([&configuration] {
		(void)iga::ResolveOneDThreeDBifurcation(configuration);
	});
	{
		const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
		const auto root = std::filesystem::temp_directory_path()
			/("tubularflowiga_graph_case_"+std::to_string(nonce));
		const auto outside = std::filesystem::temp_directory_path()
			/("tubularflowiga_graph_outside_"+std::to_string(nonce));
		std::filesystem::create_directories(root/"domains/upstream");
		std::filesystem::create_directories(root/"domains/roi");
		std::filesystem::create_directories(root/"domains/downstream");
		std::filesystem::create_directories(outside);
		std::ofstream(root/"domains/roi/roi.ntiga") << "fixture";
		const auto assets = iga::ResolveGraphDomainAssets(configuration, root);
		assert(assets.at("upstream").case_directory
			== std::filesystem::canonical(root/"domains/upstream"));
		assert(assets.at("roi").database
			== std::filesystem::canonical(root/"domains/roi/roi.ntiga"));
		std::filesystem::create_directory_symlink(outside, root/"escape");
		const auto escaped = iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"case\": \"domains/upstream\"", "\"case\": \"escape\""));
		RequireRejected([&escaped, &root] {
			(void)iga::ResolveGraphDomainAssets(escaped, root);
		});
		std::ofstream(outside/"network.swc") << "fixture";
		std::filesystem::create_symlink(outside/"network.swc",
			root/"domains/upstream/network.swc");
		RequireRejected([&root] {
			(void)iga::ResolveContainedCaseFile(
				std::filesystem::canonical(root/"domains/upstream"), "network.swc",
				"upstream geometry");
		});
		std::filesystem::remove_all(root);
		std::filesystem::remove_all(outside);
	}

	RequireRejected([] {
		iga::ParseMultidomainConfiguration(
			Replace(ValidConfiguration(), "\"schema_version\": 5", "\"schema_version\": 4"));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(
			Replace(ValidConfiguration(), "\"steps\": 8", "\"steps\": 0"));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"case\": \"domains/upstream\"", "\"case\": \"../upstream\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"database\": \"domains/roi/roi.ntiga\"", "\"database\": \"/tmp/roi.ntiga\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"dimension\": \"3d\", \"kind\": \"body_fitted_iga_flow\"",
			"\"dimension\": \"1d\", \"kind\": \"body_fitted_iga_flow\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"start_domain\": \"upstream\"", "\"start_domain\": \"missing\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"kind\": \"aitken\"", "\"kind\": \"invented\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"kind\": \"aitken\"", "\"kind\": \"explicit\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"maximum_iterations\": 12", "\"maximum_iterations\": 0"));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"minimum_relaxation\": 0.1", "\"minimum_relaxation\": 0.6"));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"provides\"", "\"unknown_key\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"native_to_outward_sign\":-1", "\"native_to_outward_sign\":0"));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"mode\":\"pressure_flow\"", "\"mode\":\"species\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"port\":\"inlet\"", "\"port\":\"missing\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"id\": \"roi\", \"dimension\": \"3d\"",
			"\"id\": \"upstream\", \"dimension\": \"3d\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"id\":\"roi_to_downstream\"", "\"id\":\"upstream_to_roi\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"flow_rate\"]",
			"\"provides\":[\"area\",\"mean_pressure\"],\"requires\":[\"flow_rate\"]"));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"database\": \"domains/roi/roi.ntiga\",", ""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"case\": \"domains/upstream\", \"inlet_policy\"",
			"\"case\": \"domains/upstream\", \"database\": \"upstream.ntiga\", \"inlet_policy\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"locator_kind\":\"runtime_port\"", "\"locator_kind\":\"boundary_label\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"locator\":\"outlet:2\",",
			"\"locator\":\"outlet:2\",\"native_to_outward_sign\":-1,"));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"locator\":\"outlet:2\"", "\"locator\":\"outlet:not-a-node\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"locator\":\"outlet:2\"", "\"locator\":\"outlet:0\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"requires\":[\"mean_pressure\"]", "\"requires\":[\"flow_rate\"]"));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"locator\":\"root\",\n         \"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[]",
			"\"locator\":\"root\",\n         \"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"flow_rate\"]"));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"locator\":\"root\",\n         \"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"flow_rate\"]",
			"\"locator\":\"root\",\n         \"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[]"));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"locator_kind\":\"boundary_label\",\"locator\":\"1\"",
			"\"locator_kind\":\"runtime_port\",\"locator\":\"1\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"locator_kind\":\"boundary_label\",\"locator\":\"1\"",
			"\"locator_kind\":\"boundary_label\",\"locator\":\"inlet\""));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			R"json("provides":["flow_rate"],"requires":[])json",
			R"json("provides":["flow_rate"],"requires":["mean_normal_traction"])json"));
	});
	RequireRejected([] {
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			R"json("provides":["flow_rate"],"requires":[])json",
			R"json("provides":["flow_rate"],"requires":["mean_pressure","mean_normal_traction"])json"));
	});
	RequireRejected([] {
		const std::string second = R"json(,
    {"id":"roi_to_downstream","a":{"domain":"roi","port":"outlet"},
     "b":{"domain":"downstream","port":"root"},"mode":"pressure_flow",
     "initial_pressure_pa":5})json";
		iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(), second, ""));
	});
	{
		const auto reverse = iga::ParseMultidomainConfiguration(Replace(ValidConfiguration(),
			"\"start_domain\": \"upstream\"", "\"start_domain\": \"downstream\""));
		RequireRejected([&reverse] {
			(void)iga::MakeSequentialPressureFlowPlan(reverse.graph,
				reverse.start_domain_id);
		});
		RequireRejected([&reverse] {
			(void)iga::ResolveSequentialOneDThreeDOneD(reverse);
		});
	}
	{
		auto branch_text = Replace(ValidConfiguration(),
			R"json({"id":"wall","locator_kind":"boundary_label","locator":"0",
         "provides":["flow_rate"],"requires":[]})json",
			R"json({"id":"wall","locator_kind":"boundary_label","locator":"0",
         "provides":["flow_rate"],"requires":[]},
        {"id":"branch_outlet","locator_kind":"boundary_label","locator":"3",
         "provides":["area","flow_rate","mean_pressure"],"requires":["mean_pressure"]})json");
		branch_text = Replace(branch_text, R"json(    }
  ],
  "couplings")json", R"json(    },
    {
      "id": "branch", "dimension": "1d", "kind": "network_flow",
      "case": "domains/branch", "inlet_policy": "coupled_root",
      "ports": [
        {"id":"root","locator_kind":"runtime_port","locator":"root",
         "provides":["area","flow_rate","mean_pressure"],"requires":["flow_rate"]},
        {"id":"terminal_state","locator_kind":"runtime_port","locator":"outlet:2",
         "provides":["area","flow_rate","mean_pressure"],"requires":[]}
      ]
    }
  ],
  "couplings")json");
		branch_text = InsertBeforeLast(branch_text, "\n  ]\n}", R"json(,
    {"id":"roi_to_branch","a":{"domain":"roi","port":"branch_outlet"},
     "b":{"domain":"branch","port":"root"},"mode":"pressure_flow",
     "initial_pressure_pa":4})json");
		const auto branch = iga::ParseMultidomainConfiguration(branch_text);
		assert(branch.graph.Domains().size() == 4 && branch.graph.Edges().size() == 3);
		RequireRejected([&branch] {
			(void)iga::MakeSequentialPressureFlowPlan(branch.graph, branch.start_domain_id);
		});
		const auto branch_plan = iga::MakeAcyclicPressureFlowPlan(branch.graph,
			branch.start_domain_id);
		assert(branch_plan.interfaces.size() == 3);
		const auto resolved_branch = iga::ResolveOneDThreeDBifurcation(branch);
		assert(resolved_branch.upstream_domain_id == "upstream");
		assert(resolved_branch.three_d_domain_id == "roi");
		assert(resolved_branch.upstream_edge_id == "upstream_to_roi");
		assert(resolved_branch.branches.size() == 2);
		assert(resolved_branch.branches[0].edge_id == "roi_to_branch");
		assert(resolved_branch.branches[1].edge_id == "roi_to_downstream");
		RequireRejected([&branch_text] {
			(void)iga::ParseMultidomainConfiguration(Replace(branch_text,
				"\"case\": \"domains/branch\", \"inlet_policy\": \"coupled_root\"",
				"\"case\": \"domains/branch\", \"inlet_policy\": \"configured_open_loop\""));
		});
	}
	std::cout << "multidomain configuration tests passed\n";
	return 0;
}
