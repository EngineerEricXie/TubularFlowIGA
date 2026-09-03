#include "BlockCsr.hpp"
#include "BoundaryFlow.hpp"
#include "CaseInput.hpp"
#include "GenericCaseInput.hpp"
#include "CudaRuntime.hpp"
#include "DeviceMesh.hpp"
#include "Gmres.hpp"
#include "GenericTransportKernels.cuh"
#include "IgaCudaKernels.cuh"
#include "IgaDatabase.hpp"
#include "FlowCheckpoint.hpp"
#include "OutletModel.hpp"
#include "OutletCheckpoint.hpp"
#include "OutletFlow.hpp"
#include "PressureTraction.hpp"
#include "PressureTractionFlow.hpp"
#include "TemporalFunction.hpp"
#include "TransportCheckpoint.hpp"
#include "VelocitySeries.hpp"
#include "PhysiologyOutput.hpp"
#include "TemporalVtkHdf.hpp"
#include "VtkOutput.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

namespace iga::cuda {

struct BezierVtkHdfOutput {
	std::unique_ptr<iga::BezierVisualizationMesh> mesh;
	std::unique_ptr<iga::TemporalVtkHdfWriter> writer;

	void Initialize(iga::Database& database, const fs::path& output, bool resume)
	{
		mesh = std::make_unique<iga::BezierVisualizationMesh>(
			iga::BuildSourceCoordinateBezierVisualizationMesh(database, false));
		const auto report = iga::BezierGeometryReportPath(output);
		iga::WriteBezierGeometryReport(report, mesh->validation);
		iga::RequireValidBezierGeometry(mesh->validation);
		writer = std::make_unique<iga::TemporalVtkHdfWriter>(
			iga::VtkHdfPath(output), *mesh, resume);
		std::cout << "bezier_geometry_points=" << mesh->points.size()
			<< " local_point_references=" << mesh->validation.local_points
			<< " geometry_report=" << report.string()
			<< " vtkhdf=" << writer->path().string() << '\n';
	}
};

class ReferenceData {
public:
	ReferenceData()
		: basis_(kQuadraturePoints*64), gradient_(kQuadraturePoints*64*3),
		  hessian_(kQuadraturePoints*64*9), weight_(kQuadraturePoints),
		  surface_basis_(6*16*64), surface_gradient_(6*16*64*3), surface_weight_(6*16)
	{
		constexpr std::array<double,4> points{{0.06943184420297371, 0.33000947820757187,
			0.6699905217924281, 0.9305681557970262}};
		constexpr std::array<double,4> weights{{0.3478548451374539, 0.6521451548625461,
			0.6521451548625461, 0.3478548451374539}};
		std::vector<double> basis(kQuadraturePoints*64);
		std::vector<double> gradient(kQuadraturePoints*64*3);
		std::vector<double> hessian(kQuadraturePoints*64*9);
		std::vector<double> weight(kQuadraturePoints);
		std::vector<double> surface_basis(6*16*64);
		std::vector<double> surface_gradient(6*16*64*3);
		std::vector<double> surface_weight(6*16);
		for (int qz = 0; qz < 4; ++qz)
			for (int qy = 0; qy < 4; ++qy)
				for (int qx = 0; qx < 4; ++qx) {
					const int q = qz*16+qy*4+qx;
					const double coordinate[3] = {points[qx], points[qy], points[qz]};
					double b[3][4], db[3][4], d2b[3][4];
					for (int d = 0; d < 3; ++d) {
						const double x = coordinate[d];
						b[d][0] = std::pow(1.0-x,3);
						b[d][1] = 3.0*std::pow(1.0-x,2)*x;
						b[d][2] = 3.0*(1.0-x)*x*x;
						b[d][3] = x*x*x;
						db[d][0] = -3.0*std::pow(1.0-x,2);
						db[d][1] = 3.0-12.0*x+9.0*x*x;
						db[d][2] = 3.0*(2.0-3.0*x)*x;
						db[d][3] = 3.0*x*x;
						d2b[d][0] = 6.0*(1.0-x);
						d2b[d][1] = -12.0+18.0*x;
						d2b[d][2] = 6.0-18.0*x;
						d2b[d][3] = 6.0*x;
					}
					int p = 0;
					for (int k = 0; k < 4; ++k)
						for (int j = 0; j < 4; ++j)
							for (int i = 0; i < 4; ++i, ++p) {
								basis[q*64+p] = b[0][i]*b[1][j]*b[2][k];
								gradient[(q*64+p)*3] = db[0][i]*b[1][j]*b[2][k];
								gradient[(q*64+p)*3+1] = b[0][i]*db[1][j]*b[2][k];
								gradient[(q*64+p)*3+2] = b[0][i]*b[1][j]*db[2][k];
								double* second = hessian.data()+(q*64+p)*9;
								second[0] = d2b[0][i]*b[1][j]*b[2][k];
								second[1] = db[0][i]*db[1][j]*b[2][k];
								second[2] = db[0][i]*b[1][j]*db[2][k];
								second[3] = second[1];
								second[4] = b[0][i]*d2b[1][j]*b[2][k];
								second[5] = b[0][i]*db[1][j]*db[2][k];
								second[6] = second[2];
								second[7] = second[5];
								second[8] = b[0][i]*b[1][j]*d2b[2][k];
							}
					weight[q] = weights[qx]*weights[qy]*weights[qz];
				}
		constexpr int fixed_axis[6] = {2, 1, 0, 1, 0, 2};
		constexpr int varying_axes[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0, 2}, {1, 2}, {0, 1}};
		constexpr double fixed_value[6] = {0.0, 0.0, 1.0, 1.0, 0.0, 1.0};
		for (int face = 0; face < 6; ++face)
			for (int qj = 0; qj < 4; ++qj)
				for (int qi = 0; qi < 4; ++qi) {
					double coordinate[3]{};
					coordinate[fixed_axis[face]] = fixed_value[face];
					coordinate[varying_axes[face][0]] = points[qi];
					coordinate[varying_axes[face][1]] = points[qj];
					double b[3][4], db[3][4];
					for (int d = 0; d < 3; ++d) {
						const double x = coordinate[d];
						b[d][0] = std::pow(1.0-x,3);
						b[d][1] = 3.0*std::pow(1.0-x,2)*x;
						b[d][2] = 3.0*(1.0-x)*x*x;
						b[d][3] = x*x*x;
						db[d][0] = -3.0*std::pow(1.0-x,2);
						db[d][1] = 3.0-12.0*x+9.0*x*x;
						db[d][2] = 3.0*(2.0-3.0*x)*x;
						db[d][3] = 3.0*x*x;
					}
					const int q = face*16+qj*4+qi;
					int p = 0;
					for (int k = 0; k < 4; ++k)
						for (int j = 0; j < 4; ++j)
							for (int i = 0; i < 4; ++i, ++p) {
								surface_basis[q*64+p] = b[0][i]*b[1][j]*b[2][k];
								surface_gradient[(q*64+p)*3] = db[0][i]*b[1][j]*b[2][k];
								surface_gradient[(q*64+p)*3+1] = b[0][i]*db[1][j]*b[2][k];
								surface_gradient[(q*64+p)*3+2] = b[0][i]*b[1][j]*db[2][k];
							}
					surface_weight[q] = 0.25*weights[qi]*weights[qj];
				}
		basis_.CopyFromHost(basis.data(), basis.size());
		gradient_.CopyFromHost(gradient.data(), gradient.size());
		hessian_.CopyFromHost(hessian.data(), hessian.size());
		weight_.CopyFromHost(weight.data(), weight.size());
		surface_basis_.CopyFromHost(surface_basis.data(), surface_basis.size());
		surface_gradient_.CopyFromHost(surface_gradient.data(), surface_gradient.size());
		surface_weight_.CopyFromHost(surface_weight.data(), surface_weight.size());
		view_ = {basis_.data(), gradient_.data(), hessian_.data(), weight_.data(),
			surface_basis_.data(), surface_gradient_.data(), surface_weight_.data()};
	}

	const ReferenceView& view() const { return view_; }

private:
	DeviceBuffer<double> basis_, gradient_, hessian_, weight_;
	DeviceBuffer<double> surface_basis_, surface_gradient_, surface_weight_;
	ReferenceView view_;
};

class GeometryData {
public:
	GeometryData(const DeviceMesh& mesh, const ReferenceData& reference, bool with_hessian)
		: count_(static_cast<std::size_t>(mesh.view().elements)*kQuadraturePoints),
		  determinant_(count_), inverse_(count_*9),
		  inverse_second_(with_hessian ? count_*27 : 0)
	{
		const int threads = 128;
		const int blocks = static_cast<int>((count_+threads-1)/threads);
		BuildGeometryKernel<<<blocks,threads>>>(mesh.view(), reference.view(),
			determinant_.data(), inverse_.data(), inverse_second_.data(), with_hessian ? 1 : 0);
		CheckKernel("BuildGeometryKernel");
		std::vector<double> determinant(count_);
		determinant_.CopyToHost(determinant.data(), determinant.size());
		minimum_ = *std::min_element(determinant.begin(), determinant.end());
		bad_samples_ = static_cast<std::size_t>(std::count_if(determinant.begin(), determinant.end(),
			[](double value) { return !std::isfinite(value) || value <= 0.0; }));
		view_ = {determinant_.data(), inverse_.data(), inverse_second_.data()};
	}

	const GeometryView& view() const { return view_; }
	double minimum() const { return minimum_; }
	std::size_t bad_samples() const { return bad_samples_; }
	std::size_t bytes() const
	{
		return determinant_.bytes()+inverse_.bytes()+inverse_second_.bytes();
	}

private:
	std::size_t count_;
	DeviceBuffer<double> determinant_, inverse_, inverse_second_;
	GeometryView view_;
	double minimum_ = 0.0;
	std::size_t bad_samples_ = 0;
};

class ElementTiles {
public:
	explicit ElementTiles(const FlatMesh& mesh)
	{
		std::vector<int> element, first;
		for (std::size_t e = 0; e < mesh.elements(); ++e) {
			const int nen = mesh.element_offsets[e+1]-mesh.element_offsets[e];
			for (int pair = 0; pair < nen*nen; pair += kPairTile) {
				element.push_back(static_cast<int>(e));
				first.push_back(pair);
			}
		}
		element_.Allocate(element.size());
		first_.Allocate(first.size());
		element_.CopyFromHost(element.data(), element.size());
		first_.CopyFromHost(first.data(), first.size());
		view_ = {static_cast<int>(element.size()), element_.data(), first_.data()};
	}

	const ElementTilesView& view() const { return view_; }
	std::size_t bytes() const { return element_.bytes()+first_.bytes(); }

private:
	DeviceBuffer<int> element_, first_;
	ElementTilesView view_;
};

template <std::size_t Width>
DeviceBuffer<double> Flatten(const std::vector<std::array<double,Width>>& values)
{
	DeviceBuffer<double> result(values.size()*Width);
	result.CopyFromHost(reinterpret_cast<const double*>(values.data()), values.size()*Width);
	return result;
}

DeviceBuffer<int> CopyLabels(const std::vector<int>& labels)
{
	DeviceBuffer<int> result(labels.size());
	result.CopyFromHost(labels.data(), labels.size());
	return result;
}

DeviceBuffer<double> CopyScalars(const std::vector<double>& values)
{
	DeviceBuffer<double> result(values.size());
	result.CopyFromHost(values.data(), values.size());
	return result;
}

void RequireGeometry(const GeometryData& geometry)
{
	if (geometry.bad_samples())
		throw std::runtime_error("mesh has "+std::to_string(geometry.bad_samples())
			+" non-positive quadrature Jacobians; minimum detJ "+std::to_string(geometry.minimum()));
}

void PrintDevice()
{
	int device = 0;
	Check(cudaGetDevice(&device), "cudaGetDevice");
	cudaDeviceProp properties{};
	Check(cudaGetDeviceProperties(&properties, device), "cudaGetDeviceProperties");
	std::size_t free = 0, total = 0;
	Check(cudaMemGetInfo(&free, &total), "cudaMemGetInfo");
	std::cout << "device=" << device << " name=\"" << properties.name << "\""
		<< " compute_capability=" << properties.major << '.' << properties.minor
		<< " fp64_ratio=1/" << properties.singleToDoublePrecisionPerfRatio
		<< " memory_free_gib=" << Gibibytes(free)
		<< " memory_total_gib=" << Gibibytes(total) << '\n';
}

void PrintMesh(const FlatMesh& mesh, const BlockPattern* pattern = nullptr)
{
	std::cout << "mesh nodes=" << mesh.nodes << " elements=" << mesh.elements()
		<< " basis_entries=" << mesh.connectivity.size()
		<< " extraction_nonzeros=" << mesh.extraction_values.size()
		<< " max_basis=" << mesh.maximum_basis;
	if (pattern)
		std::cout << " block_nonzeros=" << pattern->columns.size()
			<< " mean_block_row=" << static_cast<double>(pattern->columns.size())/mesh.nodes;
	std::cout << '\n';
}

void CheckAvailableMemory(std::size_t required)
{
	std::size_t free = 0, total = 0;
	Check(cudaMemGetInfo(&free, &total), "cudaMemGetInfo");
	std::cout << "gpu_memory planned_gib=" << Gibibytes(required)
		<< " currently_free_gib=" << Gibibytes(free) << '\n';
	if (required > static_cast<std::size_t>(0.9*static_cast<double>(free)))
		throw std::runtime_error("planned GPU allocations exceed 90% of currently free memory");
}

void WriteTransport(const fs::path& path, const std::vector<double>& values)
{
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create "+path.string());
	output << std::setprecision(17);
	for (std::size_t node = 0; node < values.size()/2; ++node)
		output << node << ' ' << values[node*2] << ' ' << values[node*2+1] << '\n';
}

void WriteTransportVtk(const fs::path& mesh_path, const fs::path& path, const std::vector<double>& values)
{
	std::ifstream mesh(mesh_path);
	std::ofstream output(path);
	if (!mesh) throw std::runtime_error("cannot open VTK mesh "+mesh_path.string());
	if (!output) throw std::runtime_error("cannot create "+path.string());
	output << mesh.rdbuf() << "\nSCALARS N0 double 1\nLOOKUP_TABLE default\n";
	output << std::setprecision(17);
	for (std::size_t node = 0; node < values.size()/2; ++node) output << values[node*2] << '\n';
	output << "SCALARS Nplus double 1\nLOOKUP_TABLE default\n";
	for (std::size_t node = 0; node < values.size()/2; ++node) output << values[node*2+1] << '\n';
}

void WriteNavierStokesVtk(const fs::path& mesh_path, const fs::path& path, const std::vector<double>& values)
{
	std::ifstream mesh(mesh_path);
	std::ofstream output(path);
	if (!mesh) throw std::runtime_error("cannot open VTK mesh "+mesh_path.string());
	if (!output) throw std::runtime_error("cannot create "+path.string());
	output << mesh.rdbuf() << "\nVECTORS Velocity double\n";
	output << std::setprecision(17);
	for (std::size_t node = 0; node < values.size()/4; ++node)
		output << values[node*4] << ' ' << values[node*4+1] << ' ' << values[node*4+2] << '\n';
	output << "SCALARS Pressure double 1\nLOOKUP_TABLE default\n";
	for (std::size_t node = 0; node < values.size()/4; ++node) output << values[node*4+3] << '\n';
}

void WriteNavierStokesVtu(const fs::path& mesh_path, const fs::path& path,
	const std::vector<double>& values, double physical_time,
	iga::VisualizationFormat visualization_format,
	iga::TemporalVtkHdfWriter* vtkhdf)
{
	std::vector<double> velocity(3*values.size()/4), pressure(values.size()/4);
	for (std::size_t node = 0; node < values.size()/4; ++node) {
		for (int component = 0; component < 3; ++component)
			velocity[3*node+component] = values[4*node+component];
		pressure[node] = values[4*node+3];
	}
	std::vector<iga::VtkPointArray> arrays{
		{"velocity", 3, std::move(velocity)},
		{"pressure", 1, std::move(pressure)}};
	if (visualization_format == iga::VisualizationFormat::Vtu)
		iga::WriteVtu(mesh_path, path, arrays, physical_time);
	else {
		if (!vtkhdf) throw std::runtime_error("VTKHDF writer is unavailable");
		vtkhdf->Append(physical_time, arrays);
	}
}

void WriteNavierStokes(const fs::path& path, const std::vector<double>& values)
{
	std::ofstream velocity(path), pressure(path.string()+".pressure");
	if (!velocity || !pressure) throw std::runtime_error("cannot create Navier-Stokes output");
	velocity << std::setprecision(17);
	pressure << std::setprecision(17);
	for (std::size_t node = 0; node < values.size()/4; ++node) {
		velocity << values[node*4] << ' ' << values[node*4+1] << ' ' << values[node*4+2] << '\n';
		pressure << values[node*4+3] << '\n';
	}
}

struct CudaFlowOptions {
	fs::path database;
	fs::path case_dir;
	std::string system;
	int maximum_newton = 12;
	double nonlinear_relative_tolerance = 1e-5;
	double nonlinear_absolute_tolerance = 1e-10;
	double mass_relative_tolerance = 1e-3;
	fs::path output;
	fs::path checkpoint;
	fs::path restart;
	int output_every = 0;
	int checkpoint_every = 0;
	int stop_after_step = 0;
	iga::VisualizationFormat visualization_format = iga::VisualizationFormat::Automatic;
};

int ParsePositiveInteger(const std::string& text, const std::string& option)
{
	std::size_t used = 0;
	int value = 0;
	try {
		value = std::stoi(text, &used);
	} catch (const std::exception&) {
		throw std::runtime_error(option+" requires a positive integer");
	}
	if (used != text.size() || value <= 0)
		throw std::runtime_error(option+" requires a positive integer");
	return value;
}

double ParsePositiveFiniteDouble(const std::string& text, const std::string& option)
{
	std::size_t used = 0;
	double value = 0.0;
	try {
		value = std::stod(text, &used);
	} catch (const std::exception&) {
		throw std::runtime_error(option+" requires a finite positive value");
	}
	if (used != text.size() || !std::isfinite(value) || value <= 0.0)
		throw std::runtime_error(option+" requires a finite positive value");
	return value;
}

CudaFlowOptions ParseCudaFlowOptions(int argc, char** argv)
{
	if (argc < 4) throw std::runtime_error(
		"usage: iga_cuda navier-stokes DATABASE.ntiga CASE_DIR [MAX_NEWTON] [OUTPUT] "
		"[--system NAME] [--max-newton N] [--output PATH] [--output-every N] "
		"[--checkpoint PREFIX] [--checkpoint-every N] [--restart PREFIX] "
		"[--stop-after-step N] [--nonlinear-rtol R] [--nonlinear-atol A] "
		"[--mass-rtol R] [--visualization-format auto|vtu|vtkhdf]");
	CudaFlowOptions options;
	options.database = argv[2];
	options.case_dir = argv[3];
	int positional = 0;
	for (int i = 4; i < argc; ++i) {
		const std::string argument(argv[i]);
		if (argument.rfind("--", 0) != 0) {
			if (positional == 0)
				options.maximum_newton = ParsePositiveInteger(argument, "MAX_NEWTON");
			else if (positional == 1) options.output = argument;
			else throw std::runtime_error("too many positional arguments");
			++positional;
			continue;
		}
		if (i+1 >= argc) throw std::runtime_error(argument+" requires a value");
		const std::string value(argv[++i]);
		if (argument == "--system") options.system = value;
		else if (argument == "--max-newton")
			options.maximum_newton = ParsePositiveInteger(value, argument);
		else if (argument == "--output") options.output = value;
		else if (argument == "--output-every")
			options.output_every = ParsePositiveInteger(value, argument);
		else if (argument == "--checkpoint") options.checkpoint = value;
		else if (argument == "--checkpoint-every")
			options.checkpoint_every = ParsePositiveInteger(value, argument);
		else if (argument == "--restart") options.restart = value;
		else if (argument == "--stop-after-step")
			options.stop_after_step = ParsePositiveInteger(value, argument);
		else if (argument == "--nonlinear-rtol")
			options.nonlinear_relative_tolerance = ParsePositiveFiniteDouble(value, argument);
		else if (argument == "--nonlinear-atol")
			options.nonlinear_absolute_tolerance = ParsePositiveFiniteDouble(value, argument);
		else if (argument == "--mass-rtol")
			options.mass_relative_tolerance = ParsePositiveFiniteDouble(value, argument);
		else if (argument == "--visualization-format")
			options.visualization_format = iga::ParseVisualizationFormat(value);
		else throw std::runtime_error("unknown option: "+argument);
	}
	if (options.output_every > 0 && options.output.empty())
		throw std::runtime_error("--output-every requires --output or legacy OUTPUT");
	if (options.checkpoint_every > 0 && options.checkpoint.empty())
		throw std::runtime_error("--checkpoint-every requires --checkpoint");
	return options;
}

struct CudaFlowConvergenceMetrics {
	double continuity_l2 = 0.0;
	double continuity_sum = 0.0;
	double net_boundary_flow = 0.0;
	double absolute_boundary_flow = 0.0;
	double relative_mass_imbalance = 0.0;
};

std::vector<iga::Element> LoadBoundaryElements(iga::Database& database)
{
	std::vector<iga::Element> result;
	for (std::uint64_t index = 0; index < database.header().elements; ++index) {
		auto element = database.Load(index);
		if (std::any_of(element.boundary_labels.begin(), element.boundary_labels.end(),
			[](std::int32_t label) { return label >= 0; }))
			result.push_back(std::move(element));
	}
	return result;
}

CudaFlowConvergenceMetrics MeasureFlowConvergence(
	const std::vector<iga::Element>& boundary_elements,
	const std::vector<double>& state, const std::vector<double>& raw_rhs)
{
	if (state.size() != raw_rhs.size() || state.size()%4 != 0)
		throw std::runtime_error("CUDA flow convergence vectors have inconsistent sizes");
	CudaFlowConvergenceMetrics result;
	for (std::size_t node = 0; node < raw_rhs.size()/4; ++node) {
		const auto continuity = raw_rhs[4*node+3];
		result.continuity_l2 += continuity*continuity;
		result.continuity_sum += continuity;
	}
	result.continuity_l2 = std::sqrt(result.continuity_l2);
	std::map<int, double> boundary_flow;
	for (const auto& element : boundary_elements) {
		std::vector<std::array<double, 4>> nodal(element.connectivity.size());
		for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
			const auto node = element.connectivity[a];
			if (node < 0 || static_cast<std::size_t>(node) >= state.size()/4)
				throw std::runtime_error("boundary element references an invalid node");
			for (int field = 0; field < 4; ++field)
				nodal[a][field] = state[4*static_cast<std::size_t>(node)+field];
		}
		for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
			const auto label = element.boundary_labels[face];
			if (label >= 0)
				boundary_flow[label] += iga::IntegrateBoundaryFlow(element, face, nodal);
		}
	}
	for (const auto& entry : boundary_flow) {
		result.net_boundary_flow += entry.second;
		result.absolute_boundary_flow += std::abs(entry.second);
	}
	if (result.absolute_boundary_flow > 0.0)
		result.relative_mass_imbalance = 2.0*std::abs(result.net_boundary_flow)
			/result.absolute_boundary_flow;
	return result;
}

std::string PreciseNumber(double value)
{
	std::ostringstream stream;
	stream << std::scientific << std::setprecision(16) << value;
	return stream.str();
}

std::vector<double> ReadCudaFlowCheckpoint(const fs::path& prefix,
	const iga::FlowCheckpointMetadata& metadata)
{
	if (metadata.state_format != "raw_float64")
		throw std::runtime_error("CUDA flow restart requires raw_float64 checkpoint state");
	fs::path path(metadata.state_file);
	if (path.is_relative()) path = iga::FlowCheckpointMetadataPath(prefix).parent_path()/path;
	std::ifstream input(path, std::ios::binary);
	if (!input) throw std::runtime_error("cannot open CUDA flow checkpoint state: "+path.string());
	std::vector<double> values(static_cast<std::size_t>(metadata.nodes)*metadata.fields);
	input.read(reinterpret_cast<char*>(values.data()),
		static_cast<std::streamsize>(values.size()*sizeof(double)));
	char extra = 0;
	if (!input || input.read(&extra, 1))
		throw std::runtime_error("CUDA flow checkpoint state is truncated or has trailing data");
	return values;
}

void WriteCudaFlowCheckpoint(const fs::path& prefix,
	iga::FlowCheckpointMetadata metadata, const DeviceBuffer<double>& state,
	std::size_t values)
{
	std::vector<double> host(values);
	state.CopyToHost(host.data(), host.size());
	const auto state_path = iga::FlowCheckpointStatePath(prefix);
	std::ofstream output(state_path, std::ios::binary);
	if (!output) throw std::runtime_error("cannot create CUDA flow checkpoint state: "+state_path.string());
	output.write(reinterpret_cast<const char*>(host.data()),
		static_cast<std::streamsize>(host.size()*sizeof(double)));
	if (!output) throw std::runtime_error("cannot write CUDA flow checkpoint state: "+state_path.string());
	output.close();
	if (!output) throw std::runtime_error("cannot close CUDA flow checkpoint state: "+state_path.string());
	metadata.state_file = state_path.filename().string();
	metadata.state_format = "raw_float64";
	iga::WriteFlowCheckpointMetadata(prefix, metadata);
}

int MeshCheck(const std::string& database_path)
{
	iga::Database database(database_path);
	FlatMesh host(database);
	PrintMesh(host);
	DeviceMesh mesh(host);
	ReferenceData reference;
	GeometryData geometry(mesh, reference, false);
	std::cout << "geometry min_detJ=" << std::setprecision(17) << geometry.minimum()
		<< " bad_samples=" << geometry.bad_samples()
		<< " device_mesh_gib=" << Gibibytes(mesh.bytes())
		<< " geometry_gib=" << Gibibytes(geometry.bytes()) << '\n';
	return geometry.bad_samples() ? 1 : 0;
}

#include "ConfiguredTransportHost.inc"

int Transport(int argc, char** argv)
{
	if (argc < 4)
		throw std::runtime_error("usage: iga_cuda transport DATABASE.ntiga CASE_DIR [STEPS] [OUTPUT] [VELOCITY]");
	const auto total_start = Clock::now();
	iga::Database database(argv[2]);
	const fs::path case_dir(argv[3]);
	auto parameters = iga::ReadTransportParameters((case_dir/"simulation_parameter.txt").string());
	if (argc >= 5) parameters.steps = std::stoi(argv[4]);
	const auto labels = iga::ReadPointLabels((case_dir/"controlmesh.vtk").string(), database.header().nodes);
	const fs::path velocity_path = argc >= 7 ? fs::path(argv[6]) : case_dir/"initial_velocityfield.txt";
	const auto velocity_host = iga::ReadVelocity(velocity_path.string(), database.header().nodes);
	const auto case_configuration = iga::ReadCaseConfiguration((case_dir/"case_config.json").string());
	const auto boundaries = iga::ResolveBoundaryConditions(case_configuration, labels, velocity_host, parameters);
	std::cout << "boundary_config=" << (case_configuration.present ? "case_config.json" : "legacy-defaults")
		<< " transport_nodes=" << boundaries.transport_nodes << '\n';

	FlatMesh host(database);
	BlockPattern pattern_host(host);
	PrintMesh(host, &pattern_host);
	const auto preprocess_end = Clock::now();
	DeviceMesh mesh(host);
	DevicePattern pattern(pattern_host);
	ReferenceData reference;
	GeometryData geometry(mesh, reference, false);
	RequireGeometry(geometry);
	ElementTiles tiles(host);
	auto boundary_mask = CopyLabels(boundaries.transport_constrained);
	auto boundary_n0 = CopyScalars(boundaries.n0);
	auto boundary_nplus = CopyScalars(boundaries.nplus);
	auto velocity = Flatten(velocity_host);
	constexpr int linear_restart = 50;
	const std::size_t matrix_bytes = static_cast<std::size_t>(pattern.view().blocks)*4*sizeof(double);
	const std::size_t vector_bytes = static_cast<std::size_t>(host.nodes)*2*sizeof(double);
	CheckAvailableMemory(2*matrix_bytes+(11+linear_restart)*vector_bytes);
	BlockMatrix<2> left(pattern), previous(pattern);
	DeviceBuffer<double> current(host.nodes*2), next(host.nodes*2), rhs(host.nodes*2);
	left.Clear(); previous.Clear(); current.Clear(); next.Clear(); rhs.Clear();

	const auto assembly_start = Clock::now();
	AssembleTransportKernel<<<tiles.view().count,kPairTile>>>(
		mesh.view(), reference.view(), geometry.view(), tiles.view(), pattern.view(),
		velocity.data(), parameters, left.values(), previous.values());
	CheckKernel("AssembleTransportKernel");
	const int node_blocks = (static_cast<int>(host.nodes)+255)/256;
	ApplyTransportBoundaryKernel<<<node_blocks,256>>>(
		pattern.view(), boundary_mask.data(), left.values(), previous.values());
	CheckKernel("ApplyTransportBoundaryKernel");
	SetTransportBoundaryVectorKernel<<<node_blocks,256>>>(
		static_cast<int>(host.nodes), boundary_mask.data(), boundary_n0.data(),
		boundary_nplus.data(), current.data());
	CheckKernel("SetTransportBoundaryVectorKernel current");
	Check(cudaDeviceSynchronize(), "transport assembly synchronize");
	const auto assembly_end = Clock::now();

	long long total_iterations = 0;
	unsigned int singular_blocks = 0;
	double peak_gpu_used = 0.0;
	GmresWorkspace<2> linear_workspace(pattern.view().nodes, linear_restart);
	const auto solve_start = Clock::now();
	for (int step = 0; step < parameters.steps; ++step) {
		BlockSpmvKernel<2><<<node_blocks,256>>>(
			pattern.view(), previous.values(), current.data(), rhs.data());
		CheckKernel("transport previous MatMult");
		SetTransportBoundaryVectorKernel<<<node_blocks,256>>>(
			static_cast<int>(host.nodes), boundary_mask.data(), boundary_n0.data(),
			boundary_nplus.data(), rhs.data());
		CheckKernel("SetTransportBoundaryVectorKernel rhs");
		if (step == 0) next.Clear();
		else Check(cudaMemcpy(next.data(), current.data(), current.bytes(), cudaMemcpyDeviceToDevice),
			"transport warm start copy");
		auto result = SolveGmres(left, linear_workspace, rhs.data(), next.data(), 10000, 1e-8, false, 0, step == 0);
		if (!result.converged)
			throw std::runtime_error("transport GMRES did not converge at step "+std::to_string(step)
				+"; residual "+std::to_string(result.residual));
		total_iterations += result.iterations;
		peak_gpu_used = std::max(peak_gpu_used, result.device_used_gib);
		singular_blocks = std::max(singular_blocks, result.singular_diagonal_blocks);
		std::swap(current, next);
		if ((step+1)%25 == 0 || step+1 == parameters.steps)
			std::cout << "transport step=" << step+1 << '/' << parameters.steps
				<< " iterations=" << result.iterations << " residual=" << result.residual << '\n';
	}
	Check(cudaDeviceSynchronize(), "transport solve synchronize");
	const auto solve_end = Clock::now();
	BlasHandle blas;
	double norm = 0.0;
	Check(cublasDnrm2(blas, static_cast<int>(host.nodes*2), current.data(), 1, &norm),
		"transport final norm");
	if (argc >= 6) {
		std::vector<double> output(host.nodes*2);
		current.CopyToHost(output.data(), output.size());
		WriteTransport(argv[5], output);
		WriteTransportVtk(case_dir/"controlmesh.vtk", std::string(argv[5])+".vtk", output);
	}
	std::cout << "transport_cuda nodes=" << host.nodes << " elements=" << host.elements()
		<< " steps=" << parameters.steps
		<< " preprocess_s=" << std::chrono::duration<double>(preprocess_end-total_start).count()
		<< " assembly_s=" << std::chrono::duration<double>(assembly_end-assembly_start).count()
		<< " solve_s=" << std::chrono::duration<double>(solve_end-solve_start).count()
		<< " total_iterations=" << total_iterations << " final_l2=" << norm
		<< " singular_diagonal_blocks=" << singular_blocks << " gpu_used_gib=" << peak_gpu_used << '\n';
	return 0;
}

int NavierStokes(int argc, char** argv)
{
	const auto options = ParseCudaFlowOptions(argc, argv);
	const auto total_start = Clock::now();
	iga::Database database(options.database.string());
	BezierVtkHdfOutput vtkhdf;
	const auto& case_dir = options.case_dir;
	const int maximum_newton = options.maximum_newton;
	const auto labels = iga::ReadPointLabels((case_dir/"controlmesh.vtk").string(), database.header().nodes);
	const auto boundary_velocity_host = iga::ReadVelocity((case_dir/"initial_velocityfield.txt").string(), database.header().nodes);
	iga::ResolvedBoundaryConditions boundaries;
	iga::SimulationConfiguration simulation_configuration;
	std::vector<iga::OutletModelState> outlet_models;
	std::map<int, double> pressure_tractions;
	double viscosity = 0.1, density = 1.0, dt = 0.0;
	int physical_steps = 1;
	bool configured = false, transient = false;
	std::string boundary_config;
	if (fs::exists(case_dir/"simulation_config.json")) {
			simulation_configuration = iga::ReadSimulationConfiguration((case_dir/"simulation_config.json").string());
			iga::RequireFlowOnlyCoupling(simulation_configuration.coupling,
				"CUDA Navier-Stokes");
			const auto& flow = iga::FindNavierStokesSystem(simulation_configuration, options.system);
		configured = true;
		transient = flow.time_integration == "backward_euler";
		viscosity = flow.viscosity;
		density = flow.density;
		dt = transient ? simulation_configuration.time.dt : 0.0;
		physical_steps = transient ? simulation_configuration.time.steps : 1;
		outlet_models = iga::InitializeOutletModels(simulation_configuration, flow);
		const auto waveform_configuration = transient
			? iga::MaterializeBoundaryWaveforms(simulation_configuration, case_dir, 0.0)
			: simulation_configuration;
		const auto initial_configuration = iga::MaterializeOutletPressures(
			waveform_configuration, outlet_models);
		pressure_tractions = iga::ExtractPressureTractions(
			initial_configuration,
			iga::FindNavierStokesSystem(initial_configuration, options.system));
		boundaries = iga::ResolveFlowBoundaries(initial_configuration,
			iga::FindNavierStokesSystem(initial_configuration, options.system), labels,
			boundary_velocity_host);
		boundary_config = "simulation_config.json";
	} else {
		const auto parameters = iga::ReadTransportParameters((case_dir/"simulation_parameter.txt").string());
		const auto case_configuration = iga::ReadCaseConfiguration((case_dir/"case_config.json").string());
		boundaries = iga::ResolveBoundaryConditions(case_configuration, labels, boundary_velocity_host, parameters);
		boundary_config = case_configuration.present ? "case_config.json" : "legacy-defaults";
	}
	if (!transient)
		for (const auto& model : outlet_models)
			if (model.kind != iga::FieldBoundaryKind::Resistance)
				throw std::runtime_error("RC/RCR outlets require backward_euler flow");
	if (options.stop_after_step > physical_steps)
		throw std::runtime_error("--stop-after-step exceeds configured physical steps");
	const auto run_end_step = options.stop_after_step > 0
		? options.stop_after_step : physical_steps;
	const auto visualization_format = iga::ResolveVisualizationFormat(
		options.visualization_format, transient);
	std::cout << "boundary_config=" << boundary_config << " viscosity=" << viscosity
		<< " density=" << density
		<< " time_integration=" << (transient ? "backward_euler" : "steady")
		<< " dt=" << dt << " steps=" << physical_steps
		<< " run_end_step=" << run_end_step
		<< " velocity_nodes=" << boundaries.velocity_nodes << " pressure_nodes=" << boundaries.pressure_nodes << '\n';
	if (!transient && (!options.restart.empty() || !options.checkpoint.empty()
		|| options.output_every > 0 || options.checkpoint_every > 0))
		throw std::runtime_error("restart, checkpoint, and time-indexed output require transient flow");

	FlatMesh host(database);
	const auto boundary_elements = LoadBoundaryElements(database);
	const auto pressure_traction_elements = iga::LoadPressureTractionElements(
		database, pressure_tractions);
	std::vector<iga::Element> outlet_elements;
	if (!outlet_models.empty()) {
		std::unordered_map<int, std::size_t> face_counts;
		for (const auto& model : outlet_models) face_counts.emplace(model.label, 0);
		for (std::uint64_t index = 0; index < database.header().elements; ++index) {
			auto element = database.Load(index);
			bool retained = false;
			for (const auto label : element.boundary_labels) {
				const auto found = face_counts.find(label);
				if (found == face_counts.end()) continue;
				++found->second;
				retained = true;
			}
			if (retained) outlet_elements.push_back(std::move(element));
		}
		for (const auto& count : face_counts)
			if (count.second == 0)
				throw std::runtime_error("outlet model label "+std::to_string(count.first)
					+" has no boundary faces in the .ntiga database; repack with iga_pack");
	}
	BlockPattern pattern_host(host);
	PrintMesh(host, &pattern_host);
	const auto preprocess_end = Clock::now();
	DeviceMesh mesh(host);
	DevicePattern pattern(pattern_host);
	ReferenceData reference;
	GeometryData geometry(mesh, reference, true);
	RequireGeometry(geometry);
	ElementTiles tiles(host);
	auto velocity_mask = CopyLabels(boundaries.velocity_constrained);
	auto pressure_mask = CopyLabels(boundaries.pressure_constrained);
	auto boundary_velocity = Flatten(boundaries.velocity);
	auto boundary_pressure = CopyScalars(boundaries.pressure);
	constexpr int linear_restart = 200;
	const std::size_t matrix_bytes = static_cast<std::size_t>(pattern.view().blocks)*16*sizeof(double);
	const std::size_t vector_bytes = static_cast<std::size_t>(host.nodes)*4*sizeof(double);
	CheckAvailableMemory(matrix_bytes+(15+linear_restart)*vector_bytes);
	BlockMatrix<4> jacobian(pattern);
	DeviceBuffer<double> state(host.nodes*4), previous(host.nodes*4), update(host.nodes*4),
		rhs(host.nodes*4), pressure_traction_rhs(host.nodes*4);
	state.Clear(); previous.Clear(); update.Clear(); rhs.Clear(); pressure_traction_rhs.Clear();
	int start_step = 0;
	if (transient && !options.restart.empty()) {
		const auto metadata = iga::ReadFlowCheckpointMetadata(options.restart);
		iga::ValidateFlowCheckpoint(metadata, database.header().nodes, physical_steps,
			dt, density, viscosity);
		const auto restart_state = ReadCudaFlowCheckpoint(options.restart, metadata);
		iga::RestoreOutletCheckpoint(metadata, outlet_models);
		state.CopyFromHost(restart_state.data(), restart_state.size());
		previous.CopyFromHost(restart_state.data(), restart_state.size());
		start_step = metadata.completed_step;
		std::cout << "restart=" << options.restart.string()
			<< " completed_step=" << start_step
			<< " physical_time=" << metadata.physical_time << '\n';
	} else if (transient) {
		std::vector<double> initial_state(host.nodes*4, 0.0);
		for (std::size_t node = 0; node < host.nodes; ++node) {
			if (boundaries.velocity_constrained[node])
				for (int field = 0; field < 3; ++field)
					initial_state[node*4+field] = boundaries.velocity[node][field];
			if (boundaries.pressure_constrained[node])
				initial_state[node*4+3] = boundaries.pressure[node];
		}
		state.CopyFromHost(initial_state.data(), initial_state.size());
		previous.CopyFromHost(initial_state.data(), initial_state.size());
	}
	if (!options.output.empty()
		&& visualization_format == iga::VisualizationFormat::BezierVtkHdf)
		vtkhdf.Initialize(database, options.output, !options.restart.empty());
	const int node_blocks = (static_cast<int>(host.nodes)+255)/256;
	BlasHandle blas;
	long long total_iterations = 0;
	double total_assembly = 0.0, total_linear = 0.0;
	double peak_gpu_used = 0.0;
	GmresWorkspace<4> linear_workspace(pattern.view().nodes, linear_restart);
	std::vector<double> convergence_state(host.nodes*4);
	std::vector<double> raw_rhs(host.nodes*4);
	std::vector<std::pair<double, fs::path>> vtk_snapshots;
	std::vector<iga::VelocitySnapshot> velocity_snapshots_written;
	if (transient && options.output_every > 0) {
		std::vector<double> initial_output(host.nodes*4);
		state.CopyToHost(initial_output.data(), initial_output.size());
		const auto text_path = iga::TimeIndexedPath(options.output, start_step);
		const auto vtk_path = iga::VtuStepPath(options.output, start_step);
		WriteNavierStokes(text_path, initial_output);
		WriteNavierStokesVtu(case_dir/"controlmesh.vtk", vtk_path,
			initial_output, start_step*dt, visualization_format,
			vtkhdf.writer.get());
		if (visualization_format == iga::VisualizationFormat::Vtu)
			vtk_snapshots.push_back({start_step*dt, vtk_path});
		velocity_snapshots_written.push_back({start_step*dt, text_path});
	}
	for (int step = start_step; step < run_end_step; ++step) {
		if (step > start_step)
			Check(cudaMemcpy(previous.data(), state.data(), vector_bytes, cudaMemcpyDeviceToDevice),
				"Navier-Stokes previous-state copy");
		const double physical_time = transient ? (step+1)*dt : 0.0;
		iga::SimulationConfiguration step_configuration;
		if (configured)
			step_configuration = transient
				? iga::MaterializeBoundaryWaveforms(
					simulation_configuration, case_dir, physical_time)
				: simulation_configuration;
		std::vector<double> previous_capacitor_pressure(outlet_models.size());
		for (std::size_t i = 0; i < outlet_models.size(); ++i)
			previous_capacitor_pressure[i] = outlet_models[i].capacitor_pressure;
		bool outlet_converged = false;
		const int maximum_outlet_iterations = outlet_models.empty() ? 1 : 12;
		for (int coupling = 0; coupling < maximum_outlet_iterations; ++coupling) {
			if (configured) {
				const auto boundary_configuration = iga::MaterializeOutletPressures(
					step_configuration, outlet_models);
				pressure_tractions = iga::ExtractPressureTractions(
					boundary_configuration,
					iga::FindNavierStokesSystem(boundary_configuration, options.system));
				boundaries = iga::ResolveFlowBoundaries(boundary_configuration,
					iga::FindNavierStokesSystem(boundary_configuration, options.system), labels,
					boundary_velocity_host);
				boundary_velocity.CopyFromHost(
					reinterpret_cast<const double*>(boundaries.velocity.data()), host.nodes*3);
				boundary_pressure.CopyFromHost(boundaries.pressure.data(), host.nodes);
			}
			const auto pressure_traction_host = iga::IntegratePressureTractionForces(
				pressure_tractions, pressure_traction_elements, host.nodes);
			pressure_traction_rhs.CopyFromHost(
				pressure_traction_host.data(), pressure_traction_host.size());
			double initial_residual = -1.0;
			double last_residual = -1.0;
			CudaFlowConvergenceMetrics last_convergence;
			bool converged = false;
			for (int nonlinear = 0; nonlinear < maximum_newton; ++nonlinear) {
			const auto assembly_start = Clock::now();
			jacobian.Clear();
			rhs.Clear();
			AssembleNavierStokesKernel<<<tiles.view().count,kPairTile>>>(
			mesh.view(), reference.view(), geometry.view(), tiles.view(), pattern.view(),
			state.data(), previous.data(), density, viscosity, dt,
			jacobian.values(), rhs.data());
			CheckKernel("AssembleNavierStokesKernel");
			const double one = 1.0;
			Check(cublasDaxpy(blas, static_cast<int>(host.nodes*4), &one,
			pressure_traction_rhs.data(), 1, rhs.data(), 1),
			"Navier-Stokes pressure traction rhs");
			rhs.CopyToHost(raw_rhs.data(), raw_rhs.size());
			state.CopyToHost(convergence_state.data(), convergence_state.size());
			const auto convergence = MeasureFlowConvergence(
				boundary_elements, convergence_state, raw_rhs);
			last_convergence = convergence;
		ApplyNavierStokesBoundaryKernel<<<node_blocks,256>>>(
			pattern.view(), velocity_mask.data(), pressure_mask.data(), jacobian.values());
		CheckKernel("ApplyNavierStokesBoundaryKernel");
		SetNavierStokesBoundaryRhsKernel<<<node_blocks,256>>>(
			static_cast<int>(host.nodes), velocity_mask.data(), pressure_mask.data(),
			boundary_velocity.data(), boundary_pressure.data(), state.data(), rhs.data());
		CheckKernel("SetNavierStokesBoundaryRhsKernel");
		Check(cudaDeviceSynchronize(), "Navier-Stokes assembly synchronize");
		const auto assembly_end = Clock::now();
		total_assembly += std::chrono::duration<double>(assembly_end-assembly_start).count();
		double residual_norm = 0.0;
		Check(cublasDnrm2(blas, static_cast<int>(host.nodes*4), rhs.data(), 1, &residual_norm),
			"Navier-Stokes residual norm");
		last_residual = residual_norm;
		if (initial_residual < 0.0) initial_residual = residual_norm;
		const auto residual_scale = std::sqrt(static_cast<double>(host.nodes*4));
		const auto residual_rms = residual_norm/residual_scale;
		const double nonlinear_tolerance = std::max(
			options.nonlinear_absolute_tolerance*residual_scale,
			options.nonlinear_relative_tolerance*initial_residual);
		if (residual_norm <= nonlinear_tolerance
			&& convergence.relative_mass_imbalance <= options.mass_relative_tolerance) {
			std::cout << "step=" << step+1 << " time=" << physical_time
				<< " converged newton=" << nonlinear << " residual_l2=" << residual_norm
				<< " residual_rms=" << residual_rms
				<< " tolerance=" << nonlinear_tolerance
				<< " absolute_rms_tolerance=" << options.nonlinear_absolute_tolerance
				<< " continuity_l2=" << convergence.continuity_l2
				<< " continuity_sum=" << convergence.continuity_sum
				<< " net_boundary_flow=" << convergence.net_boundary_flow
				<< " relative_mass_imbalance=" << convergence.relative_mass_imbalance
				<< " mass_tolerance=" << options.mass_relative_tolerance << '\n';
			converged = true;
			break;
		}
		update.Clear();
		const auto linear_start = Clock::now();
		auto result = SolveGmres(jacobian, linear_workspace, rhs.data(), update.data(), 5000, 1e-8, false);
		Check(cudaDeviceSynchronize(), "Navier-Stokes linear synchronize");
		const auto linear_end = Clock::now();
		total_linear += std::chrono::duration<double>(linear_end-linear_start).count();
		if (!result.converged)
			throw std::runtime_error("Navier-Stokes GMRES failed at nonlinear iteration "
				+std::to_string(nonlinear)+" (residual="+PreciseNumber(result.residual)
				+", nonlinear_residual="+PreciseNumber(residual_norm)
				+", continuity_l2="+PreciseNumber(convergence.continuity_l2)
				+", net_boundary_flow="+PreciseNumber(convergence.net_boundary_flow)
				+", relative_mass_imbalance="
				+PreciseNumber(convergence.relative_mass_imbalance)+")");
		total_iterations += result.iterations;
		peak_gpu_used = std::max(peak_gpu_used, result.device_used_gib);
		double update_norm = 0.0;
		Check(cublasDnrm2(blas, static_cast<int>(host.nodes*4), update.data(), 1, &update_norm),
			"Navier-Stokes update norm");
		Check(cublasDaxpy(blas, static_cast<int>(host.nodes*4), &one, update.data(), 1,
			state.data(), 1), "Navier-Stokes state update");
		std::cout << "step=" << step+1 << " time=" << physical_time
			<< " newton=" << nonlinear << " residual_l2=" << residual_norm
			<< " residual_rms=" << residual_rms
			<< " continuity_l2=" << convergence.continuity_l2
			<< " continuity_sum=" << convergence.continuity_sum
			<< " net_boundary_flow=" << convergence.net_boundary_flow
			<< " relative_mass_imbalance=" << convergence.relative_mass_imbalance
			<< " update_l2=" << update_norm << " linear_iterations=" << result.iterations
			<< " linear_residual=" << result.residual
			<< " singular_diagonal_blocks=" << result.singular_diagonal_blocks
			<< " assembly_s=" << std::chrono::duration<double>(assembly_end-assembly_start).count()
			<< " linear_s=" << std::chrono::duration<double>(linear_end-linear_start).count() << '\n';
		}
			if (!converged)
				throw std::runtime_error("Navier-Stokes nonlinear solve reached MAX_NEWTON at physical step "
					+std::to_string(step+1)+" (nonlinear_residual="
					+PreciseNumber(last_residual)+", continuity_l2="
					+PreciseNumber(last_convergence.continuity_l2)+", net_boundary_flow="
					+PreciseNumber(last_convergence.net_boundary_flow)
					+", relative_mass_imbalance="
					+PreciseNumber(last_convergence.relative_mass_imbalance)+")");
			if (outlet_models.empty()) {
				outlet_converged = true;
				break;
			}
			std::vector<double> host_state(host.nodes*4);
			state.CopyToHost(host_state.data(), host_state.size());
			const auto flows = iga::IntegrateOutletModelFlows(
				outlet_models, outlet_elements, host_state);
			const auto evaluated = iga::EvaluateOutletCoupling(
				outlet_models, previous_capacitor_pressure, flows, dt);
			const auto outlet_tolerance = iga::OutletCouplingTolerance(evaluated);
			std::cout << "step=" << step+1 << " time=" << physical_time
				<< " outlet_iteration=" << coupling
				<< " pressure_change=" << evaluated.maximum_pressure_change
				<< " tolerance=" << outlet_tolerance << '\n';
			if (evaluated.maximum_pressure_change <= outlet_tolerance) {
				iga::CommitOutletCoupling(outlet_models, evaluated);
				for (std::size_t i = 0; i < outlet_models.size(); ++i) {
					std::cout << "outlet label=" << outlet_models[i].label
						<< " flow=" << outlet_models[i].flow
						<< " pressure=" << outlet_models[i].pressure
						<< " capacitor_pressure=" << outlet_models[i].capacitor_pressure << '\n';
				}
				outlet_converged = true;
				break;
			}
			iga::RelaxOutletCoupling(outlet_models, evaluated);
		}
		if (!outlet_converged)
			throw std::runtime_error("outlet fixed-point iteration did not converge at physical step "
				+std::to_string(step+1));
		const auto completed_step = step+1;
		if (options.output_every > 0 && completed_step%options.output_every == 0) {
			std::vector<double> values(host.nodes*4);
			state.CopyToHost(values.data(), values.size());
			const auto text_path = iga::TimeIndexedPath(options.output, completed_step);
			const auto vtk_path = iga::VtuStepPath(options.output, completed_step);
			WriteNavierStokes(text_path, values);
			WriteNavierStokesVtu(case_dir/"controlmesh.vtk", vtk_path, values,
				completed_step*dt, visualization_format, vtkhdf.writer.get());
			if (visualization_format == iga::VisualizationFormat::Vtu)
				vtk_snapshots.push_back({completed_step*dt, vtk_path});
			velocity_snapshots_written.push_back({completed_step*dt, text_path});
		}
		if (!options.checkpoint.empty()
			&& (completed_step == run_end_step
				|| (options.checkpoint_every > 0 && completed_step%options.checkpoint_every == 0))) {
			iga::FlowCheckpointMetadata metadata;
			metadata.nodes = database.header().nodes;
			metadata.completed_step = completed_step;
			metadata.physical_time = completed_step*dt;
			metadata.dt = dt;
			metadata.density = density;
			metadata.viscosity = viscosity;
			iga::AppendOutletCheckpoint(outlet_models, metadata);
			WriteCudaFlowCheckpoint(options.checkpoint, metadata, state, host.nodes*4);
			std::cout << "checkpoint=" << options.checkpoint.string()
				<< " completed_step=" << completed_step << '\n';
		}
	}
	double state_norm = 0.0;
	Check(cublasDnrm2(blas, static_cast<int>(host.nodes*4), state.data(), 1, &state_norm),
		"Navier-Stokes final norm");
	std::vector<double> output(host.nodes*4);
	state.CopyToHost(output.data(), output.size());
	double velocity_squared = 0.0, pressure_squared = 0.0;
	for (std::size_t node = 0; node < host.nodes; ++node) {
		for (int field = 0; field < 3; ++field)
			velocity_squared += output[node*4+field]*output[node*4+field];
		pressure_squared += output[node*4+3]*output[node*4+3];
	}
	if (!options.output.empty()) {
		WriteNavierStokes(options.output, output);
		if (visualization_format == iga::VisualizationFormat::Vtu)
			WriteNavierStokesVtk(
				case_dir/"controlmesh.vtk", options.output.string()+".vtk", output);
		const double final_time = transient ? run_end_step*dt : 0.0;
		WriteNavierStokesVtu(case_dir/"controlmesh.vtk",
			iga::VtuFinalPath(options.output), output, final_time,
			visualization_format, vtkhdf.writer.get());
		if (visualization_format == iga::VisualizationFormat::Vtu) {
			if (vtk_snapshots.empty())
				vtk_snapshots.push_back({final_time, iga::VtuFinalPath(options.output)});
			iga::WritePvd(iga::PvdPath(options.output), vtk_snapshots);
		}
		if (!velocity_snapshots_written.empty())
			iga::WriteVelocityManifest(iga::VelocityManifestPath(options.output),
				velocity_snapshots_written);
	}
	std::cout << "navier_stokes_cuda nodes=" << host.nodes << " elements=" << host.elements()
		<< " steps=" << physical_steps << " run_end_step=" << run_end_step
		<< " preprocess_s=" << std::chrono::duration<double>(preprocess_end-total_start).count()
		<< " assembly_s=" << total_assembly << " linear_s=" << total_linear
		<< " total_linear_iterations=" << total_iterations << " state_l2=" << state_norm
		<< " velocity_l2=" << std::sqrt(velocity_squared)
		<< " pressure_l2=" << std::sqrt(pressure_squared) << " gpu_used_gib=" << peak_gpu_used << '\n';
	return 0;
}

} // namespace iga::cuda

int main(int argc, char** argv)
{
	try {
		if (argc < 2)
			throw std::runtime_error("usage: iga_cuda device-info|mesh-check|solve|transport|navier-stokes ...");
		const std::string command(argv[1]);
		if (command == "device-info") {
			iga::cuda::PrintDevice();
			return 0;
		}
		iga::cuda::PrintDevice();
		if (command == "mesh-check") {
			if (argc != 3) throw std::runtime_error("usage: iga_cuda mesh-check DATABASE.ntiga");
			return iga::cuda::MeshCheck(argv[2]);
		}
		if (command == "transport") return iga::cuda::Transport(argc, argv);
		if (command == "solve") return iga::cuda::SolveConfigured(argc, argv);
		if (command == "navier-stokes") return iga::cuda::NavierStokes(argc, argv);
		throw std::runtime_error("unknown command: "+command);
	} catch (const std::exception& error) {
		std::cerr << "iga_cuda: " << error.what() << '\n';
		return 1;
	}
}
