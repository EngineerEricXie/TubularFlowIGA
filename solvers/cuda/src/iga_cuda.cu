#include "BlockCsr.hpp"
#include "CaseInput.hpp"
#include "CudaRuntime.hpp"
#include "DeviceMesh.hpp"
#include "Gmres.hpp"
#include "IgaCudaKernels.cuh"
#include "IgaDatabase.hpp"

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
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

namespace iga::cuda {

class ReferenceData {
public:
	ReferenceData()
		: basis_(kQuadraturePoints*64), gradient_(kQuadraturePoints*64*3),
		  hessian_(kQuadraturePoints*64*9), weight_(kQuadraturePoints)
	{
		constexpr std::array<double,4> points{{0.06943184420297371, 0.33000947820757187,
			0.6699905217924281, 0.9305681557970262}};
		constexpr std::array<double,4> weights{{0.3478548451374539, 0.6521451548625461,
			0.6521451548625461, 0.3478548451374539}};
		std::vector<double> basis(kQuadraturePoints*64);
		std::vector<double> gradient(kQuadraturePoints*64*3);
		std::vector<double> hessian(kQuadraturePoints*64*9);
		std::vector<double> weight(kQuadraturePoints);
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
		basis_.CopyFromHost(basis.data(), basis.size());
		gradient_.CopyFromHost(gradient.data(), gradient.size());
		hessian_.CopyFromHost(hessian.data(), hessian.size());
		weight_.CopyFromHost(weight.data(), weight.size());
		view_ = {basis_.data(), gradient_.data(), hessian_.data(), weight_.data()};
	}

	const ReferenceView& view() const { return view_; }

private:
	DeviceBuffer<double> basis_, gradient_, hessian_, weight_;
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
	auto labels_device = CopyLabels(labels);
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
		pattern.view(), labels_device.data(), left.values(), previous.values());
	CheckKernel("ApplyTransportBoundaryKernel");
	SetTransportBoundaryVectorKernel<<<node_blocks,256>>>(
		static_cast<int>(host.nodes), labels_device.data(), parameters.n0_bc,
		parameters.nplus_bc, current.data());
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
			static_cast<int>(host.nodes), labels_device.data(), parameters.n0_bc,
			parameters.nplus_bc, rhs.data());
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
	if (argc < 4)
		throw std::runtime_error("usage: iga_cuda navier-stokes DATABASE.ntiga CASE_DIR [MAX_NEWTON] [OUTPUT]");
	const auto total_start = Clock::now();
	iga::Database database(argv[2]);
	const fs::path case_dir(argv[3]);
	const int maximum_newton = argc >= 5 ? std::stoi(argv[4]) : 8;
	const auto parameters = iga::ReadTransportParameters((case_dir/"simulation_parameter.txt").string());
	const auto labels = iga::ReadPointLabels((case_dir/"controlmesh.vtk").string(), database.header().nodes);
	const auto boundary_velocity_host = iga::ReadVelocity(
		(case_dir/"initial_velocityfield.txt").string(), database.header().nodes);

	FlatMesh host(database);
	BlockPattern pattern_host(host);
	PrintMesh(host, &pattern_host);
	const auto preprocess_end = Clock::now();
	DeviceMesh mesh(host);
	DevicePattern pattern(pattern_host);
	ReferenceData reference;
	GeometryData geometry(mesh, reference, true);
	RequireGeometry(geometry);
	ElementTiles tiles(host);
	auto labels_device = CopyLabels(labels);
	auto boundary_velocity = Flatten(boundary_velocity_host);
	constexpr int linear_restart = 200;
	const std::size_t matrix_bytes = static_cast<std::size_t>(pattern.view().blocks)*16*sizeof(double);
	const std::size_t vector_bytes = static_cast<std::size_t>(host.nodes)*4*sizeof(double);
	CheckAvailableMemory(matrix_bytes+(13+linear_restart)*vector_bytes);
	BlockMatrix<4> jacobian(pattern);
	DeviceBuffer<double> state(host.nodes*4), update(host.nodes*4), rhs(host.nodes*4);
	state.Clear(); update.Clear(); rhs.Clear();
	const int node_blocks = (static_cast<int>(host.nodes)+255)/256;
	BlasHandle blas;
	double initial_residual = -1.0;
	long long total_iterations = 0;
	bool converged = false;
	double total_assembly = 0.0, total_linear = 0.0;
	double peak_gpu_used = 0.0;
	GmresWorkspace<4> linear_workspace(pattern.view().nodes, linear_restart);
	for (int nonlinear = 0; nonlinear < maximum_newton; ++nonlinear) {
		const auto assembly_start = Clock::now();
		jacobian.Clear();
		rhs.Clear();
		AssembleNavierStokesKernel<<<tiles.view().count,kPairTile>>>(
			mesh.view(), reference.view(), geometry.view(), tiles.view(), pattern.view(),
			state.data(), 0.1, jacobian.values(), rhs.data());
		CheckKernel("AssembleNavierStokesKernel");
		ApplyNavierStokesBoundaryKernel<<<node_blocks,256>>>(
			pattern.view(), labels_device.data(), jacobian.values());
		CheckKernel("ApplyNavierStokesBoundaryKernel");
		SetNavierStokesBoundaryRhsKernel<<<node_blocks,256>>>(
			static_cast<int>(host.nodes), labels_device.data(), boundary_velocity.data(),
			parameters.vplus, state.data(), rhs.data());
		CheckKernel("SetNavierStokesBoundaryRhsKernel");
		Check(cudaDeviceSynchronize(), "Navier-Stokes assembly synchronize");
		const auto assembly_end = Clock::now();
		total_assembly += std::chrono::duration<double>(assembly_end-assembly_start).count();
		double residual_norm = 0.0;
		Check(cublasDnrm2(blas, static_cast<int>(host.nodes*4), rhs.data(), 1, &residual_norm),
			"Navier-Stokes residual norm");
		if (initial_residual < 0.0) initial_residual = residual_norm;
		const double nonlinear_tolerance = std::max(1e-10, 1e-5*initial_residual);
		if (residual_norm <= nonlinear_tolerance) {
			std::cout << "converged newton=" << nonlinear << " residual_l2=" << residual_norm
				<< " tolerance=" << nonlinear_tolerance << '\n';
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
				+std::to_string(nonlinear)+"; residual "+std::to_string(result.residual));
		total_iterations += result.iterations;
		peak_gpu_used = std::max(peak_gpu_used, result.device_used_gib);
		double update_norm = 0.0;
		Check(cublasDnrm2(blas, static_cast<int>(host.nodes*4), update.data(), 1, &update_norm),
			"Navier-Stokes update norm");
		const double one = 1.0;
		Check(cublasDaxpy(blas, static_cast<int>(host.nodes*4), &one, update.data(), 1,
			state.data(), 1), "Navier-Stokes state update");
		std::cout << "newton=" << nonlinear << " residual_l2=" << residual_norm
			<< " update_l2=" << update_norm << " linear_iterations=" << result.iterations
			<< " linear_residual=" << result.residual
			<< " singular_diagonal_blocks=" << result.singular_diagonal_blocks
			<< " assembly_s=" << std::chrono::duration<double>(assembly_end-assembly_start).count()
			<< " linear_s=" << std::chrono::duration<double>(linear_end-linear_start).count() << '\n';
		if (update_norm < 1e-10) {
			converged = true;
			break;
		}
	}
	if (!converged)
		throw std::runtime_error("Navier-Stokes nonlinear solve reached MAX_NEWTON without convergence");
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
	if (argc >= 6) {
		WriteNavierStokes(argv[5], output);
		WriteNavierStokesVtk(case_dir/"controlmesh.vtk", std::string(argv[5])+".vtk", output);
	}
	std::cout << "navier_stokes_cuda nodes=" << host.nodes << " elements=" << host.elements()
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
			throw std::runtime_error("usage: iga_cuda device-info|mesh-check|transport|navier-stokes ...");
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
		if (command == "navier-stokes") return iga::cuda::NavierStokes(argc, argv);
		throw std::runtime_error("unknown command: "+command);
	} catch (const std::exception& error) {
		std::cerr << "iga_cuda: " << error.what() << '\n';
		return 1;
	}
}
