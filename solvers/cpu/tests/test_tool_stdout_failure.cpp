#ifdef IGA_TEST_CUDA_STDOUT
#define main NativeToolMain
#include "../../cuda/src/iga_cuda.cu"
#undef main
#elif defined(IGA_TEST_MESH_STDOUT)
#define main NativeToolMain
#include "../src/iga_mesh_check.cpp"
#undef main
#elif defined(IGA_TEST_FLOW_STDOUT)
#define main NativeToolMain
#include "../src/iga_navier_stokes.cpp"
#undef main
#elif defined(IGA_TEST_TRANSPORT_STDOUT)
#define main NativeToolMain
#include "../src/iga_solve.cpp"
#undef main
#elif defined(IGA_TEST_LEGACY_STDOUT)
#define main NativeToolMain
#include "../src/iga_transport.cpp"
#undef main
#elif defined(IGA_TEST_ONE_D_STDOUT)
#define main NativeToolMain
#include "../../one_d/src/iga_1d.cpp"
#undef main
#elif defined(IGA_TEST_GRAPH_STDOUT)
#define main NativeToolMain
#include "../../coupling/src/iga_1d_3d_bifurcation.cpp"
#undef main
#elif defined(IGA_TEST_SEQUENTIAL_STDOUT)
#define main NativeToolMain
#include "../../coupling/src/iga_1d_3d_explicit.cpp"
#undef main
#else
#define main NativeToolMain
#include "../src/iga_assembly_smoke.cpp"
#undef main
#endif

#include <cstdlib>
#include <cstring>
#include <new>
#include <streambuf>

namespace {
class FailedOutput : public std::streambuf {
public:
	FailedOutput(std::streambuf* original, std::string target, std::string mode)
		: original_(original), target_(std::move(target)), mode_(std::move(mode)) {}
	bool Injected() const { return injected_; }
protected:
	std::streamsize xsputn(const char* data, std::streamsize size) override
	{
		if (std::string(data, static_cast<std::size_t>(size)).find(target_) != std::string::npos)
			armed_ = true;
		if (armed_ && mode_ != "flush") {
			injected_ = true;
			if (mode_ == "exception") throw std::runtime_error("injected stdout write error");
			if (mode_ == "allocation") throw std::bad_alloc();
			if (mode_ == "nonstandard") throw 7;
			return 0;
		}
		return original_->sputn(data, size);
	}
	int_type overflow(int_type value) override
	{
		if (traits_type::eq_int_type(value, traits_type::eof())) return traits_type::not_eof(value);
		const auto c = traits_type::to_char_type(value);
		return xsputn(&c, 1) == 1 ? value : traits_type::eof();
	}
	int sync() override
	{
		if (armed_ && mode_ == "flush") { injected_ = true; return -1; }
		return original_->pubsync();
	}
private:
	std::streambuf* original_;
	std::string target_, mode_;
	bool armed_ = false, injected_ = false;
};
}

int main(int argc, char** argv)
{
	const char* rank = std::getenv("OMPI_COMM_WORLD_RANK");
	const char* target = std::getenv("IGA_TEST_STDOUT_TARGET");
	const char* mode = std::getenv("IGA_TEST_STDOUT_MODE");
	if (!rank || !target || !mode) return 3;
	auto* original = std::cout.rdbuf();
	FailedOutput failed(original, target, mode);
	if (std::strcmp(rank, "0") == 0) std::cout.rdbuf(&failed);
	const int status = NativeToolMain(argc, argv);
	std::cout.rdbuf(original);
	std::cout.clear();
	std::cerr << "stdout_test rank=" << rank << " injected=" << failed.Injected()
		<< " native_status=" << status << '\n';
	return status;
}
