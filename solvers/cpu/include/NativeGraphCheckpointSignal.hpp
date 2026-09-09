#ifndef IGA_NATIVE_GRAPH_CHECKPOINT_SIGNAL_HPP
#define IGA_NATIVE_GRAPH_CHECKPOINT_SIGNAL_HPP
#include <csignal>
#include <stdexcept>

namespace iga {
namespace native_checkpoint_signal_detail {
inline volatile std::sig_atomic_t requested = 0;
inline void Request(int) noexcept
{
	requested = 1;
}
}
// Install only when checkpoint output is enabled. The handler does no MPI, I/O
// or allocation; the runner agrees on the request after an accepted macro-step.
class NativeGraphCheckpointSignal {
public:
	NativeGraphCheckpointSignal()
	{
		struct sigaction action{}; action.sa_handler = native_checkpoint_signal_detail::Request;
		sigemptyset(&action.sa_mask); action.sa_flags = SA_RESTART;
		native_checkpoint_signal_detail::requested = 0;
		if (::sigaction(SIGUSR1, &action, &previous_) != 0) throw std::runtime_error("cannot install checkpoint SIGUSR1 handler");
	}
	NativeGraphCheckpointSignal(const NativeGraphCheckpointSignal&) = delete;
	NativeGraphCheckpointSignal& operator=(const NativeGraphCheckpointSignal&) = delete;
	~NativeGraphCheckpointSignal()
	{
		::sigaction(SIGUSR1, &previous_, nullptr);
	}
	bool Requested() const noexcept
	{
		return native_checkpoint_signal_detail::requested != 0;
	}
private:
	struct sigaction previous_{};
};
} // namespace iga
#endif
