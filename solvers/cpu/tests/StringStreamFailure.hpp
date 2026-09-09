#ifndef IGA_TEST_STRING_STREAM_FAILURE_HPP
#define IGA_TEST_STRING_STREAM_FAILURE_HPP

#include "CollectiveFailure.hpp"
#include <iostream>
#include <locale>
#include <new>
#include <sstream>

namespace iga_test {

// A test-only facet throws while an actual ostringstream formats a number.
// Existing file streams and MPI/PETSc allocations are unaffected. Installing
// this locale is restricted to the test's main thread with no active workers.
class StringStreamFailure : public std::num_put<char> {
public:
	StringStreamFailure(bool enabled, int mode, std::string prefix = {}, long number = -1)
		: enabled_(enabled), mode_(mode), prefix_(std::move(prefix)), number_(number) {}
	bool Injected() const { return injected_; }
protected:
	void MaybeFail(std::ios_base& output, long double number) const
	{
		if (!enabled_ || injected_) return;
		if (number_ >= 0 && number != number_) return;
		const auto* stream = dynamic_cast<std::ostringstream*>(&output);
		if (!stream || stream->str().compare(0, prefix_.size(), prefix_) != 0) return;
		injected_ = true;
		if (mode_ == 0) throw std::runtime_error("injected stringstream formatting failure");
		if (mode_ == 1) throw std::bad_alloc();
		throw 7;
	}

#define IGA_TEST_NUM_PUT(Type) \
	iter_type do_put(iter_type iterator, std::ios_base& output, char fill, Type value) const override \
	{ MaybeFail(output, value); return std::num_put<char>::do_put(iterator, output, fill, value); }
	IGA_TEST_NUM_PUT(bool)
	IGA_TEST_NUM_PUT(long)
	IGA_TEST_NUM_PUT(unsigned long)
	IGA_TEST_NUM_PUT(long long)
	IGA_TEST_NUM_PUT(unsigned long long)
	IGA_TEST_NUM_PUT(double)
	IGA_TEST_NUM_PUT(long double)
#undef IGA_TEST_NUM_PUT
private:
	bool enabled_;
	int mode_;
	std::string prefix_;
	long number_;
	mutable bool injected_ = false;
};

class ScopedStringStreamFailure {
public:
	ScopedStringStreamFailure(bool enabled, int mode, const std::string& prefix = {}, long number = -1)
		: facet_(new StringStreamFailure(enabled, mode, prefix, number)),
		  installed_(std::locale(), facet_), previous_(std::locale::global(installed_)) {}
	~ScopedStringStreamFailure() { std::locale::global(previous_); }
	bool Injected() const { return facet_->Injected(); }
	ScopedStringStreamFailure(const ScopedStringStreamFailure&) = delete;
	ScopedStringStreamFailure& operator=(const ScopedStringStreamFailure&) = delete;
private:
	StringStreamFailure* facet_;
	std::locale installed_, previous_;
};

template<class Work>
void ExpectStringStreamFailure(MPI_Comm comm, const char* stage, Work&& work)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	for (int mode = 0; mode < 3; ++mode) {
		std::string diagnostic;
		bool injected = false;
		{
			ScopedStringStreamFailure fault(rank == ranks-1, mode);
			try { work(); }
			catch (const std::exception& error) { diagnostic = error.what(); }
			injected = fault.Injected();
		}
		iga::CollectiveLocalStage(comm, "stringstream failure checks", [&] {
			if (injected != (rank == ranks-1)) throw std::runtime_error("stringstream injection missed target");
			if (diagnostic.find(std::string(stage)+": rank ") == std::string::npos)
				throw std::runtime_error("wrong stringstream failure boundary: "+diagnostic);
		});
		iga::RequireCollectiveSameText(comm, "stringstream common diagnostic", diagnostic);
	}
	if (rank == 0) std::cout << "stringstream stage=" << stage << " ranks=" << ranks << " faults=3 passed\n";
}

} // namespace iga_test
#endif
