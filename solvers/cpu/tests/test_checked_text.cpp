#include "CheckedText.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <new>
#include <sstream>
#include <streambuf>

namespace {

void Check(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}

class FailingBuffer : public std::streambuf {
public:
	FailingBuffer(std::string prefix, int mode) : prefix_(std::move(prefix)), mode_(mode)
	{
		setg(prefix_.data(), prefix_.data(), prefix_.data()+prefix_.size());
	}
protected:
	int_type underflow() override
	{
		if (mode_ == 0) throw std::runtime_error("injected read error");
		if (mode_ == 1) throw std::bad_alloc();
		throw 7;
	}
private:
	std::string prefix_;
	int mode_;
};

} // namespace

int main()
{
	try {
		for (const std::size_t size : {0, 1, 4095, 4096, 4097, 8192, 12000}) {
			std::string expected(size, 'x');
			if (size > 5) expected[5] = '\0';
			std::istringstream input(expected);
			Check(iga::ReadCheckedText(input) == expected, "text read changed bytes or EOF behavior");
		}
		const std::string prefix = "{\"valid\":true}"+std::string(9000, ' ');
		for (int mode = 0; mode < 3; ++mode) {
			// Demonstrate that testing only input.good() after streambuf
			// insertion can accept a valid JSON prefix despite a read failure.
			FailingBuffer old_buffer(prefix, mode);
			std::istream old_input(&old_buffer);
			std::ostringstream old_copy;
			old_copy << old_input.rdbuf();
			Check(old_input.good() && old_copy.fail() && old_copy.str() == prefix,
				"unchecked-copy failure reproduction changed");
			for (const auto& initial : {std::string{}, prefix}) {
				FailingBuffer buffer(initial, mode);
				std::istream input(&buffer);
				bool failed = false;
				try { (void)iga::ReadCheckedText(input); }
				catch (...) { failed = true; }
				Check(failed, "checked reader returned partial text after a read failure");
				std::istringstream retry(prefix);
				Check(iga::ReadCheckedText(retry) == prefix, "healthy retry changed text");
			}
		}
		std::istringstream failed_input("unused");
		failed_input.setstate(std::ios::failbit);
		bool failed = false;
		try { (void)iga::ReadCheckedText(failed_input); }
		catch (const std::exception&) { failed = true; }
		Check(failed, "reader accepted a stream already in failure state");
		std::cout << "checked text: 7 byte/EOF cases, 3 old-copy reproductions, 6 read faults/retries, prefailed stream passed\n";
		return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
