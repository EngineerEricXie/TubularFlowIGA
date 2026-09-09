#ifndef IGA_CHECKED_TEXT_HPP
#define IGA_CHECKED_TEXT_HPP

#include <array>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace iga {

// Standard output streams do not throw by default. Flush within the caller's
// error boundary so a failed buffered write cannot escape as apparent success.
// This does not change the caller's exception mask or stream ownership.
inline void FlushCheckedText(std::ostream& output)
{
	output.flush();
	if (!output) throw std::runtime_error("cannot write complete text stream");
}

// Read through the istream so a filebuf error updates the state we check.
// Direct streambuf insertion can leave input.good() true after a failed copy.
// Empty input and the final short read are normal; no partial result escapes
// a read or allocation failure. Callers retain their existing open/path policy.
inline std::string ReadCheckedText(std::istream& input)
{
	std::string result;
	std::array<char, 4096> buffer{};
	for (;;) {
		input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
		if (input.bad() || (input.fail() && !input.eof()))
			throw std::runtime_error("cannot read complete text stream");
		result.append(buffer.data(), static_cast<std::size_t>(input.gcount()));
		if (input.eof()) return result;
	}
}

} // namespace iga

#endif
