#include "Sha256.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

void Check(const char* input, const char* expected)
{
	iga::Sha256 hash;
	hash.Append(input, std::strlen(input));
	assert(hash.Hex() == expected);
}

template <class Function>
void RequireRejected(Function&& function)
{
	bool rejected = false;
	try {
		function();
	} catch (const std::exception&) {
		rejected = true;
	}
	assert(rejected);
}

} // namespace

int main()
{
	Check("", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
	Check("abc", "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
	Check("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
		"248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
	iga::Sha256 million;
	for (int i = 0; i < 1000000; ++i) million.Append("a", 1);
	assert(million.Hex() == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");

	const char payload[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
	iga::Sha256 whole;
	whole.Append(payload, sizeof(payload)-1);
	const auto expected = whole.Hex();
	for (const auto split : {1U, 55U, 56U, 63U, 64U, 65U}) {
		iga::Sha256 partial;
		for (std::size_t offset = 0; offset < sizeof(payload)-1; offset += split)
			partial.Append(payload+offset,
				std::min<std::size_t>(split, sizeof(payload)-1-offset));
		assert(partial.Hex() == expected);
	}
	const unsigned char embedded_input[] = {0, 1, 2, 0, 255};
	iga::Sha256 embedded;
	embedded.Append(embedded_input, sizeof(embedded_input));
	assert(embedded.Hex() ==
		"ef7e301027f931dfba06c7ded4ef305797f43cc115a664f9af9b57d08c3172c2");

	iga::Sha256 canonical;
	canonical.AppendLittleEndian32(0x01020304);
	canonical.AppendLittleEndian64(0x0102030405060708ULL);
	canonical.AppendNormalizedDouble(0.0);
	const std::array<unsigned char, 20> canonical_bytes{{4, 3, 2, 1, 8, 7, 6, 5,
		4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0}};
	iga::Sha256 explicit_bytes;
	explicit_bytes.Append(canonical_bytes.data(), canonical_bytes.size());
	assert(canonical.Hex() == explicit_bytes.Hex());
	iga::Sha256 negative_zero;
	negative_zero.AppendNormalizedDouble(-0.0);
	iga::Sha256 positive_zero;
	positive_zero.AppendNormalizedDouble(0.0);
	assert(negative_zero.Hex() == positive_zero.Hex());
	RequireRejected([] { iga::Sha256 hash; hash.Append(nullptr, 1); });
	RequireRejected([] {
		iga::Sha256 hash;
		hash.AppendNormalizedDouble(std::numeric_limits<double>::infinity());
	});
	iga::Sha256 empty_null;
	empty_null.Append(nullptr, 0);
	assert(empty_null.Hex() ==
		"e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
	std::cout << "SHA-256 tests passed\n";
}
