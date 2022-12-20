#include <iostream>
#include <algorithm>

#include "Codeword.hpp"
#include "util/simd.hpp"

namespace Mastermind {

/// Tests whether the codeword contains any color more than once.
bool Codeword::has_repetition() const
{
#if 1
	int mask = util::simd::byte_mask(
		*(util::simd::simd_t<int8_t,16>*)this > (int8_t)1);
	return (mask & ((1 << MM_MAX_COLORS) - 1)) != 0;
#else	
	for (int i = 0; i < MM_MAX_COLORS; ++i)
	{
		if (_counter[i] > 1)
			return true;
	}
	return false;
#endif
}

bool Codeword::conforming(const Rules &rules) const
{
	if (!rules)
		return false;

	// Check pegs.
	for (int p = 0; p < rules.pegs(); ++p)
	{
		if (_digit[p] < 0)
			return false;
	}
	for (int p = rules.pegs(); p < MM_MAX_PEGS; ++p)
	{
		if (_digit[p] >= 0)
			return false;
	}

	// Check colors.
	if (!rules.repeatable())
	{
		for (int c = 0; c < rules.colors(); ++c)
		{
			if (_counter[c] > 1)
				return false;
		}
	}
	for (int c = rules.colors(); c < MM_MAX_COLORS; ++c)
	{
		if (_counter[c] > 0)
			return false;
	}
	return true;
}

static const char *codeword_alphabet = "1234567890abcdef";

std::ostream& operator << (std::ostream &os, const Codeword &c)
{
	// Build a string.
	char s[MM_MAX_PEGS + 1] = {0};
	for (int k = 0; k < MM_MAX_PEGS; k++)
	{
		int d = c[k];
		if (d == Codeword::EmptyColor)
			break;
		assert(d < 16);
		s[k] = codeword_alphabet[(int)d];
	}
	return os << s;
}

std::istream& operator >> (std::istream &is, Codeword &codeword)
{
	Codeword ret;
	const char *first = codeword_alphabet;
	const char *last = first + 16;

	// Read the first character, skipping leading whitespaces.
	char c = 0;
	if (!(is >> c))
		return is;

	// Process characters.
	int i = 0;
	bool ok = true;
	do
	{
		const char *p = std::find(first, last, c);
		if (p == last) // invalid input
		{
			// Unget character.
			is.unget();
			break;
		}
		if (i < MM_MAX_PEGS && (p - first) < MM_MAX_COLORS)
		{
			ret.set(i, (int)(p - first));
		}
		else
		{
			ok = false;
		}
		++i;
	} while ((c = (char)is.get()) != EOF);

	// Try get the rules associated with this stream.
	Rules rules = getrules(is);

	// If the codeword is successfully read, return it.
	// Otherwise, sets the error flag.
	if (ok && !ret.IsEmpty() && (!rules || ret.conforming(rules)))
	{
		if (c == EOF)
			is.clear(std::ios_base::eofbit);
		codeword = ret;
	}
	else
	{
		is.setstate(std::ios_base::failbit);
	}
	return is;
}

} // namespace Mastermind
