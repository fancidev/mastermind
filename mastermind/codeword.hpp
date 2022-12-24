#ifndef MASTERMIND_CODEWORD_HPP
#define MASTERMIND_CODEWORD_HPP

#include <iostream>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <cassert>

#include "Rules.hpp"
#include "util/intrinsic.hpp"

namespace Mastermind {

/// <summary>Represents a codeword (such as 2587).</summary>
/// <remarks>
/// For performance reasons, a codeword must be aligned on a 16-byte boundary.
/// </remarks>
class 
#ifdef _MSC_VER
	__declspec(align(16))
#else
	__attribute__ ((aligned (16)))
#endif
	Codeword
{
	/// The number of occurrences of each color. If a color is not included,
	/// the corresponding value is zero.
	int8_t _counter[MM_MAX_COLORS];

	/// The (zero-based) color on each peg. If a peg is empty, the value is
	/// <code>(int8_t)(-1)</code>.
	int8_t _digit[MM_MAX_PEGS];

public:

	/// Constant representing an 'empty' color.
	static const int EmptyColor = -1;
		
	/// <summary>Creates an empty codeword.</summary>
	Codeword()
	{
		std::memset(_counter, 0, sizeof(_counter));
		std::memset(_digit, -1, sizeof(_digit));
	}

	/// <summary>Checks whether the codeword is empty.</summary>
    /// <return><code>true</code> if the codeword is empty, otherwise
    /// <code>false</code>.</return>
	bool IsEmpty() const { return _digit[0] < 0; }

#if 0
	/// Tests whether the codeword is non-empty.
	operator void* () const { return empty()? 0 : (void*)this; }

	/// Tests whether the codeword is empty.
	bool operator ! () const { return empty(); }
#endif

	/// <summary>Gets the color in the given peg.</summary>
    /// <param name="peg">Zero-based index of the peg.</param>
    /// <return>Color in the given peg.</return>
	int operator [] (int peg) const
	{
		assert(peg >= 0 && peg < MM_MAX_PEGS);
		return _digit[peg];
	}

	/// Sets the color on a given peg.
	void set(int peg, int color)
	{
		assert(peg >= 0 && peg < MM_MAX_PEGS);
		assert(color == -1 || (color >= 0 && color < MM_MAX_COLORS));
		if (_digit[peg] >= 0)
			--_counter[(int)_digit[peg]];
		if ((_digit[peg] = (char)color) >= 0)
			++_counter[color];
	}

	/// Returns the number of occurrences of a given color.
	int count(int color) const
	{
		assert(color >= 0 && color < MM_MAX_COLORS);
		return _counter[color];
	}

	/// Tests whether the codeword contains any color more than once.
	bool has_repetition() const;

#if 0
	/// Returns the number of pegs in the codeword.
	int pegs() const
	{
		int n = 0;
		for (int i = 0; i < MM_MAX_PEGS; ++i)
		{
			if (_digit[i] >= 0)
				++n;
		}
		return n;
	}
#endif

	/// Checks whether this codeword conforms to the supplied rules.
	bool conforming(const Rules &rules) const;

	/// Type of the packed value.
	typedef uint32_t compact_type;

	/// Packs a codeword into a 4-byte representation.
	compact_type pack() const
	{
		uint32_t w = 0xffffffff;
		for (int i = 0; i < MM_MAX_PEGS; i++)
		{
			char d = _digit[i];
			if (d < 0)
				break;
			w <<= 4;
			w |= d;
		}
		return w;
	}

	/// Unpacks a codeword from a 4-byte representation.
	static Codeword unpack(compact_type w)
	{
		Codeword c;
		int i = 0;
		for (int nibble = 7; nibble >= 0; --nibble)
		{
			char d = (w >> (nibble*4)) & 0xF;
			if (d != 0xF)
				c.set(i++, d);
		}
		return c;
	}
};

/// Tests whether two codewords are equal.
/// @ingroup Codeword
inline bool operator == (const Codeword &a, const Codeword &b)
{
	return memcmp(&a, &b, sizeof(Codeword)) == 0;
}

/// Tests whether two codewords are not equal.
/// @ingroup Codeword
inline bool operator != (const Codeword &a, const Codeword &b)
{
	return ! operator == (a, b);
}

/// Outputs a codeword to a stream.
/// @ingroup Codeword
std::ostream& operator << (std::ostream &os, const Codeword &c);

/// Inputs a codeword from a stream.
/// @ingroup Codeword
std::istream& operator >> (std::istream &is, Codeword &c);

} // namespace Mastermind

#endif // MASTERMIND_CODEWORD_HPP
