#ifndef MASTERMIND_PERMUTATION_HPP
#define MASTERMIND_PERMUTATION_HPP

#include <iostream>
#include <cassert>
#include <cstdint>
#include <numeric>
#include "Codeword.hpp"

namespace Mastermind {

/**
 * Permutation that permutes the pegs and colors in a codeword.
 * It is the composition of a peg permutation and a color permutation.
 * @ingroup type
 */
struct CodewordPermutation
{
	/// Permuted value of a color. For example, <code>color[1] = 2</code>
	/// means color @c 1 is changed to color @c 2 after the permutation.
	int8_t color[MM_MAX_COLORS];

	/// Permuted value of a peg. For example, <code>peg[0] = 3</code>
	/// means the peg in position @c 0 is moved to position @c 3 after
	/// the permutation.
	int8_t peg[MM_MAX_PEGS];

public:

	/// Creates an identity permutation.
	CodewordPermutation()
	{
		std::iota(color + 0, color + MM_MAX_COLORS, (int8_t)0);
		std::iota(peg + 0, peg + MM_MAX_PEGS, (int8_t)0);
	}

#if 0
	/// Returns the inverse of the permutation.
	CodewordPermutation inverse() const
	{
		CodewordPermutation ret(_rules);
		for (int i = 0; i < _rules.pegs(); ++i)
		{
			if (peg[i] >= 0)
				ret.peg[peg[i]] = i;
		}
		for (int i = 0; i < _rules.colors(); ++i)
		{
			if (color[i] >= 0)
				ret.color[color[i]] = i;
		}
		return ret;
	}
#endif

	/// Permutes the pegs and colors in a codeword.
	Codeword permute(const Codeword &w) const
	{
		Codeword ret;
		for (int i = 0; i < MM_MAX_PEGS && w[i] != Codeword::EmptyColor; ++i)
		{
			ret.set(peg[i], color[w[i]]);
		}
		return ret;
	}

	/// Permutes the pegs in a codeword.
	Codeword permute_pegs(const Codeword &w) const
	{
		Codeword ret;
		for (int i = 0; i < MM_MAX_PEGS /* && w[i] != 0xFF */; ++i)
			ret.set(peg[i], w[i]);
		return ret;
	}

#if 0
	Codeword permute_colors(const Codeword &w) const
	{
		Codeword ret;
		for (int i = 0; i < _rules.pegs(); ++i)
			ret.set(i, color[w[i]]);
		return ret;
	}
#endif
};

/// Outputs a codeword permutation to a stream.
inline std::ostream& operator << (std::ostream& os, const CodewordPermutation &p)
{
	// Output peg permutation.
	os << "(";
	for (int i = 0; i < MM_MAX_PEGS; ++i)
	{
		if (i > 0)
			os << ' ';
		os << (size_t)p.peg[i];
	}
	os << ") o (";

	// Output color permutation.
	for (int i = 0; i < MM_MAX_COLORS; ++i)
	{
		if (i > 0)
			os << ' ';
		os << (size_t)p.color[i];
	}
	os << ")";
	return os;
}

} // namespace Mastermind

#endif // MASTERMIND_PERMUTATION_HPP
