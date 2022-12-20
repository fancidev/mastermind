/// @defgroup Choose Combinatorics Routines
/// @ingroup util

#ifndef UTILITIES_CHOOSE_HPP
#define UTILITIES_CHOOSE_HPP

#include <cassert>

namespace util {

/// Returns the number of ways to choose @c r elements out of @c n 
/// elements, optionally with order and/or repetition.
/// @ingroup Choose
template <typename TVal, typename TArg>
TVal choice(TArg n, TArg r, bool order = false, bool rep = false)
{
	assert(n >= 0 && r >= 0);

	if (order)
	{
		if (rep) // c = n^r
		{
			TVal c = 1;
			for (TArg i = 0; i < r; ++i)
				c *= TVal(n);
			return c;
		}
		else // c = P(n,r) = n*(n-1)*...*(n-r+1)
		{
			assert(n >= r);
			TVal c = 1;
			for (TArg i = n; i > n - r; --i)
				c *= TVal(i);
			return c;
		}
	}
	else
	{
		assert(0);
		return 0;
	}
}

} // namespace util

#endif // UTILITIES_CHOOSE_HPP

#if 0
template <size_t MaxSize, class Iter1, class Iter2, class Func>
void generate_permutations_recursion(
	Iter1 from, Iter2 to, 
	size_t n, size_t r, size_t level,
	std::bitset<MaxSize> &mask, Func f)
{
	for (size_t i = 0; i < n; ++i)
	{
		if (!mask[i])
		{
			if (level == r - 1)
			{
				to[level] = from[i];
				f();
			}
			else
			{
				mask[i] = true;
				to[level] = from[i];
				generate_permutations_recursion(from, to, n, r, level+1, mask, f);
				mask[i] = false;
			}
		}
	}
}

/// Generates all permutatations of choosing @c r objects from @c n 
/// objects. The objects to choose from is supplied in 
/// <code>[from_first,from_last)</code>. The objects chosen are
/// stored in <code>[to_first,to_last)</code>. 
///
/// For each permutation generated, the callback function @c f
/// is invoked with no parameters. If <code>r = 0</code>, @c f
/// is invoked once.
///
/// @c MaxSize specifies an upper bound of @c n.
template <size_t MaxSize, class Iter1, class Iter2, class Func>
void generate_permutations(
	Iter1 from_first, Iter1 from_last, 
	Iter2 to_first, Iter2 to_last,
	Func f)
{
	size_t n = from_last - from_first;
	size_t r = to_last - to_first;
	assert(n >= r && r >= 0 && MaxSize >= n);

	if (r == 0)
	{
		f();
	}
	else
	{
		std::bitset<MaxSize> mask;
		generate_permutations_recursion(from_first, to_first, n, r, 0, mask, f);
	}
}
#endif

#if 0
/// Utility class that computes the lexicographical index of a codeword.
/// For performance reason, we ALWAYS computes the index as if repetition
/// was allowed in the colors. This improves speed at the cost of a little
/// waste of memory.
class CodewordIndexer
{
	Rules _rules;
	int _weights[MM_MAX_COLORS];

public:

	CodewordIndexer(const Rules &rules) : _rules(rules)
	{
#if 0
		if (_rules.repeatable())
		{
#endif
			int w = 1;
			for (int i = _rules.pegs() - 1; i >= 0; --i)
			{
				_weights[i] = w;
				w *= _rules.colors();
			}
#if 0
		}
		else
		{
			int w = 1;
			for (int i = _rules.pegs() - 1; i >= 0; --i)
			{
				_weights[i] = w;
				w *= (_rules.colors() - i);
			}
		}
#endif
	}

	int size() const { return _weights[0] * _rules.colors(); }

	int operator()(const Codeword &c) const 
	{
#if 0
		if (_rules.repeatable())
		{
#endif
			int index = 0;
			for (int i = 0; i < _rules.pegs(); ++i)
				index += c[i] * _weights[i];
			return index;
#if 0
		}
		else
		{
			unsigned short bitmask = 0;
			int index = 0;
			for (int i = 0; i < _rules.pegs(); ++i)
			{
				int t = c[i] - util::intrinsic::pop_count(bitmask & ((1 << c[i]) - 1));
				index += t * _weights[i];
				bitmask |= (1 << c[i]);
			}
			return index;
		}
#endif
	}
};
#endif
