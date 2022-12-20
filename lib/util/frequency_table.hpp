/// @defgroup FreqTable Frequency Table
/// @ingroup util

#ifndef UTILITIES_FREQUENCY_TABLE_HPP
#define UTILITIES_FREQUENCY_TABLE_HPP

#include <cassert>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cstring>

namespace util {

/// Represents a frequency table of (value, frequency) pairs where
/// the values are mapped to and from zero-based indices.
/// @ingroup FreqTable
template <class TKey, class TVal, size_t Capacity>
class frequency_table
{
	TVal _freq[Capacity];
	size_t _count;

public:

	/// Constructs an empty frequency table.
	frequency_table() : _count(0) { }

	/// Constructs a frequency table of the given size.
	frequency_table(size_t n) : _count(n) 
	{
		assert(n <= Capacity);
		// Note: gcc does not optimize std::fill() with memset.
		// So we need to use memset explicitly. Note that this
		// requires TVal to be suitably typed.
#if 0
		std::fill(_freq + 0, _freq + n, TVal(0));
#else
		memset(_freq, 0, sizeof(TVal)*n);
#endif
	}

	/// Returns the size of the frequency table.
	size_t size() const { return _count; }

	/// Sets the size of the frequency table and optionally resets all 
	/// frequency counts to zero.
	void resize(size_t n, bool zero_memory = true)
	{
		assert(n <= Capacity);
		_count = n;
		if (zero_memory)
		{
#if 0
			std::fill(_freq + 0, _freq + n, TVal(0));
#else
			memset(_freq, 0, sizeof(TVal)*n);
#endif
		}
	}

	/// Returns a pointer to the first frequency value.
	TVal* data() { return _freq; }

	/// Returns a const pointer to the first frequency value.
	const TVal* data() const { return _freq; }

	/// Returns an iterator that points to the beginning of the table.
	TVal* begin() { return _freq + 0; }

	/// Returns a const iterator that points to the beginning of the table.
	const TVal* begin() const { return _freq + 0; }

	/// Returns an iterator that points past the end of the table.
	TVal* end() { return _freq + _count; }

	/// Returns a const iterator that points past the end of the table.
	const TVal* end() const { return _freq + _count; }

	/// Returns a reference to the frequency of the <code>k</code>-th element.
	TVal& operator [] (size_t k)
	{
		assert(k < _count);
		return _freq[k];
	}

	/// Returns the frequency of the <code>k</code>-th element.
	TVal operator [] (size_t k) const
	{
		assert(k < _count);
		return _freq[k];
	}

	/// Returns the maximum frequency value.
	TVal max() const
	{
		return _count == 0? 0 : *std::max_element(begin(), end());
	}

	/// Returns the number of non-zero categories.
	size_t nonzero_count() const
	{
		return std::count_if(begin(), end(), [](TVal f) -> bool {
			return f > 0;
		});
	}

};

/// Outputs a frequency table to a stream.
/// @ingroup FreqTable
template <class TKey, class TVal, size_t Capacity>
std::ostream&
operator << (std::ostream& os, const frequency_table<TKey,TVal,Capacity> &f)
{
	size_t total = 0;
	for (size_t i = 0; i < f.size(); i++)
	{
		if (f[i] != 0)
		{
			os << TKey(i) << " => " << f[i] << std::endl;
			total += f[i];
		}
	}
	os << "Total: " << total << std::endl;
	return os;
}

/// Compares two frequency tables lexicographcially.
/// @ingroup FreqTable
template <class TKey, class TVal, size_t Capacity>
inline bool operator < (
	const frequency_table<TKey,TVal,Capacity> &t1,
	const frequency_table<TKey,TVal,Capacity> &t2)
{
	return std::lexicographical_compare(t1.begin(), t1.end(), t2.begin(), t2.end());
}

} // namespace util

#endif // UTILITIES_FREQUENCY_TABLE_HPP
