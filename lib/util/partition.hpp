#ifndef UTILITIES_PARTITION_HPP
#define UTILITIES_PARTITION_HPP

#include <cassert>
#include "range.hpp"
#include "frequency_table.hpp"

namespace util {

/// Fixed-size array of consecutive ranges that represents the cells in the
/// partition of a larger range.
/// @tparam Iter Iterator type for the underlying range.
/// @tparam Capacity The maximum number of cells supported.
template <class Iter, int Capacity>
class partition_cells
{
	size_t _size;
	Iter _begin[Capacity+1];

public:

	/// Creates an empty partition, which contains no cells.
	partition_cells() : _size(0) { }

	/// Creates a unit partition, which contains only one cell.
	partition_cells(range<Iter> cell) : _size(1) 
	{
		_begin[0] = cell.begin();
		_begin[1] = cell.end();
	}

	/// Creates a partition of a range from an associated frequency table.
	/// @tparam TKey Key type of the frequency table.
	/// @tparam TVal Value type of the frequency table.
	/// @param all Range to partition.
	/// @param freq Frequency table which contains the sizes of the cells
	///      in order.
	template <class TKey, class TVal>
	partition_cells(range<Iter> all, const frequency_table<TKey,TVal,Capacity> &freq)
		: _size(freq.size())
	{
		_begin[0] = all.begin();
		for (size_t j = 0; j < freq.size(); ++j)
		{
			_begin[j+1] = _begin[j] + (size_t)freq[j];
		}
		assert(_begin[_size] == all.end());
	}

	/// Returns the number of cells in this partition.
	size_t size() const { return _size; }

	/// Returns the given cell.
	range<Iter> operator [] (size_t i) const
	{
		assert(i < _size);
		return range<Iter>(_begin[i], _begin[i+1]);
	}
};

} // namespace util

#endif // UTILITIES_PARTITION_HPP
