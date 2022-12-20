/// @defgroup Range Iterator Range
/// @ingroup util

#ifndef UTILITIES_RANGE_HPP
#define UTILITIES_RANGE_HPP

#include <utility>
#include <iterator>

namespace util {

/// Represents a range specified by a pair of <code>[begin,end)</code> 
/// iterators.
/// @ingroup Range
template <class Iter>
struct range : private std::pair<Iter,Iter>
{
	/// Type of the value stored in this range.
	typedef typename std::iterator_traits<Iter>::value_type value_type;

	/// Type of reference to the value stored in this range.
	typedef typename std::iterator_traits<Iter>::reference reference;

	/// Constructs a range from a pair of iterators.
	range(Iter first, Iter last)
		: std::pair<Iter,Iter>(first, last) { }

	/// Constructs a range from another range.
	template <class Iter2>
	range(const range<Iter2> &r)
		: std::pair<Iter,Iter>(r.begin(), r.end()) { }

	/// Constructs the entire range of a container.
	template <class Container>
	range(Container &c)
		: std::pair<Iter,Iter>(c.begin(), c.end()) { }

	/// Returns the begin iterator of the range.
	Iter begin() const { return std::pair<Iter,Iter>::first; }

	/// Returns the end iterator of the range.
	Iter end() const { return std::pair<Iter,Iter>::second; }

	/// Tests whether the range is empty.
	bool empty() const { return begin() == end(); }

	/// Returns the number of elements in the range.
	/// The underlying iterator type must support random access.
	size_t size() const { return end() - begin(); }

	/// Returns the element at a given index. 
	/// The underlying iterator type must support random access.
	reference operator[](size_t i) { return std::pair<Iter,Iter>::first[i]; }
};

} // namespace util


#endif // UTILITIES_RANGE_HPP
