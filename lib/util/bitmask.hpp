/// @defgroup BitMask Bit-Mask of fixed size
/// @ingroup util

#ifndef UTILITIES_BITMASK_HPP
#define UTILITIES_BITMASK_HPP

#include <cassert>
#include "intrinsic.hpp"

namespace util {

/**
 * Represents a bitmask of fixed size.
 * This class serves as a simpler and faster alternative to
 * <code>std::bitset<N></code>.
 * @ingroup BitMask
 */
template <class T, size_t Bits>
class bitmask
{
public:

	static_assert(Bits <= sizeof(T)*8, "The storage type does not have enough bits.");
	
	typedef T value_type;

private:

	value_type _value;

public:

	/// Creates an empty bitmask.
	bitmask() : _value(0) { }

	/// Creates a bitmask using the supplied mask.
	explicit bitmask(value_type value) : _value(value) { }

	/// Gets the internal value of the mask.
	value_type value() const { return _value; }

	/// Tests a given bit.
	bool operator [] (size_t bit) const
	{
		assert(bit <= Bits);
		return (_value & ((value_type)1 << bit)) != 0;
	}

	/// Resets a given bit to zero.
	void reset(int bit)
	{
		assert(bit >= 0 && bit <= (int)Bits);
		_value &= ~((value_type)1 << bit);
	}

	/// Resets the bits corresponding to the set bits in the given mask.
	void reset(const bitmask<T,Bits> &m)
	{
		_value &= ~ m.value();
	}

	/// Resets all bits to zero.
	void reset() { _value = 0; }

	/// Returns @c true if all bits are reset.
	bool empty() const { return _value == 0; }

	/// Tests whether the bitmask is empty.
	bool operator ! () const { return _value == 0; }

	/// Returns @c true if there is exactly one bit set.
	bool unique() const
	{
		return _value && (_value & (_value - 1)) == 0;
	}

	/// Returns @c true if there are no more than one bit set.
	bool empty_or_unique() const
	{
		return (_value & (_value - 1)) == T(0);
	}

	/// Returns the index of the least significant bit set.
	/// If no bit is set, returns @c -1.
	int smallest() const
	{
		return _value == 0 ? -1 : util::intrinsic::bit_scan_forward(_value);
	}

#if 0
	int count() const
	{
		int n = 0;
		for (int i = 0; i < MM_MAX_COLORS; ++i)
		{
			if (value & (1 << i))
				++n;
		}
		return n;
	}

	void set_count(int count)
	{
		assert(count >= 0 && count <= MM_MAX_COLORS);
		value = (1 << count) - 1;
	}
#endif

	/// Returns a bitmask with the least significant @c count bits set.
	static bitmask<T,Bits> fill(size_t count)
	{
		assert(count >= 0 && count <= Bits);
		return bitmask<T,Bits>(((value_type)1 << count) - 1);
	}
};

template <class T, size_t Bits>
inline bitmask<T,Bits> operator & (
	const bitmask<T,Bits> &x, 
	const bitmask<T,Bits> &y)
{
	return bitmask<T,Bits>(x.value() & y.value());
}

} // namespace util

#endif // UTILITIES_BITMASK_HPP
