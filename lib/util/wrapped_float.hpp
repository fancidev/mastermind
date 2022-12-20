/// @defgroup WrappedFloat Wrapped Floating Number
/// @ingroup util

#ifndef UTILITIES_WRAPPED_FLOAT_HPP
#define UTILITIES_WRAPPED_FLOAT_HPP

#include <limits>

namespace util {

/**
 * Represents a wrapped floating number suitable for inexact comparison.
 * @ingroup WrappedFloat
 */
template <class T, unsigned int NEps>
class wrapped_float
{
	T _value;

public:
	
	/// Wraps the floating value <code>0.0</code>.
	wrapped_float() : _value(0) { }

	/// Wraps the given floating value.
	explicit wrapped_float(const T &value) : _value(value) { }

	/// Returns the exact value being wrapped.
	T value() const { return _value; }

	/// Implicitly convert the wrapped value to its exact value.
	operator T () const { return _value; }
};

/// Compares two wrapped floating numbers.
/// @ingroup WrappedFloat
template <class T, unsigned int NEps>
bool operator < (const wrapped_float<T,NEps> &x, const wrapped_float<T,NEps> &y)
{
	return x.value() < y.value() - NEps * std::numeric_limits<T>::epsilon();
}

/// Compares two wrapped floating numbers.
/// @ingroup WrappedFloat
template <class T, unsigned int NEps>
bool operator > (const wrapped_float<T,NEps> &x, const wrapped_float<T,NEps> &y)
{
	return (y < x);
}

/// Compares two wrapped floating numbers.
/// @ingroup WrappedFloat
template <class T, unsigned int NEps>
bool operator <= (const wrapped_float<T,NEps> &x, const wrapped_float<T,NEps> &y)
{
	return !(y < x);
}

/// Compares two wrapped floating numbers.
/// @ingroup WrappedFloat
template <class T, unsigned int NEps>
bool operator >= (const wrapped_float<T,NEps> &x, const wrapped_float<T,NEps> &y)
{
	return !(x < y);
}

/// Compares two wrapped floating numbers.
/// @ingroup WrappedFloat
template <class T, unsigned int NEps>
bool operator == (const wrapped_float<T,NEps> &x, const wrapped_float<T,NEps> &y)
{
	return !(x < y) && !(y < x);
}

/// Compares two wrapped floating numbers.
/// @ingroup WrappedFloat
template <class T, unsigned int NEps>
bool operator != (const wrapped_float<T,NEps> &x, const wrapped_float<T,NEps> &y)
{
	return (x < y) || (y < x);
}

} // namespace util

#endif // UTILITIES_WRAPPED_FLOAT_HPP
