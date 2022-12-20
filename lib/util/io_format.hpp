/// @defgroup IOFormat IO Stream Formatter Helper
/// @ingroup util

#ifndef UTILITIES_IO_FORMAT_HPP
#define UTILITIES_IO_FORMAT_HPP

#include <iostream>

namespace util {

/// Defines the file format to serialize an object.
/// @ingroup IOFormat
enum FileFormat
{
	DefaultFormat = 0,
	TextFormat = 1,
	XmlFormat = 2,
	BinaryFormat = 3,
	CompressedFormat = 4
};

#if 0
class withheader
{
	bool _value;

	static int Index() 
	{ 
		static int i = std::ios_base::xalloc();
		return i;
	}

public:
	withheader() : _value(true) { }
	withheader(bool value) : _value(value) { }

	friend std::ostream& operator << (std::ostream &os, const withheader
};
#endif

namespace details {

/// Returns the index to a custom ios format field.
inline int header_index() 
{
	static int i = std::ios_base::xalloc();
	return i;
}

} // namespace details

/// Indicates that a header should be (de)serialized.
/// @ingroup IOFormat
inline std::ostream& header(std::ostream &os)
{
	os.iword(details::header_index()) = 1;
	return os;
}

/// Indicates that no header should be (de)serialized.
/// @ingroup IOFormat
inline std::ostream& noheader(std::ostream &os)
{
	os.iword(details::header_index()) = 0;
	return os;
}

/// Checks whether a header should be (de)serialized.
/// By default, no header is (de)serialized.
/// @ingroup IOFormat
inline bool hasheader(std::ostream &os)
{
	return os.iword(details::header_index()) != 0;
}

} // namespace util

#endif // UTILITIES_IO_FORMAT_HPP
