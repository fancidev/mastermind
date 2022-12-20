#include "Codeword.hpp"
#include "Algorithm.hpp"
#include "util/simd.hpp"

namespace Mastermind {

// Routine to compute the color mask of an array of codewords.

// This implementation works with "long" codeword.
// The performance is slightly (5-10%) worse than v2,
// but the implementation is much more readable!
// Plus, this routine is not the bottleneck!
// So we'll use v1 in production.
static unsigned short ScanDigitMask_v1(
	const Codeword *_first,
	const Codeword *_last)
{
	typedef util::simd::simd_t<uint8_t,16> simd_t;

	simd_t mask = simd_t::zero();
	const simd_t *first = (const simd_t *)_first;
	const simd_t *last = (const simd_t *)_last;
	for (const simd_t *it = first; it != last; ++it)
	{
		mask |= *it;
	}

	mask = (mask == simd_t::zero());
	unsigned short result = (unsigned short)util::simd::byte_mask(mask);
	result = (~result) & ((1 << MM_MAX_COLORS) - 1);
	return result;
}

///////////////////////////////////////////////////////////////////////////
// Interface Routines
//

unsigned short GetPresentColors(const Codeword *first, const Codeword *last)
{
	return ScanDigitMask_v1(first, last);
}

} // namespace Mastermind
