///////////////////////////////////////////////////////////////////////////
// Codeword comparison routines.

#include <cassert>
#include <utility>
#include "util/simd.hpp"
#include "util/intrinsic.hpp"
#include "Algorithm.hpp"

//#include "util/call_counter.hpp"

/// Define the following macro to 1 to enable a specialized comparison routine
/// when one of the codewords (specifically, the secret) contains no repeated
/// colors. This marginally improves the performance at the cost of increased
/// code complexity and memory footprint. Therefore, it is not enabled by 
/// default.
#define SPECIALIZE_NONREPEATED_GENERIC 0

namespace Mastermind {

/// Codeword comparer for generic codewords (with or without repetition).
class GenericComparer
{
	// Lookup table that converts (nA<<4|nAB) -> feedback.
	// Both nA and nAB must be >= 0 and <= 15.
	struct lookup_table_t
	{
		Feedback table[0x100];

		lookup_table_t()
		{
			for (int i = 0; i < 0x100; i++)
			{
				int nA = i >> 4;
				int nAB = i & 0xF;
				table[i] = Feedback(nA, nAB - nA);
			}
		}
	};

	static const lookup_table_t lookup;

	typedef util::simd::simd_t<uint8_t,16> simd_t;

	simd_t secret;
	simd_t secret_colors;

public:

	GenericComparer(const Codeword &_secret) 
		: secret(*reinterpret_cast<const simd_t *>(&_secret))
	{
		// Change 0xff in secret to 0x0f.
		secret &= (uint8_t)0x0f;

		// Cache the color part of the secret.
		secret_colors = util::simd::keep_right<MM_MAX_COLORS>(secret);
	}

	/**
	 * In VC++, it's extremely subtle how the C++ code will compile into 
	 * assembly. Effectively equivalent code may lead to significantly 
	 * different performance. Usually the performance difference comes
	 * from redundant memory reads/writes. Hence it's always a good idea
	 * to examine the ASM after compiling the code.
	 *
	 * The following C++ code is hand-crafted so that the generated ASM
	 * by VC++ is almost optimal. The generated instructions for a loop
	 * that compares codewords and increment the frequencies are:
	 *
	 *   movdqa      xmm0,xmmword ptr [r10]  ; xmm0 := *guesses
	 *   add         r10,10h                 ; ++guesses
	 *   movdqa      xmm1,xmm0               
	 *   pminub      xmm0,xmm3              
	 *   pcmpeqb     xmm1,xmm4               
	 *   psadbw      xmm0,xmm5     
	 *   pand        xmm1,xmm2  
	 *   pextrw      edx,xmm0,4  
	 *   pextrw      eax,xmm0,0  
	 *   add         edx,eax  
	 *   psadbw      xmm1,xmm5  
	 *   pextrw      ecx,xmm1,4  
	 *   mov         ecx,ecx                 ; *
	 *   or          rdx,rcx  
	 *   movsx       rax,byte ptr [rdx+r11]  ; rax := lookup_table[mask]
	 *   inc         dword ptr [r9+rax*4]    ; ++freq[Feedback]
	 *   dec         r8                      ; --count
	 *   jne         ...                     ; continue loop
	 *
	 * The instruction marked with a (*) is redundant and would have been
	 * eliminated if the intrinsic function returns @c size_t instead of 
	 * @c int.
	 */
	Feedback operator () (const Codeword &_guess) const
	{
		// Create mask for use later.
		const simd_t mask_pegs = util::simd::fill_left<MM_MAX_PEGS>((uint8_t)0x10);
		// const simd_t mask_colors = util::simd::fill_right<MM_MAX_COLORS>((uint8_t)0xff);

		const simd_t guess(*reinterpret_cast<const simd_t *>(&_guess));

		// Compute nA.
		// It turns out that there's a significant (20%) performance hit if
		// <code>sum<8,16>(...)</code> is changed to <code>sum<0,16>(...)</code>
		// In the latter case, VC++ generates an extra (redundant) memory read.
		// Therefore we try to be save every instruction possible.
		unsigned int nA_shifted = util::simd::sum(util::simd::make_slice
			<MM_MAX_PEGS <= 8 ? 8 : 0, 16>((guess == secret) & mask_pegs));

		// Compute nAB.
		unsigned int nAB = util::simd::sum(util::simd::make_slice
			<0, MM_MAX_COLORS <= 8 ? 8 : 16>	(min(guess, secret_colors)));

		Feedback nAnB = lookup.table[nA_shifted|nAB];
		return nAnB;
	}
};

const GenericComparer::lookup_table_t GenericComparer::lookup;

/// Specialized codeword comparer for codewords without repetition.
class NoRepeatComparer
{
	// Pre-computed table that converts a comparison bitmask of
	// non-repeatable codewords into a feedback.
	//
	// Note that the performance impact of the table lookup can be subtle.
	// The VC++ compiler tends to generate a MOVZX (32-bit) or MOVSXD (64-bit)
	// instruction if the integer sizes do not match. What's worse, because
	// an SSE2 intrinsic function in use returns an <code>int</code>, there
	// is no way to get rid of the MOVSXD instruction which deteriates
	// performance by 30%!
	//
	// Hence, the best we could do is probably to compile under 32-bit, and 
	// access the table directly.
	struct lookup_table_t
	{
		Feedback table[0x10000];

		lookup_table_t()
		{
			for (int i = 0; i < 0x10000; i++)
			{
				int nA = util::intrinsic::pop_count((unsigned short)(i >> MM_MAX_COLORS));
				int nAB = util::intrinsic::pop_count((unsigned short)(i & ((1<<MM_MAX_COLORS)-1)));
				table[i] = Feedback(nA, nAB - nA);
			}
		}
	};

	static const lookup_table_t lookup;

	typedef util::simd::simd_t<int8_t,16> simd_t;
	
	simd_t secret;

public:

	NoRepeatComparer(const Codeword &_secret)
		: secret(*reinterpret_cast<const simd_t *>(&_secret))
	{
		// Change 0xFF in secret to 0x0F
		secret &= (int8_t)0x0f;

		// Set zero counters in secret to 0xFF, so that if a counter in the
		// guess and secret are both zero, they won't compare equal.
		secret |= util::simd::keep_right<MM_MAX_COLORS>(secret == simd_t::zero());
	}

	/**
	 * The following code translates into the following ASM:
	 *
	 *   movdqa      xmm0,xmmword ptr [rbx]  ; xmm0 := *guesses
	 *   movdqa      xmm1,xmm2               ; xmm1 := this->secret
	 *   lea         rbx,[rbx+10h]           ; ++guesses
	 *   pcmpeqb     xmm1,xmm0  
	 *   pmovmskb    eax,xmm1  
	 *   movsxd      rcx,eax                 ; *
	 *   movsx       rax,byte ptr [rcx+rdx]  ; rax := lookup_table[mask]
	 *   inc         dword ptr [rbp+rax*4]   ; ++freq[Feedback]
	 *   dec         rsi                     ; --count
	 *   jne         ...                     ; continue loop
	 *
	 * Note that we could save the instruction (*) if the intrinsic
	 * function for @c pmovmskb has return type int64 instead of int.
	 */
	Feedback operator () (const Codeword &guess) const 
	{
		return lookup.table
			[ util::simd::byte_mask(*reinterpret_cast<const simd_t *>(&guess) == secret) ];
	}
};

const NoRepeatComparer::lookup_table_t NoRepeatComparer::lookup;

/// Function object that appends a feedback to a feedback list.
class FeedbackUpdater
{
	Feedback * feedbacks;

public:
	
	explicit FeedbackUpdater(Feedback *fbs) : feedbacks(fbs) { }

	void operator () (const Feedback &fb) 
	{
		*(feedbacks++) = fb;
	}
};

/// Function object that increments the frequency statistic of a feedback.
class FrequencyUpdater
{
	unsigned int * freq;

public:
	
	// We do not zero the memory here. It must be initialized by the caller.
	explicit FrequencyUpdater(unsigned int *_freq) : freq(_freq) { }

	void operator () (const Feedback &fb) 
	{
		++freq[fb.value()];
	}
};

/// Function object that invokes two functions.
template <class T1, class T2>
class CompositeUpdater
{
	T1 u1;
	T2 u2;

public:

	CompositeUpdater(T1 updater1, T2 updater2) 
		: u1(updater1), u2(updater2) { }

	void operator () (const Feedback &fb)
	{
		u1(fb);
		u2(fb);
	}
};

/// Compares a secret to a list of codewords using @c Comparer, and 
/// processes each feedback using @c Updater.
template <class Comparer, class Updater>
static inline void compare_codewords(
	const Codeword &secret,
	const Codeword *_guesses,
	size_t _count,
	Updater _update)
{
	// The following redundant copy is to make VC++ happy treat the value
	// of @c _update as constant in the loop.
	Updater update(_update);

	Comparer compare(secret);
	size_t count = _count;
	const Codeword *guesses = _guesses;
	for (; count > 0; --count)
	{
		Feedback nAnB = compare(*guesses++);
		update(nAnB);
	}
}

/// Compares a secret to a list of codewords using @c Comparer, and 
/// update the feedback and/or frequencies.
template <class Comparer>
static inline void compare_codewords(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	Feedback *result,
	unsigned int *freq)
{
	if (freq && result)
	{
		FeedbackUpdater u1(result);
		FrequencyUpdater u2(freq);
		CompositeUpdater<FeedbackUpdater,FrequencyUpdater> update(u1,u2);
		compare_codewords<Comparer>(secret, guesses, count, update);
	}
	else if (freq)
	{
		FrequencyUpdater update(freq);
		compare_codewords<Comparer>(secret, guesses, count, update);
	}
	else if (result)
	{
		FeedbackUpdater update(result);
		compare_codewords<Comparer>(secret, guesses, count, update);
	}
}

#if 0
void compare_codewords(
	const Rules &rules,
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	unsigned int *freq)
{
	FrequencyUpdater update(freq);
	if (rules.repeatable())
		compare_codewords<GenericComparer>(secret, guesses, count, update);
	else
		compare_codewords<NoRepeatComparer>(secret, guesses, count, update);
}
#endif

#if SPECIALIZE_NONREPEATED_GENERIC
/// Specialized codeword comparer for generic codewords where the secret 
/// doesn't contain repeated colors.
class SpecializedComparer
{
	// Pre-computed table that converts a comparison bitmask into a feedback.
	struct lookup_table_t
	{
		Feedback table[0x10000];

		lookup_table_t()
		{
			for (int i = 0; i < 0x10000; i++)
			{
				int nA = util::intrinsic::pop_count((unsigned short)(i >> MM_MAX_COLORS));
				int nAB = util::intrinsic::pop_count((unsigned short)(i & ((1<<MM_MAX_COLORS)-1)));
				table[i] = Feedback(nA, nAB - nA);
			}
		}
	};

	static const lookup_table_t lookup;

	typedef util::simd::simd_t<uint8_t,16> simd_t;
	
	simd_t secret;
	simd_t barrier;
	static const simd_t _barrier;

public:

	SpecializedComparer(const Codeword &_secret)
		: secret(_secret.value()), barrier(_barrier)
	{
		// Change 0xFF in secret to 0x0F
		secret &= (uint8_t)0x0f;

		// Set zero counters in secret to 0xFF, so that if a counter in the
		// guess and secret are both zero, they won't compare equal.
		secret |= util::simd::keep_right<MM_MAX_COLORS>(secret == simd_t::zero());

		// Make a copy of barrier so that VC++ will not read the static memory
		// in a tight loop every time.
	}

	/**
	 * The following code translates into the following ASM:
	 * 
	 *   movdqa      xmm0,xmmword ptr [rdx]  ; xmm0 := *guesses
	 *   add         rdx,10h                 ; ++guesses
	 *   pminub      xmm0,xmm1               ; xmm0 := min(xmm0, barrier)
	 *   pcmpeqb     xmm0,xmm2               ; xmm0 := cmp(xmm0, secret)
	 *   pmovmskb    eax,xmm0                ; eax := movemask(xmm0)
	 *   movsxd      rcx,eax                 ; *
	 *   movsx       rax,byte ptr [rcx+r10]  ; feedback := map(eax)
	 *   inc         dword ptr [r11+rax*4]   ; ++freq[feedback]
	 *   dec         r8                      ; --count
	 *   jne         ...                     ; continue loop
	 *
	 * Note that we could save the instruction (*) if the intrinsic function 
	 * for @c pmovmskb had return type int64 instead of int.
	 */
	Feedback operator () (const Codeword &guess) const 
	{
		// Create a mask to convert color counters greater than one to one.
		// register simd_t barrier = _barrier;

		// Create a mask with the equal bytes in (augmented) secret and 
		// (augmented) guess set to 0xFF.
		simd_t mask = util::simd::min(simd_t(guess.value()), barrier) == secret;

		// Maps the mask to a feedback value and return.
		return lookup.table[ util::simd::byte_mask(mask) ];
	}
};

const SpecializedComparer::simd_t SpecializedComparer::_barrier = 
	util::simd::fill_left<MM_MAX_PEGS>(0xFF) | 
	util::simd::fill_right<MM_MAX_COLORS>(0x01);

const SpecializedComparer::lookup_table_t SpecializedComparer::lookup;
#endif

#if 0
static void compare_codewords_generic(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	Feedback *result,
	unsigned int *freq)
{
#if SPECIALIZE_NONREPEATED_GENERIC
	if (secret.has_repetition())
		return compare_codewords<GenericComparer>(secret, guesses, count, result, freq);
	else
		return compare_codewords<SpecializedComparer>(secret, guesses, count, result, freq);
#else
	return compare_codewords<GenericComparer>(secret, guesses, count, result, freq);
#endif
}

static void compare_codewords_norepeat(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	Feedback *result,
	unsigned int *freq)
{
	return compare_codewords<NoRepeatComparer>(secret, guesses, count, result, freq);
}
#endif

#if 0
/// [Test] Codeword comparer for generic codewords where the secret doesn't
/// contain repetitive colors.
class TestComparer
{
	// Pre-computed table that converts a comparison bitmask into a feedback.
	struct lookup_table_t
	{
		Feedback table[0x10000];

		lookup_table_t()
		{
			for (int i = 0; i < 0x10000; i++)
			{
				int nA = util::intrinsic::pop_count((unsigned short)(i >> MM_MAX_COLORS));
				int nAB = util::intrinsic::pop_count((unsigned short)(i & ((1<<MM_MAX_COLORS)-1)));
				table[i] = Feedback(nA, nAB - nA);
			}
		}
	};

	static const lookup_table_t lookup;

	typedef util::simd::simd_t<uint8_t,16> simd_t;
	
	simd_t secret;

public:

	TestComparer(const Codeword &_secret)
		: secret(_secret.value()) 
	{
		// Change 0xFF in secret to 0x0F
		secret &= (uint8_t)0x0f;

		// Set zero counters in secret to 0xFF, so that if a counter in the
		// guess and secret are both zero, they won't compare equal.
		secret |= util::simd::keep_right<MM_MAX_COLORS>(secret == simd_t::zero());
	}

	/**
	 * The following code translates into the following ASM:
	 * 
	 *   movdqa      xmm0,xmmword ptr [rdx]  ; xmm0 := *guesses
	 *   add         rdx,10h                 ; ++guesses
	 *   pminub      xmm0,xmm1               ; xmm0 := min(xmm0, barrier)
	 *   pcmpeqb     xmm0,xmm2               ; xmm0 := cmp(xmm0, secret)
	 *   pmovmskb    eax,xmm0                ; eax := movemask(xmm0)
	 *   movsxd      rcx,eax                 ; *
	 *   movsx       rax,byte ptr [rcx+r10]  ; feedback := map(eax)
	 *   inc         dword ptr [r11+rax*4]   ; ++freq[feedback]
	 *   dec         r8                      ; --count
	 *   jne         ...                     ; continue loop
	 *
	 * Note that we could save the instruction (*) if the intrinsic function 
	 * for @c pmovmskb had return type int64 instead of int.
	 */
	Feedback operator () (const Codeword &guess) const 
	{
		// Create a mask to convert color counters greater than one to one.
		const simd_t barrier = 
			util::simd::fill_left<MM_MAX_COLORS>(0xFF) |
			util::simd::fill_right<MM_MAX_COLORS>(0x01);

		// Create a mask with the equal bytes in (augmented) secret and 
		// (augmented) guess set to 0xFF.
		simd_t mask = util::simd::min(simd_t(guess.value()), barrier) == secret;

		// Maps the mask to a feedback value and return.
		return lookup.table[ util::simd::byte_mask(mask) ];
	}
};

const TestComparer::lookup_table_t TestComparer::lookup;

static void compare_codewords_test(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	Feedback *result,
	unsigned int *freq)
{
	if (secret.has_repetition())
		return compare_codewords<GenericComparer>(secret, guesses, count, result, freq);
	else
		return compare_codewords<TestComparer>(secret, guesses, count, result, freq);
}

REGISTER_ROUTINE(ComparisonRoutine, "test", compare_codewords_test)
#endif

/// Compares generic codewords and returns feedbacks.
void CompareGeneric1(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	Feedback *result)
{
	FeedbackUpdater update(result);
	compare_codewords<GenericComparer>(secret, guesses, count, update);
}

/// Compares generic codewords and returns frequencies.
void CompareGeneric2(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	unsigned int *freq)
{
	FrequencyUpdater update(freq);
	compare_codewords<GenericComparer>(secret, guesses, count, update);
}

/// Compares generic codewords and returns feedbacks and frequencies.
void CompareGeneric3(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	Feedback *result,
	unsigned int *freq)
{
	FeedbackUpdater u1(result);
	FrequencyUpdater u2(freq);
	CompositeUpdater<FeedbackUpdater,FrequencyUpdater> update(u1,u2);
	compare_codewords<GenericComparer>(secret, guesses, count, update);
}

/// Compares norepeat codewords and returns feedbacks.
void CompareNorepeat1(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	Feedback *result)
{
	FeedbackUpdater update(result);
	compare_codewords<NoRepeatComparer>(secret, guesses, count, update);
}

/// Compares norepeat codewords and returns frequencies.
void CompareNorepeat2(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	unsigned int *freq)
{
	FrequencyUpdater update(freq);
	compare_codewords<NoRepeatComparer>(secret, guesses, count, update);
}

/// Compares norepeat codewords and returns feedbacks and frequencies.
void CompareNorepeat3(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	Feedback *result,
	unsigned int *freq)
{
	FeedbackUpdater u1(result);
	FrequencyUpdater u2(freq);
	CompositeUpdater<FeedbackUpdater,FrequencyUpdater> update(u1,u2);
	compare_codewords<NoRepeatComparer>(secret, guesses, count, update);
}

} // namespace Mastermind
