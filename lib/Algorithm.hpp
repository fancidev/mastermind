#ifndef MASTERMIND_ALGORITHM_HPP
#define MASTERMIND_ALGORITHM_HPP

#include <string>
#include "Rules.hpp"
#include "Codeword.hpp"
#include "Feedback.hpp"
#include "Registry.hpp"

namespace Mastermind {

/// Type of a function that compares a codeword to a list of codewords and
/// returns the feedbacks.
typedef void ComparisonRoutine1(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	Feedback *result);

/// Type of a function that compares a codeword to a list of codewords and
/// increments the feedback frequencies. The caller is responsible for 
/// allocating the frequencies and initializing them to zero.
typedef void ComparisonRoutine2(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	unsigned int *freq);

/// Type of a function that compares a codeword to a list of codewords and
/// returns the feedbacks as well as increments the feedback frequencies. 
/// The caller is responsible for allocating the frequencies and initializing
/// them to zero.
typedef void ComparisonRoutine3(
	const Codeword &secret,
	const Codeword *guesses,
	size_t count,
	Feedback *result,
	unsigned int *freq);

/// Comparison functions for generic codewords.
extern ComparisonRoutine1 CompareGeneric1;
extern ComparisonRoutine2 CompareGeneric2;
extern ComparisonRoutine3 CompareGeneric3;

/// Comparison functions for norepeat codewords.
extern ComparisonRoutine1 CompareNorepeat1;
extern ComparisonRoutine2 CompareNorepeat2;
extern ComparisonRoutine3 CompareNorepeat3;

/// Generates all codewords conforming to the given set of rules. 
/// The caller is responsible for allocating memory for the results.
extern void GenerateCodewords(const Rules &rules, Codeword *results);

/// <summary>
/// Gets a bit-mask of the colors present in a list of codewords.
/// </summary>
/// <param name="first">Pointer to the first element of the codeword list.</param>
/// <param name="last">Pointer to one past the last element of the codeword list.</param>
/// <remarks>
/// In the bit-mask returned, a bit is set if the corresponding color 
/// is present in at least one of the codewords. A bit is cleared if 
/// the corresponding color never appears in any of the codewords. 
/// The bits are numbered from LSB to MSB. 
///
/// For example, if the codeword is <code>4169</code>, then bits 1, 4, 
/// 6, and 9 are set, and the rest are unset. The mask returned will  
/// be <code>0000-0010-0101-0010</code>, or <code>0x0252</code>.
///
/// Note that the highest bit (corresponding to color 0xF) is never set
/// in the returned mask, even if it exists in the codeword, because
/// <code>0xF</code> is reserved for special use.
/// </remarks>
extern unsigned short GetPresentColors(const Codeword *first, const Codeword *last);

template <class Routine>
struct RoutineRegistry : public Utilities::Registry<std::string, Routine>
{
};

#define REGISTER_ROUTINE(type,id,item) REGISTER_ITEM(type,id,item)

} // namespace Mastermind

#endif // MASTERMIND_ALGORITHM_HPP
