#ifndef MASTERMIND_FILTER_HPP
#define MASTERMIND_FILTER_HPP

#include <vector>
#include <bitset>
#include "Engine.hpp"

namespace Mastermind {

/**
 * Represents a constraint of the form <code>(guess,response)</code>.
 */
#if 0
struct Constraint
{
	Codeword guess;
	Feedback response;
};
#endif

/// Returns a subset of the codewords consisting of elements satisfy
/// the given filter.
// not used
//CodewordList match(CodewordConstRange codewords, const Filter &filter);

#if 0
/**
 * Represents a state of the game. 
 *
 * Conceptually, a <i>state</i> of the game in an intermediate status
 * of the guessing process. A state is uniquely defined by a filter. 
 * However, the @c State structure contains some more book-keeping
 * members.
 */
struct State
{
	/// Set of remaining possible secrets given the feedbacks so far.
	/// Since this is a pair of iterators, the actual elements must
	/// be allocated separately.
	CodewordRange possibilities;

	/// Bit-mask of colors that have not been guessed so far.
	ColorMask fresh_colors;

	/// Bit-mask of colors that are excluded from the possibilities.
	/// These colors are sure to be not present in the secret.
	ColorMask excluded_colors;

	/// Initializes the state at the beginning of a game.
	State(Engine &e, CodewordRange all)
		: possibilities(all),
		fresh_colors((1 << e.rules().colors()) - 1), 
		excluded_colors(0) { }
	
	/// Updates the status.
	void udpate(	Engine &e, 
		const Codeword &guess, Feedback /*feedback*/, CodewordRange remaining)
	{
		ColorMask allmask((1 << e.rules().colors()) - 1);
		fresh_colors &= ~ e.colorMask(guess);
		excluded_colors = allmask & ~ e.colorMask(remaining);
		possibilities = remaining;
	}

};
#endif

} // namespace Mastermind

#endif // MASTERMIND_FILTER_HPP
