#ifndef MASTERMIND_OBVIOUS_STRATEGY_HPP
#define MASTERMIND_OBVIOUS_STRATEGY_HPP

#include <cassert>
#include "Engine.hpp"
#include "Strategy.hpp"

namespace Mastermind {

Codeword make_obvious_guess(
	const Engine *e,
	CodewordConstRange possibilities,
	int max_depth,
	StrategyObjective min_obj,
	StrategyCost &cost,
	StrategyObjective &obj);

/**
 * Strategy that makes an obviously-optimal guess when one exists.
 *
 * This is an example (but useful) strategy that implements the @c
 * Strategy interface. This strategy returns an obviously-optimal
 * guess if one exists. If such a guess is not found, it returns
 * failure. Therefore, this strategy cannot be used alone to break
 * a code; instead it serves as a quick check before other more
 * sophisticated (and time consuming) strategies are applied.
 *
 * An obviously-optimal guess is a guess (which may or may not be
 * in the remaining possibilities) that partitions the possibility
 * set into discrete cells, i.e. where every cell contains exactly
 * one element. If such a guess exists and comes from the possibility
 * set, then it is optimal because it reveals one secret in the
 * immediate step and reveals all the rest in two steps. If no such
 * guess exists in the possibility set but one exists outside the
 * possibility set, then that one is optimal because it reverals
 * all secrets in two steps.
 *
 * A necessary condition for an obviously-optimal guess to exist
 * is that the number of remaining possibilities is no more than
 * the number of distinct feedbacks. For a game with @c p pegs,
 * the number of distinct feedbacks is <code>p*(p+3)/2</code>.
 *
 * If the maximum depth (i.e. number of guesses) to reveal each secret
 * is not restricted, an additional observation can be made. If there
 * exists no guess in the remaining possibilities that can partition
 * @c n remaining possibilities into @c n singleton cells, but there
 * exists one guess in the remaining possibilities that partitions the
 * remaining possibilities into <code>(n-2)</code> singleton cells
 * and one cell with 2 secrets, then this guess achieves the minimum
 * number of total steps, which is equal to @c n. However, it would
 * require two extra guesses in the worst case; there might exist a
 * non-possible guess that partitions the remaining possibilities into
 * singleton cells and thus requires only one extra guess.
 *
 * @ingroup Obvious
 */
class ObviousStrategy : public Strategy
{
	const Engine *e;

public:

	/// Constructs the strategy.
	ObviousStrategy(const Engine *engine) : e(engine) { }

	virtual std::string name() const
	{
		return "obvious";
	}

	virtual std::string description() const
	{
		return "makes an obviously optimal guess";
	}

	/*
	 * Returns an obviously-optimal guess if one exists.
	 *
	 * @param possibilities List of remaining possibilities.
	 * @param max_depth If an obvious guess is found, stores on return
	 *      the maximum number of guesses required to reveal all secrets,
	 *      including the first guess.
	 * @timecomplexity No more than <code>K=p*(p+3)/2</code> passes.
	 *      In each pass, a candidate guess is compared to all @c N
	 *      possibilities, and the frequencies of the feedbacks are
	 *      checked. The overall complexity is <code>O(KN)</code>.
	 *
	 * @spacecomplexity Constant.
	 */
	Codeword make_guess(
		CodewordConstRange possibilities,
		int *max_depth) const
	{
		assert(max_depth != NULL);
		size_t count = possibilities.size();

		// If there are no remaining possibilities, no guess is needed.
		if (count == 0)
		{
			*max_depth = 0;
			return Codeword();
		}

		// If there is only one possibility left, guess it.
		if (count == 1)
		{
			*max_depth = 1;
			return *possibilities.begin();
		}

		// If there are only two possibilities left, guess the first one.
		if (count == 2)
		{
			*max_depth = 2;
			return *possibilities.begin();
		}

		// If the number of possibilities is greater than the number of
		// distinct feedbacks, there will be no obvious guess.
		size_t p = e->rules().pegs();
		if (count > p*(p+3)/2)
			return Codeword();

		// Returns the first obviously optimal guess in the possibility
		// set (if any). If a less-obviously optimal guess is found,
		// store it temporarily.
		Codeword less_obvious_guess;
		for (size_t i = 0; i < count; ++i)
		{
			Codeword guess = possibilities[i];
			FeedbackFrequencyTable freq = e->compare(guess, possibilities);
			size_t nonzero = freq.nonzero_count();
			if (nonzero == count)
			{
				*max_depth = 2;
				return guess;
			}
			if (nonzero == count - 1 && less_obvious_guess.IsEmpty())
			{
				less_obvious_guess = guess;
			}
		}

		// Returns the less obvious guess if one exists.
		if (!less_obvious_guess.IsEmpty())
		{
			*max_depth = 3;
		}
		return less_obvious_guess;
	}

	/*
	 * Returns an obviously-optimal guess if one exists.
	 *
	 * @param possibilities List of remaining possibilities.
	 * @param candidates    List of (suggestive) candidate guesses.
	 * @returns An obviously-optimal guess if one is found; otherwise
	 *      <code>Codeword::emptyValue()</code>.
	 * @timecomplexity No more than <code>K=p*(p+3)/2</code> passes.
	 *      In each pass, a candidate guess is compared to all @c N
	 *      possibilities, and the frequencies of the feedbacks are
	 *      checked. The overall complexity is <code>O(KN)</code>.
	 * @spacecomplexity Constant.
	 * @remarks The possibility set MUST be a subset of candidates.
	 *
	 * @todo Currently we only support candidate == possibilities.
	 *      This should be improved in the future.
	 */
	virtual Codeword make_guess(
		CodewordConstRange possibilities,
		CodewordConstRange /* candidates */) const
	{
		int max_depth = -1;
		Codeword guess = make_guess(possibilities, &max_depth);
		if (!guess.IsEmpty() && max_depth <= 2)
			return guess;
		else
			return Codeword();
	}
};

} // namespace Mastermind

#endif // MASTERMIND_OBVIOUS_STRATEGY_HPP
