#include <cassert>
#include <functional>

#include "Engine.hpp"
#include "ObviousStrategy.hpp"

namespace Mastermind {

/**
 * Checks whether two codewords contain exactly the same set of colors.
 */
static bool contain_same_colors(const Codeword &a, const Codeword &b)
{
	for (int c = 0; c < MM_MAX_COLORS; ++c)
	{
		if (a.count(c) != b.count(c))
			return false;
	}
	return true;
}

/**
 * Estimates an obvious lower-bound of the cost of making a non-possible
 * guess against a set of remaining possibilities. If such a lower-bound
 * is found, returns the minimum steps, depth, and worst count. If no
 * lower-bound can be found easily, returns zero.
 */
StrategyCost estimate_obvious_lowerbound(
	const Rules &rules,
	CodewordConstRange possibilities)
{
	// We are only concerned with a handful of remaining possibilities.
	// If there are too many, we won't make an attempt.
	int p = rules.pegs();
	int n = (int)possibilities.size();
	if (n > p*(p+3)/2)
		return StrategyCost();

	// Partition the possibilities by the colors they contain. Only the
	// size of the partitions matter; the particular colors don't.
	// @todo Since we will only be processing with MM_MAX_PEGS+1 groups,
	// we can return -1 if there are more secret groups.
	const int M = MM_MAX_PEGS * (MM_MAX_PEGS + 3) / 2;
	bool visited[M] = { false };
	int group[M] = { 0 };
	int ngroup = 0;
	for (int i = 0; i < n; ++i)
	{
		if (!visited[i])
		{
			for (int j = i + 1; j < n; ++j)
			{
				if (!visited[j] && contain_same_colors(possibilities[i], possibilities[j]))
				{
					++group[ngroup];
					visited[j] = true;
				}
			}
			++ngroup;
		}
	}

	// Now we have classified the remaining possibilities into ngroup groups
	// according to the colors they contain. For any given guess, the secrets
	// in the same group must have the same number of common colors with
	// this guess.

	// We next classify the possible feedback values by the number of common
	// colors, like below:
	//   0: 0A0B
	//   1: 0A1B 1A0B
	//   2: 0A2B 1A1B 2A0B
	//   3: 0A3B 1A2B 2A1B 3A0B
	//   4: 0A4B 1A3B 2B2B - -
	// Note that 3A1B is impossible and 4A0B is excluded because the guess is
	// assumed to be outside the remaining possibilities.

	// Next, we assign each secret group a feedback group. Note that multiple
	// secret groups may be assigned the same feedback group, but two secrets
	// in the same secret group may not be split into different feedback groups.
	// We try to find a simple assignment that is guaranteed to minimize the
	// total number of steps needed to reveal all secrets. If a simple
	// assignment cannot be found easily, we give up and return failure.

	// The algorithm is to assign each secret group, in order of decreasing
	// size, the largest remaining feedback group if this feedback group is
	// no larger than the secret group. If this can be done, the resulting
	// assignment is guaranteed (??? proof needed) to yield a lower bound
	// of the total number of steps.
	if (ngroup > p + 1)
		return StrategyCost();

	std::sort(group + 0, group + ngroup, std::greater<int>());
	int extra = 0;
	for (int i = 0; i < ngroup; ++i)
	{
		int avail = p - i - ((i < 2)? 1 : 0);
		if (group[i] >= avail)
		{
			extra += (group[i] - avail);
		}
		else
		{
			return StrategyCost();
		}
	}
	return StrategyCost(
		extra + n * 2,          // total steps
		extra > 0 ? 3 : 2,      // max depth
		(unsigned short)(extra > 0 ? extra : n)); // worst count
}

/**
 * Returns an obviously-optimal guess if one exists.
 * We will only make a guess from the list of remaining possibilities.
 *
 * @param e Algorithm engine.
 * @param possibilities List of remaining possibilities.
 * @param max_depth Maximum number of guesses allowed, including the
 *      initial guess.
 * @param min_obj Minimum objective that an obvious guess must meet.
 * @param cost If an obvious guess is found, stores the cost of an optimal
 *      strategy that begins with this guess.
 * @param obj If an obvious guess is found, stores in what sense is this
 *      guess optimal.
 *
 * @returns An obvious guess such that an optimal strategy starting with
 *      this guess has the smallest number of total steps, will not take
 *      more than max_depth depth, and achieves at least min_obj.
 */
Codeword make_obvious_guess(
	const Engine *e,
	CodewordConstRange possibilities,
	int max_depth,
	StrategyObjective min_obj,
	StrategyCost &cost,
	StrategyObjective &obj)
{
	int count = (int)possibilities.size();

	// If there are no remaining possibilities, no guess is needed.
	if (count == 0)
		return Codeword();

	// Now we need to make at least one guess. If that's not allowed, return.
	if (max_depth < 1)
		return Codeword();

	// If there is only one possibility left, guess it.
	if (count == 1)
	{
		cost = StrategyCost(1, 1, 1);
		obj = MinWorst;
		return possibilities[0];
	}

	// Now we need to make at least two guesses to reveal every secret.
	// If that's not allowed, return.
	if (max_depth < 2)
		return Codeword();

	// If there are only two possibilities left, guess the first one.
	if (count == 2)
	{
		cost = StrategyCost(3, 2, 1);
		obj = MinWorst;
		return possibilities[0];
	}

	// If the number of possibilities is greater than the number of
	// distinct feedbacks, there will (very likely) be no obvious guess.
	// So we won't make any further attempt.
	int p = e->rules().pegs();
	if (count > p*(p+3)/2)
		return Codeword();

	// Check each remaining possiblity as guess in turn.
	// - If a guess partitions the remaining possibilities into singleton
	//   cells, return it immediately.
	// - If no such guess exists but one exists that partitions them into
	//   (n-1) singleton cells and one cell with 2 secrets, return it.
	//   Note however that such a guess may only be optimal in terms of
	//   total steps; a non-possible guess may yield lower depth (i.e. 1).
	// - Otherwise, keep the best guess that partitions the remaining
	//   possibilities into cells with no more than two secrets. Compare
	//   this cost with an estimated obvious lower bound of a non-possible
	//   guess. If the cost < lower bound, return this guess firmly. If
	//   cost = lower bound, then this guess is only optimal in terms of
	//   total number of steps, but not in depth.
	Codeword best_guess;
	int best_extra = -1;
	for (int i = 0; i < count; ++i)
	{
		Codeword guess = possibilities[i];
		FeedbackFrequencyTable freq = e->compare(guess, possibilities);
		unsigned int nonzero = 1;  // 4A0B
		unsigned int maxfreq = 0;
		for (size_t j = 0; j < freq.size() - 2; ++j) // skip 3A0B and 4A0B
		{
			if (freq[j])
			{
				maxfreq = std::max(maxfreq, freq[j]);
				++nonzero;
			}
		}

		if (maxfreq == 1) // all cells are singleton cells
		{
			cost = StrategyCost(2*count-1, 2, (unsigned short)(count-1));
			obj = MinWorst;
			return guess;
		}
		if (maxfreq > 2)
		{
			continue;
		}

		int extra = (int)(count - nonzero); // number of cells with 2 secrets
		if (best_extra < 0 || extra < best_extra)
		{
			best_extra = extra;
			best_guess = guess;
		}
	}

	// If no guess partitions the remaining possibilities into cells with
	// no more than two secrets, there is no obvious strategy.
	if (best_extra < 0)
		return Codeword();

	// Update the cost.
	cost = StrategyCost(2*count-1+best_extra, 3, 1);

	// If exactly one cell contains two secrets and all the rest are
	// singleton, then this guess is guaranteed to be optimal in steps,
	// and we needn't check further if that's what's required.
	if (best_extra == 1 && min_obj == MinSteps)
	{
		obj = MinSteps;
		return best_guess;
	}

	// Now it's not obvious whether our guess is optimal. We make an attempt
	// to estimate the lower bound of cost of any non-possible guess.
#if 1
	StrategyCost lb = estimate_obvious_lowerbound(e->rules(), possibilities);
	if (!superior(lb, cost, min_obj))
	{
		obj = min_obj;
		return best_guess;
	}
#endif
	return Codeword();
}

#if 0
// Searches for an obviously optimal strategy.
// If one exists, returns the total cost of the strategy.
// Otherwise, returns -1.
// @todo Change return type to StrategyCost.
static int fill_obviously_optimal_strategy(
	Engine &e,
	CodewordRange secrets,
	bool min_depth,    // whether to minimize the worst-case depth
	int max_depth,     // maximum number of extra guesses, not counting
	                   // the initial guess
	StrategyTree &tree // Strategy tree that stores the best strategy
	)
{
	int extra;
	//Codeword guess = ObviousStrategy(e).make_guess(secrets, &extra);
	Codeword guess = make_less_obvious_guess(e, secrets, &extra);
	--extra;
	// @todo Take into account the lower-bound estimate of non-possible
	// guesses.
	if (!guess || (extra > max_depth) || (min_depth && extra > 1))
		return -1;

	//	VERBOSE_COUT << "Found obvious guess: " << obvious << std::endl;

	// Automatically fill the strategy tree using this guess.This requires
	// all cells in the partition to have no more than two possibilities.
	// This is equivalent to Knuth's 'x' notation in writing a strategy.
	Feedback perfect = Feedback::perfectValue(e.rules());
	FeedbackList fbs;
	e.compare(guess, secrets, fbs);
	size_t n = secrets.size();
	int depth = tree.currentDepth();
	int cost = 0;

	for (size_t j = 0; j < Feedback::size(e.rules()); ++j)
	{
		Codeword first;
		for (size_t i = 0; i < n; ++i)
		{
			if (fbs[i] == Feedback(j))
			{
				if (!first)
				{
					++cost;
					first = secrets[i];
					StrategyTree::Node node(depth + 1, guess, fbs[i]);
					tree.append(node);
					if (fbs[i] != perfect)
					{
						++cost;
						StrategyTree::Node leaf(depth + 2, first, perfect);
						tree.append(leaf);
					}
				}
				else
				{
					cost += 3;
					tree.append(StrategyTree::Node(depth + 2,
						first, e.compare(secrets[i], first)));
					tree.append(StrategyTree::Node(depth + 3,
						secrets[i], perfect));
				}
			}
		}
	}
	return cost;
	// @todo: try guesses from non-secrets.
}
#endif

} // namespace Mastermind
