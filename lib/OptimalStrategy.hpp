#ifndef MASTERMIND_OPTIMAL_STRATEGY_HPP
#define MASTERMIND_OPTIMAL_STRATEGY_HPP

#include <cassert>
#include <vector>
#include <numeric>

#include "Engine.hpp"
#include "Strategy.hpp"
#include "util/call_counter.hpp"
#include "util/intrinsic.hpp"

namespace Mastermind {

/**
 * Define MINIMIZE_WORSTCASE_COUNT = 1 to include code to minimize the 
 * number of secrets that are revealed by the maximum number of steps.
 * This will produce a better-looking result, but will take more time
 * to run.
 */
#ifndef MINIMIZE_WORSTCASE_COUNT
#define MINIMIZE_WORSTCASE_COUNT 0
#endif

namespace Heuristics {

/**
 * Special-purpose heuristic used by an optimal strategy to score a
 * candidate guess by the lower bound of the cost if this guess is
 * made. This heuristic could also be used by a heuristic strategy.
 *
 * @ingroup Optimal
 * @todo Improve the lower-bound estimate.
 */
class MinimizeLowerBound
{
public:

	/// Type of the heuristic score.
	typedef StrategyCost score_t;

	/// Returns a simple estimate of minimum total number of steps
	/// required to reveal @c n secrets given a branching factor of @c b,
	/// including the initial guess.
	static score_t simple_estimate(
		int n, // Number of remaining secrets
		int b  // Branching factor: number of distinct non-perfect feedbacks
		)
	{
		score_t cost;
		for (int remaining = n, count = 1; remaining > 0; )
		{
			cost.steps += remaining;
			cost.depth++;
#if MINIMIZE_WORSTCASE_COUNT
			cost.worst2 = cost.worst1;
			cost.worst1 = std::min(count, remaining);
#endif
			remaining -= count;
			count *= b;
		}
		return cost;
	}

private:

	//Engine &e;
	std::vector<score_t> _cache;

public:

	/// Constructs the heuristic.
	MinimizeLowerBound(const Engine *engine)
		: /* e(engine), */ _cache(engine->rules().size()+1)
	{
		// Build a cache of simple estimates.
		int p = engine->rules().pegs();
		int b = p*(p+3)/2-1;
		for (size_t n = 0; n < _cache.size(); ++n)
		{
			_cache[n] = simple_estimate((int)n, b);
		}
	}

	/// Returns the name of the heuristic.
	std::string name() const { return "Min-LB"; }

	/// Returns a simple estimate of minimum total number of steps
	/// required to reveal @c n secrets, including the initial guess,
	/// assuming the maximum branching factor for the game.
	score_t simple_estimate(int n) const
	{
		assert(n >= 0 && n < (int)_cache.size());
		return _cache[n];
	}

	/// Computes the heuristic score. The score consist of two parts:
	/// - The total number of steps needed to reveal all secrets, excluding
	///   the initial guess
	/// - The maximum depth (i.e. number of extra guesses) needed to reveal
	///   every secret.
	score_t compute(const FeedbackFrequencyTable &freq) const
	{
		// Note: we make the critical assumptions that:
		// - Feedback::size()-1 is the perfect feedback, and
		// - Feedback::size()-2 is the impossible (p-1,1) feedback
		// Both cases can be ignored for the purpose of computing this
		// heuristic score. This reduces branching in the tight loop
		// and hence improves performance.

		// In addition, in order to find the maximum depth of all paritions,
		// we use a bitset for quicker operation.
		// Note: We could further save two instructions by storing tmp.depth
		// as a bit-mask directly. However, this has two drawbacks:
		// 1) it obscures the source code;
		// 2) it makes it harder to implement 'worst' optimization in the
		//    future.
		// Therefore, we will not use that approach for now.
		unsigned int depth_bitset = 0;
		int steps = 0;
		size_t m = freq.size() - 2;
		for (size_t j = 0; j < m; ++j)
		{
#if 0
			if (freq[j] == 0) 
				continue;
#endif
			score_t tmp = simple_estimate(freq[j]);
			steps += tmp.steps;
#if MINIMIZE_WORSTCASE_COUNT
			if (tmp.depth > lb.depth)
			{
				lb.depth = tmp.depth;
				lb.worst1 = tmp.worst1;
				lb.worst2 = tmp.worst2;
			}
			else if (tmp.depth == lb.depth)
			{
				lb.worst1 += tmp.worst1;
				lb.worst2 += tmp.worst2;
			}
			else if (tmp.depth == lb.depth - 1)
			{
				lb.worst2 += tmp.worst1;
			}
#else
			//lb.depth = std::max(lb.depth, tmp.depth);
			depth_bitset |= (1 << tmp.depth);
#endif
		}
		score_t lb;
		lb.steps = steps;
		lb.depth = (unsigned short)(depth_bitset == 0 ? 
			0 : util::intrinsic::bit_scan_reverse(depth_bitset));

		UPDATE_CALL_COUNTER("ComputeLowerBound_Steps", lb.steps);
		UPDATE_CALL_COUNTER("ComputeLowerBound_Depth", lb.depth);
		return lb;
	}

#if 0
	// Computes a lower bound for each candidate guess.
	void estimate_candidates(
		CodewordConstRange guesses,
		CodewordConstRange secrets,
		int lower_bound[]) const
	{
		const int secret_count = secrets.size();
		const int guess_count = guesses.size();
		const int table_size = (int)Feedback::maxValue(e.rules()) + 1;
		std::vector<unsigned int> frequency_cache(table_size*guess_count);

		// Partition the secrets using each candidate guess.
		//int max_b = 0; // maximum branching factor of non-perfect feedbacks
		Feedback perfect = Feedback::perfectValue(e.rules());
		for (int i = 0; i < guess_count; ++i)
		{
			Codeword guess = guesses[i];
			FeedbackFrequencyTable freq = e.frequency(e.compare(guess, secrets));

			std::copy(freq.begin(), freq.end(), frequency_cache.begin() + i*table_size);

#if 0
			int b = freq.nonzero_count();
			if (freq[perfect.value()] > 0)
				--b;
			if (b > max_b)
				max_b = b;
#endif
		}
		// assert(frequency_cache.size() == table_size*guess_count);

		// For each guess, compute the lower bound of total number of steps
		// required to reveal all secrets starting from that guess.
		for (int i = 0; i < guess_count; ++i)
		{
			int lb = secret_count;
			int j0 = i * table_size;
			for (int j = 0; j < table_size; ++j)
			{
				if (j != perfect.value())
				{
					// Note: we must NOT use max_b because canonical guesses
					// can change after we make a guess. So max_b no longer
					// works.
					//lb += simple_estimate(frequency_cache[j0+j], max_b);
					lb += simple_estimate(frequency_cache[j0+j]);
				}
			}
			lower_bound[i] = lb;
		}
	}
#endif
};

} // namespace Mastermind::Heuristics

/// Real-time optimal strategy. To be practical, the search space
/// must be small. For example, it works with Mastermind rules (p4c10r),
/// but probably not larger.
/// @ingroup Optimal
class OptimalStrategy
{
public:

	/// Returns the name of the strategy.
	virtual std::string name() const { return "optimal"; }

	/// Makes a guess.
	virtual Codeword make_guess(
		CodewordConstRange possibilities,
		CodewordConstRange candidates) const;
};

} // namespace Mastermind

#endif // MASTERMIND_OPTIMAL_STRATEGY_HPP
