#ifndef MASTERMIND_HEURISTICS_HPP
#define MASTERMIND_HEURISTICS_HPP

#include <string>
#include <cmath>
#include <algorithm>
#include <functional>

#include "Engine.hpp"
#include "util/wrapped_float.hpp"

namespace Mastermind {
namespace Heuristics {

/**
 * Heuristic that scores a guess by the worst-case number of remaining 
 * possibilities (Knuth, 1976). If two guesses produce the same number of
 * of remaining possibilities in the worst case, the second-to-worst number
 * is compared, and so on.
 *
 * If a correction is applied (which is the default), a cell corresponding
 * to a perfect match is assumed to be empty, because it will not "remain"
 * after this guess.
 *
 * @ingroup Heuristic
 */ 
struct MinimizeWorstCase
{
	/// Flag indicating whether to make an adjustment to the score
	/// if the guess is among the remaining possibilities.
	bool apply_correction;

	/// Constructs the heuristic using the given policy.
	MinimizeWorstCase(bool _apply_correction = true)
		: apply_correction(_apply_correction) { }

	/// Type of the score (a sorted array of partition sizes in descending
	/// order). 
	typedef FeedbackFrequencyTable score_t;

	/// Short identifier of the heuristic function.
	std::string name() const { return "minmax"; }

	/// Returns a sorted array of partition sizes in descending order.
	score_t compute(const FeedbackFrequencyTable &freq) const
	{
		FeedbackFrequencyTable score(freq);
		if (apply_correction)
		{
			score[score.size()-1] = 0;
		}
		std::sort(score.begin(), score.end(), std::greater<unsigned int>());
		return score;
	}
};

/// Heuristic that scores a guess by the expected number of remaining
/// possibilities (Irving, 1979). The score to minimize is
/// <code>Sum{ n[i] * (n[i]/N) }</code>, or equivalently
/// <code>Sum{ n[i]^2 }</code>, where <code>n[i]</code> is the number
/// of elements in partition <code>i</code>.
/// @ingroup Heuristic
struct MinimizeAverage
{
	/// Type of the score (64-bit integer to avoid overflow if the
	/// partition size is larger than 2^16).
	typedef unsigned long long score_t;

	/// Flag indicating whether to make an adjustment to the score
	/// if the guess is among the remaining possibilities.
	bool apply_correction;

	/// Constructs the heuristic using the given policy.
	MinimizeAverage(bool _apply_correction = true)
		: apply_correction(_apply_correction) { }

	/// Name of the heuristic.
	std::string name() const
	{
		return apply_correction? "minavg" : "minavg~";
	}

	/// Computes the heuristic score - sum of squares of the size
	/// of each partition.
	score_t compute(const FeedbackFrequencyTable &freq) const
	{
		score_t s = 0;
		for (size_t i = 0; i < freq.size(); ++i)
		{
			score_t f = freq[i];
			s += f * f;
		}
		if (apply_correction && freq[freq.size()-1])
			--s;
		return s;
	}
};

/// A theoretically advanced heuristic that scores a guess as roughly
/// the expected number of further guesses needed (Neuwirth, 1982).
/// The precise definition is to maximize the entropy, defined as
/// <code>-Sum{ (n[i]/N) * log(n[i]/N) }</code>, which is equivalent
/// to minimizing <code>Sum{ (n[i]/N) * log(n[i]) }</code>. If we
/// interpret <code>log(n[i])</code> as an estimate of the number of
/// further guesses needed for a partition of size <code>n[i]</code>,
/// then we can interpret the the objective function as an
/// estimate of the expected number of further guesses needed.
/// As a side note, note that the base of the logrithm doesn't matter
/// in computing the score.
/// @ingroup Heuristic
struct MaximizeEntropy
{
#if 0
	/// Type of the score (double precision).
	typedef double score_t;
#else
	/// Type of the score (wrapped double to avoid numerical instability
	/// during comparison).
	typedef util::wrapped_float<double, 100> score_t;
#endif

	/// Flag indicating whether to make an adjustment to the score
	/// if the guess is among the remaining possibilities.
	bool apply_correction;

	/// Constructs the heuristic using the given policy.
	MaximizeEntropy(bool _apply_correction = true)
		: apply_correction(_apply_correction) { }

	/// Short identifier of the heuristic function.
	std::string name() const
	{
		return apply_correction? "entropy" : "entropy~";
	}

	/// Computes the heuristic score - negative of the entropy.
	score_t compute(const FeedbackFrequencyTable &freq) const
	{
		double s = 0.0;
		for (size_t i = 0; i < freq.size(); ++i)
		{
			unsigned int f = freq[i];
			if (f > 1)
			{
				s += std::log((double)f) * (double)f;
			}
		}
		if (apply_correction && freq[freq.size()-1]) // 4A0B
		{
			s -= 2.0 * std::log(2.0);
		}
		return score_t(s);
	}
};

/// An aggressive heuristic that scores a guess as the number of
/// partitions it produces. The rationale is that more partitions,
/// fewer steps.
/// @ingroup Heuristic
struct MaximizePartitions
{
	/// Type of the score (signed integer because we need to minimize
	/// the negative value of the number of partitions).
	typedef int score_t;

	/// Flag indicating whether to make an adjustment to the score
	/// if the guess is among the remaining possibilities.
	bool apply_correction;

	/// Constructs the heuristic using the given policy.
	MaximizePartitions(bool _apply_correction = true)
		: apply_correction(_apply_correction) { }

	/// Short identifier of the heuristic function.
	std::string name() const
	{
		return apply_correction? "parts" : "parts~";
	}

	/// Computes the heuristic score - negative of the number of
	/// partitions.
	score_t compute(const FeedbackFrequencyTable &freq) const
	{
		int score = 2 * (int)freq.nonzero_count();
		if (apply_correction && freq[freq.size()-1])
		{
			++score;
		}
		return -score;
	}
};

} // namespace Heuristics
} // namespace Mastermind

#endif // MASTERMIND_HEURISTICS_HPP
