#ifndef MASTERMIND_STRATEGY_HPP
#define MASTERMIND_STRATEGY_HPP

#include <iostream>
#include <string>
#include "Engine.hpp"

namespace Mastermind {

/// <summary>
/// Defines the objectives of a strategy. Here we implement three types of
/// objectives, in order of their strength (i.e. an objective with a larger
/// numeric value if strictly more optimal than an objective with a smaller
/// numeric value).
/// </summary>
enum StrategyObjective
{
	/// <summary>
    /// Minimize the total number of guesses needed to reveal all secrets.
    /// </summary>
	MinSteps = 1,

	/// In addition to @c MinSteps, also minimize the maximum number of
	/// guesses required to reveal any given secret.
 	MinDepth = 2,

	/// In addition to @c MinDepth, also minimize the number of secrets
	/// revealed by the most number of guesses.
	MinWorst = 3
};

/// <summary>
/// Defines a set of constraints that must be satisfied by a strategy.
/// </summary>
struct StrategyConstraints
{
	/// <summary>
    /// Gets or sets the maximum number of guesses allowed to reveal a secret.
    /// </summary>
	unsigned char max_depth;

	/// Flag indicating whether to make a guess only from the remaining
	/// possibilities.
	bool pos_only;

	/// Flag indicating if an obvious guess can be used if available.
	bool use_obvious;

	/// Flag indicating whether to find the last one among all optimal strategies.
	bool find_last;

	/// <summary>
    /// Creates a set of default constraints, which puts no restrictions
    /// on a strategy.
    /// </summary>
	StrategyConstraints()
		: max_depth(100), pos_only(false), use_obvious(true), find_last(false)
	{ }

	/// Creates a set of constraints using the provided parameters.
	//StrategyConstraints(unsigned char _max_depth, bool _pos_only, bool _use_obvious)
	//	: max_depth(_max_depth), pos_only(_pos_only), use_obvious(_use_obvious) { }
};

#ifndef USE_UNNAMED_STRUCT 
#define USE_UNNAMED_STRUCT 0
#endif

/**
 * Represents the cost of a strategy in terms of the number of guesses
 * required to reveal the secrets.
 *
 * The cost of a strategy consists of the following parts, in order of
 * decreasing priority:
 * - @c steps: the total number of guesses needed to reveal all
 *   secrets, excluding the initial guess.
 * - @c depth: the maximum number of guesses needed to reveal a
 *   secret, excluding the initial guess.
 * - @c worst: the number of secrets revealed by the worst number
 *   of steps.
 *
 * @remarks The machine must be little-endian for the comparison routine
 *      to work correctly.
 * @ingroup Optimal
 */
struct StrategyCost
{
#if USE_UNNAMED_STRUCT
	union
	{
		unsigned long long value;
		struct
		{
#endif
			unsigned short worst; // number of secrets revealed using max depth
			unsigned short depth; // number of guesses needed in the worst case
			unsigned int   steps; // total number of steps to reveal all secrets
#if USE_UNNAMED_STRUCT
		};
	};
#endif

	StrategyCost() : worst(0), depth(0), steps(0) { }
	StrategyCost(unsigned int _steps, unsigned short _depth, unsigned short _worst)
		: worst(_worst),  depth(_depth), steps(_steps) { }

	StrategyCost& operator += (const StrategyCost &c)
	{
		steps += c.steps;
		depth = std::max(depth, c.depth);
		return *this;
	}

	StrategyCost& operator -= (const StrategyCost &c)
	{
		assert(steps >= c.steps);
		// @todo The following assertion often fails; we need to amend the
		// optimal strategy code to avoid the use of -=.
		//assert(depth >= c.depth);
		steps -= c.steps;
		return *this;
	}

	/// Tests whether the cost is zero.
	bool operator ! () const { return steps == 0; }

	/// Tests whether the cost is not zero.
	//operator void* () const { return value == 0 ? 0 : (void*)this; }
};

/// Checks if strategy cost @c a is strictly superior to (i.e. lower than)
/// strategy cost @c b with regard to the objective @c obj.
inline bool superior(const StrategyCost &a, const StrategyCost &b, StrategyObjective obj)
{
	if (a.steps < b.steps)
		return true;
	if (b.steps < a.steps)
		return false;
	if (obj <= MinSteps)
		return false;

	if (a.depth < b.depth)
		return true;
	if (b.depth < a.depth)
		return false;
	if (obj <= MinDepth)
		return false;

	return a.worst < b.worst;
}

/// Compares the costs of two strategies.
/// @ingroup strat
inline bool operator < (const StrategyCost &c1, const StrategyCost &c2)
{
#if USE_UNNAMED_STRUCT
	return c1.value < c2.value;
#else
	return superior(c1, c2, MinWorst);
	//return *(const long long *)&c1 < *(const long long *)&c2;
#endif
}

/// Tests whether two strategy costs are equal.
/// @ingroup strat
inline bool operator == (const StrategyCost &c1, const StrategyCost &c2)
{
#if USE_UNNAMED_STRUCT
	return c1.value == c2.value;
#else
	return c1.steps == c2.steps && c1.depth == c2.depth && c1.worst == c2.worst;
	//return memcmp(&c1, &c2, sizeof(StrategyCost)) == 0;
#endif
}

/// Outputs strategy cost to a stream.
/// @ingroup strat
inline std::ostream& operator << (std::ostream &os, const StrategyCost &c)
{
	return os << c.steps << ':' << c.depth;
}

/// Substracts a component from an accumulated cost.
/// @ingroup strat
inline StrategyCost operator - (const StrategyCost &c1, const StrategyCost &c2)
{
	return StrategyCost(c1) -= c2;
}

/// Function object that compares the costs of two strategies according to an
/// objective.
class StrategyCostComparer
{
	StrategyObjective obj;

public:

	/// Constructs a strategy cost comparer with the given objective.
	StrategyCostComparer(StrategyObjective objective) : obj(objective) { }

	/// Returns the objective of the strategy cost comparer.
	StrategyObjective objective() const { return obj; }

	/// Compares the costs of two strategies.
	bool operator () (const StrategyCost &a, const StrategyCost &b) const
	{
		if (a.steps < b.steps)
			return true;
		if (b.steps < a.steps)
			return false;
		if (obj <= MinSteps)
			return false;

		if (a.depth < b.depth)
			return true;
		if (b.depth < a.depth)
			return false;
		if (obj <= MinDepth)
			return false;

		return a.worst < b.worst;
	}
};

#if 0
/// Checks if strategy cost @c a is strictly superior to (i.e. lower than)
/// strategy cost @c b with regard to the objective @c obj.
inline bool superior(const StrategyCost &a, const StrategyCost &b, StrategyObjective obj)
{
	if (a.steps < b.steps)
		return true;
	if (b.steps < a.steps)
		return false;
	if (obj <= MinSteps)
		return false;

	if (a.depth < b.depth)
		return true;
	if (b.depth < a.depth)
		return false;
	if (obj <= MinDepth)
		return false;

	return a.worst < b.worst;
}
#endif

/**
 * Interface for a Mastermind strategy.
 *
 * On the micro-scope, a strategy is a function that receives as input
 *   - a list of possibilities
 *   - a list of candidate guesses
 * and returns as output a guess to make.
 *
 * A few strategies implemented (or to be implemented) include:
 *   - simple strategy
 *   - heuristic strategy
 *   - optimal strategy
 *   - playback strategy (reads from a strategy tree)
 *
 * @ingroup strat
 */
struct Strategy
{
	/// Returns the name of the strategy.
	virtual std::string name() const = 0;

	/**
	 * Makes a guess.
	 *
	 * @param possibilities List of remaining possibilities.
	 * @param candidates    List of (suggestive) candidate guesses.
	 *      Note that the function is not required to return a guess
	 *      from this list.
	 * @returns If successful, returns the guess to make. The guess
	 *      must be from @c possibility OR @c candidates. If failed,
	 *      returns <code>Codeword::emtpyValue()</code>. The condition
	 *      under which the call fails is implementation-specific.
	 */
	virtual Codeword make_guess(
		CodewordConstRange possibilities,
		CodewordConstRange candidates) const = 0;
};

} // namespace Mastermind

#endif // MASTERMIND_STRATEGY_HPP
