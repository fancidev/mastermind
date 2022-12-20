#ifndef MASTERMIND_CODE_BREAKER_HPP
#define MASTERMIND_CODE_BREAKER_HPP

#include <string>
#include <memory>
#include "Engine.hpp"
#include "Strategy.hpp"
#include "ObviousStrategy.hpp"
#include "StrategyTree.hpp"
// #include "State.hpp"
#include "Equivalence.hpp"

namespace Mastermind
{

/// Stores options to control the behavior of a code breaker routine.
struct CodeBreakerOptions
{
	/// Indicates whether to make an obvious guess first.
	bool optimize_obvious; 

	/// Indicates whether to make a guess only from the remaining
	/// possibilities.
	bool possibility_only; 
	
	/// Creates a default set of "best" options.
	CodeBreakerOptions()
		: optimize_obvious(true), possibility_only(false) { }
};

/// Stores statistics about a code breaker routine.
struct CodeBreakerStatistics
{
};

// Free-standing function that makes a guess.
Codeword MakeGuess(
	Engine &e,
	CodewordConstRange secrets,
	// State &state,
	Strategy *strat,
	const EquivalenceFilter *filter,
	const CodeBreakerOptions &options);

// Free-standing function that builds a strategy tree.
StrategyTree BuildStrategyTree(
	const Engine *e, 
	Strategy *strat, 
	const EquivalenceFilter *filter,
	const CodeBreakerOptions &options);

/// Helper class that uses a given strategy to break a code.
class CodeBreaker
{
	/// Algorithm engine.
	Engine &e;

	/// The strategy used to make guesses.
	std::unique_ptr<Strategy> _strategy;

	/// The equivalence filter used to find canonical guesses.
	std::unique_ptr<	EquivalenceFilter> _filter;

	/// Options.
	CodeBreakerOptions _options;

	/// Set of possibilities. This set is updated on the way.
	CodewordList _possibilities;

public:

	/// Creates a code breaker using the given engine and strategy.
	CodeBreaker(
		Engine &engine, 
		std::unique_ptr<Strategy> strategy, 
		std::unique_ptr<EquivalenceFilter> filter,
		const CodeBreakerOptions &options)
		: e(engine), 
		_strategy(std::move(strategy)),
		_filter(std::move(filter)),
		_options(options),
		_possibilities(e.universe().begin(), e.universe().end())
	{ 
	}

	/// Returns the strategy used.
	const Strategy* strategy() const { return _strategy.get(); }

	/// Adds a constraint (guess:feedback pair). Updates the state of
	/// the game.
	///
	/// Note that the guess provided here may not be the guess that the
	/// code breaker suggests previously (via <code>MakeGuess()</code>).
	/// It might even happen that the code breaker has not been asked
	/// to make a guess at all. If a particular implementation does not
	/// support such behavior (e.g. one that makes guesses according to
	/// a pre-built strategy tree), it must throw an exception.
	void AddConstraint(const Codeword &guess, Feedback feedback)
	{
		_possibilities = e.filterByFeedback(_possibilities, guess, feedback);
		_filter->add_constraint(guess, feedback, _possibilities);
	}

	/// Makes a guess.
	Codeword MakeGuess()
	{
		return Mastermind::MakeGuess(
			e, _possibilities, _strategy.get(), _filter.get(), _options);
	}
};

}

#endif // MASTERMIND_CODE_BREAKER_HPP
