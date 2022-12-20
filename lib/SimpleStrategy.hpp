#ifndef MASTERMIND_SIMPLE_STRATEGY_HPP
#define MASTERMIND_SIMPLE_STRATEGY_HPP

#include "Strategy.hpp"

namespace Mastermind {

/// Simple strategy that always guesses the first codeword from the
/// set of remaining possibilities.
/// @ingroup Simple
class SimpleStrategy : public Strategy
{
public:

	/// Constructs the strategy.
	SimpleStrategy() { }

	/// Returns the name of the strategy.
	virtual std::string name() const
	{
		return "simple";
	}

	/// Returns the first codeword from the possibility set as the
	/// guess. If the possibility set is empty, returns an empty 
	/// codeword.
	virtual Codeword make_guess(
		CodewordConstRange possibilities, 
		CodewordConstRange /* candidates */) const
	{
		if (possibilities.empty())
			return Codeword();
		else
			return *possibilities.begin();
	}
};

} // namespace Mastermind

#endif // MASTERMIND_SIMPLE_STRATEGY_HPP
