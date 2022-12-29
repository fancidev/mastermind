#pragma once

#include "codeword.hpp"
#include <memory>

namespace mastermind {

/// Represents an agent that makes successive guesses to reveal a secret.
class CodeBreaker
{
public:
    virtual ~CodeBreaker() = default;

    /// Returns the name of the code breaker.  The returned pointer is
    /// valid during the lifetime of the code breaker object.
    virtual const char *name() const = 0;

    /// Computes and returns the guess to make in the current situation.
    ///
    /// Throws `std::runtime_error` if there are no admissible secrets.
    virtual Codeword make_guess() = 0;

    /// Advances the state of the game with the response for a given guess.
    ///
    /// The guess is not necessarily the one returned by a prior call
    /// to `make_guess`, and `make_guess` is not required to have been
    /// called.  If the implementation does not support arbitrary guess,
    /// it should throw an exception.
    ///
    /// The behavior is undefined if the (guess, response) pair is
    /// inconsistent, i.e. if it leads to no admissible secrets.
    virtual void step(const Codeword &guess, Feedback response) = 0;

protected:
    CodeBreaker() = default;
};

/// Creates a code breaker that picks an arbitrary potential secret as
/// the next guess.  (In the current implementation, 'arbitrary' means
/// the lexicographically least.)
std::unique_ptr<CodeBreaker> create_simple_code_breaker(
    const CodewordRules &rules);

} // namespace mastermind
