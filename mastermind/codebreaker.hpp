#pragma once

#include "codeword.hpp"
#include <memory>
#include <string_view>

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

    /// Advances the state of the game with the given constraint.
    ///
    /// The guess is not necessarily the one returned by a prior call
    /// to `make_guess`, and `make_guess` is not required to have been
    /// called.  If the implementation does not support arbitrary guess,
    /// it should throw an exception.
    ///
    /// The behavior is undefined if the (guess, response) pair is
    /// inconsistent, i.e. if it leads to no admissible secrets.
    virtual void step(const Constraint &constraint) = 0;

protected:
    CodeBreaker() = default;
};

/// Creates a code breaker that picks an arbitrary potential secret as
/// the next guess.  (In the current implementation, 'arbitrary' means
/// the lexicographically least.)
std::unique_ptr<CodeBreaker> create_simple_code_breaker(
    const CodewordRules &rules);

/// Creates a code breaker that picks a guess according to a heuristic score.
///
/// For each candidate guess, the set of potential secrets is partitioned
/// according to their feedback value against this guess.
///
/// The size of the partitions (but not their contents) are then passed
/// to a *heuristic function*, which generates a *score*.  The candidate
/// with the lowest score is chosen as the next guess.  To break ties, a
/// candidate that is also admissible (i.e. one that could be the secret)
/// is preferred over one that is not.  If there are still ties, an
/// arbitrary one is chosen.
///
/// Note: Ideally, the heuristic score should have adequately captured
/// the information in the partition, such that favoring an admissible
/// candidate is redundant.  However, in practice this does improve the
/// result (slightly) in some cases, e.g. "-r bc -s entropy".  It is
/// also compatible with prior art that had this treatment.
///
std::unique_ptr<CodeBreaker> create_heuristic_breaker(
    const CodewordRules &rules, std::string_view name);

} // namespace mastermind
