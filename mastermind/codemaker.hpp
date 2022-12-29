#pragma once

#include "codeword.hpp"
#include <memory>
#include <optional>

namespace mastermind {

/// Represents an agent that keeps the secret and responds to guesses with
/// feedbacks.
///
/// The code maker is required to give *consistent* responses, i.e. if
/// challenged with every codeword in the population, it must give exactly
/// one perfect response, and the rest responses must be consistent with
/// this secret.
class CodeMaker
{
public:
    virtual ~CodeMaker() {}

    /// Makes a response to the given guess.
    ///
    /// Note: The response must be consistent with all previous responses
    /// given by this code breaker.
    virtual Feedback respond(const Codeword &guess) = 0;

    /// Returns the secret, if one is available.
    ///
    /// Note: The standard code maker fixes the secret at the beginning,
    /// so the return value is always available.  A 'tricky' code maker
    /// could instead choose to fix the secret later and return a missing
    /// value.
    virtual std::optional<Codeword> secret() const = 0;

protected:
    CodeMaker() = default;
};

/// Creates a standard code maker, optionally using the given secret.
/// A random secret is generated if secret is not provided.
std::unique_ptr<CodeMaker> create_code_maker(
    const CodewordRules &rules,
    std::optional<Codeword> secret = std::optional<Codeword>());

} // namespace mastermind
