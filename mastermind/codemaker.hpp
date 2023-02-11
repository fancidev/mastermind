#pragma once

#include "codeword.hpp"
#include <memory>
#include <optional>

namespace mastermind {

/// Represents an agent that keeps the secret and gives feedbacks to guesses.
///
/// The codemaker is required to give *consistent* responses, i.e. for any
/// challenge sequence, it must respond with a feedback sequence such that
/// there exists at least one codeword that satisfies all the constraints.
class Codemaker
{
public:
    virtual ~Codemaker() {}

    /// Makes a response to the given guess.
    ///
    /// Note: The response must be consistent with all previous responses
    /// given by this code breaker.
    virtual Feedback respond(const Codeword &guess) = 0;

protected:
    Codemaker() = default;
};

/// Creates a codemaker with the given secret.
std::unique_ptr<Codemaker> create_static_codemaker(const Codeword &secret);

/// Creates a codemaker that dynamically adjusts the secret after each guess.
std::unique_ptr<Codemaker> create_dynamic_codemaker(
    const CodewordRules &rules);

Codeword sample(const CodewordRules &rules);

} // namespace mastermind
