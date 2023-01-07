#include "canonical.hpp"

namespace mastermind {

std::vector<AutomorphismGroup>
get_canonical_guesses(const AutomorphismGroup &group,
                      std::span<const Codeword> candidate_guesses)
{
    std::vector<AutomorphismGroup> result;

    AutomorphismGroup current(group);
    const std::vector<Codeword> &history = group.guess_sequence();
    for (Codeword guess : candidate_guesses)
    {
        if (std::find(history.begin(), history.end(), guess) != history.end())
            continue; // already guessed this!

        if (current.refine(guess))
        {
            result.emplace_back(std::move(current));
            // std::cout << "Canonical: " << guess << std::endl;
            current = group;
        }
    }
    return result;
}

} // namespace mastermind
