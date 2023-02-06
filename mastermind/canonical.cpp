#include "canonical.hpp"
#include <iostream>

namespace mastermind {

std::vector<CanonicalCodewordSequence>
get_canonical_guesses(const CanonicalCodewordSequence &sequence,
                      std::span<const Codeword> candidate_guesses)
{
    std::vector<CanonicalCodewordSequence> result;

    CanonicalCodewordSequence current(sequence);
    for (Codeword guess : candidate_guesses)
    {
        if (!current.contains(guess))
        {
            if (current.extend(guess))
            {
                result.push_back(current);
                // std::cout << "Canonical: " << guess << std::endl;
                current = sequence;
            }
        }
    }
    return result;
}

} // namespace mastermind
