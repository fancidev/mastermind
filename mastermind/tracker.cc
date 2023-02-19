#include "tracker.hpp"

#include <algorithm>
#include <compare>
#include <iterator>

namespace mastermind {

/// Returns a vector v with (m + 1) elements (m := codeword length), where
/// v[j] := the size of the sub-population where the first j letters in the
/// codeword are fixed.
template <size_t M>
static constexpr std::array<size_t, M + 1>
get_sub_population_sizes(const CodewordRules &rules) noexcept
{
    std::array<size_t, M + 1> sizes {};
    size_t count = 1;
    const size_t n = rules.alphabet_size();
    const size_t m = rules.codeword_size();
    for (Position j = m; j > 0; --j)
    {
        sizes[j] = count;
        if (rules.is_heterogram())
            count *= (n - j + 1);
        else
            count *= n;
    }
    sizes[0] = count;
    return sizes;
}

template <size_t M>
static constexpr Codeword
get_codeword_at(size_t index,
                const CodewordRules &rules,
                const std::array<size_t, M + 1> &sub_population_sizes) noexcept
{
    const size_t n = rules.alphabet_size();
    const size_t m = rules.codeword_size();

    std::array<Letter, MAX_CODEWORD_SIZE> letters {};
    if (rules.is_heterogram())
    {
        std::array<Letter, MAX_ALPHABET_SIZE> alphabet {};
        std::iota(alphabet.begin(), alphabet.end(), Letter(0));
        for (Position j = 0; j < m; j++)
        {
            size_t i = index / sub_population_sizes[j + 1];
            letters[j] = alphabet[i];
            std::copy(alphabet.cbegin() + i + 1,
                      alphabet.cbegin() + n - j,
                      alphabet.begin() + i);
            index %= sub_population_sizes[j + 1];
        }
    }
    else
    {
        for (Position j = 0; j < m; j++)
        {
            size_t i = index / sub_population_sizes[j + 1];
            letters[j] = static_cast<Letter>(i);
            index %= sub_population_sizes[j + 1];
        }
    }
    return Codeword(letters.data(), letters.data() + m);
}

CodewordSet::CodewordSet(const CodewordRules &rules) : _rules(rules)
{
    std::array<size_t, MAX_CODEWORD_SIZE + 1> sub_population_sizes(
        get_sub_population_sizes<MAX_CODEWORD_SIZE>(rules));

    size_t n = sub_population_sizes[0];
    _list.reserve(n);
    for (size_t index = 0; index < n; ++index)
        _list.push_back(get_codeword_at<MAX_CODEWORD_SIZE>(index, rules, sub_population_sizes));

    const size_t m = rules.codeword_size();
    CodewordMorphism::PositionMap position_map(m);
    CodewordMorphism::LetterMap letter_map;
    do
    {
        _morphisms.push_back({position_map, letter_map});
    }
    while (std::next_permutation(position_map.begin(), position_map.begin() + m));
}

void CodewordSet::push_constraint(const Constraint &constraint)
{
    _constraints.push_back(constraint);

    auto it = std::copy_if(_list.begin(),
                           _list.end(),
                           _list.begin(),
                           constraint);
    _list.resize(it - _list.begin());

    const size_t num_used = _used.count();

    // Update the canonical mappings.
    const LetterSequence letters(constraint.guess);
    LetterSequence canonical;
    const size_t m = _rules.codeword_size();
    size_t out = 0;
    for (size_t in = 0; in < _morphisms.size(); in++)
    {
        // Make copy of morphism because we may extend its letter map.
        CodewordMorphism morph(_morphisms[in]);

        // Permute positions
        LetterSequence permuted(letters);
        for (size_t j = 0; j < m; j++)
            permuted[morph.position_map[j]] = letters[j];

        // Permute letters
        Letter next = static_cast<Letter>(num_used);
        for (size_t j = 0; j < m; j++)
        {
            Letter i = permuted[j];
            Letter ii = morph.letter_map[i];
            if (ii == CodewordMorphism::LetterMap::not_mapped)
                ii = morph.letter_map[i] = next++;
            permuted[j] = ii;
        }

        // Keep this morphism if the permuted sequence is canonical.
        if (in == 0)
        {
            canonical = permuted;
            _morphisms[out++] = morph;
        }
        else
        {
            auto cmp = (permuted <=> canonical);
            if (cmp < 0)
            {
                out = 0;
                _morphisms[out++] = morph;
                canonical = permuted;
            }
            else if (cmp == 0)
            {
                _morphisms[out++] = morph;
            }
        }
    }
    assert(out > 0);
    _morphisms.resize(out);

    for (size_t j = 0; j < m; j++)
        _used.set(letters[j]);
}

std::vector<Codeword> CodewordSet::get_canonical_guesses() const
{
    // Build a list of all automorphisms of the canonical sequence.
    // It is equal to each morphism composed with the inverse of an
    // arbitrary fixed morphism.
    assert(_morphisms.size() > 0);
    CodewordMorphism fixed_inverse { _morphisms[0].inverse() };
    std::vector<CodewordMorphism> automorphisms;
    automorphisms.reserve(_morphisms.size());
    for (const CodewordMorphism &morph : _morphisms)
    {
        automorphisms.push_back(morph.composed_with(fixed_inverse));
    }

    // "Complete" the fixed inverse morphism to map canonical guesses
    // back into original space later.
    fixed_inverse.letter_map = fixed_inverse.letter_map.complete();

    // Check each codeword in the population.  A codeword is canonical
    // if it is mapped to itself or a lexicographically larger codeword
    // under every automorphisms.
    std::vector<Codeword> canonical_guesses;
    CodewordSet population(_rules);
    const CodewordSize m = _rules.codeword_size();
    for (const Codeword &guess : population)
    {
        std::array<Letter, MAX_CODEWORD_SIZE> letters;
        std::copy(guess.begin(), guess.end(), letters.begin());

        std::strong_ordering cmp = std::strong_ordering::equal;
        for (CodewordMorphism morph : automorphisms)
        {
            Letter next = static_cast<Letter>(_used.count());

            // To compare h := morphism(g) against g in lexical order,
            // we compare h[j] vs g[j] for each j = 0, ..., m-1.
            // Note that h[j] = letter_map[g[position_inv[j]].
            for (Position j = 0; j < m; j++)
            {
                Letter i = letters[j];
                Letter ii = letters[morph.position_map.inverse()[j]];
                if (ii >= next)
                    ii = morph.letter_map[ii] = next++;
                else
                    ii = morph.letter_map[ii];
                cmp = (ii <=> i);
                if (cmp != std::strong_ordering::equal)
                    break;
            }

            if (cmp < 0) // not canonical
                break;
        }

        if (cmp >= 0) // guess is canonical
        {
            // Convert canonical guess back to the original space.
            std::array<Letter, MAX_CODEWORD_SIZE> original;
            for (Position j = 0; j < m; j++)
            {
                Letter ii = letters[fixed_inverse.position_map.inverse()[j]];
                original[j] = fixed_inverse.letter_map[ii];

            }
            canonical_guesses.push_back(
                Codeword(original.begin(), original.begin() + m));
        }
    }
    return canonical_guesses;
}

} // namespace mastermind
