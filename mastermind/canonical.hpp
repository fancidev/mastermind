#pragma once

#include "codeword.hpp"

#include <bit>
#include <bitset>
#include <cassert>
#include <concepts>

#include "thirdparty/fixed_capacity_vector"

#include <iostream>
#include <cstdint>
#include <numeric>
#include <vector>

namespace mastermind {

/// Represents a permutation of the index set {0, ..., n-1} of type T,
/// supporting up to Capacity indices.
template <std::integral T, size_t Capacity>
class Permutation
{
public:
    /// Type of the index.
    using Index = T;

    /// Creates an identity permutation of the given size.
    constexpr explicit Permutation(size_t size) noexcept : _map(size)
    {
        std::iota(_map.begin(), _map.end(), Index(0));
    }

    /// Updates the permutation to the next one in lexicographical order.
    ///
    /// Returns `false` if the updated permutation is the identity
    /// permutation, or `true` otherwise.
    constexpr bool next() noexcept
    {
        return std::next_permutation(_map.begin(), _map.end());
    }

    /// Returns the number of elements in the index set being permuted.
    constexpr size_t size() const noexcept { return _map.size(); }

    /// Returns the image of `index` under this permutation.
    constexpr Index map(const Index &index) const noexcept
    {
        assert(index >= 0 && index < _map.size());
        return _map[index];
    }

    /// Returns `true` if this is the identity permutation.
    constexpr bool is_identity() const noexcept
    {
        const size_t m = _map.size();
        for (Index j = 0; j < m; j++)
        {
            if (_map[j] != j)
                return false;
        }
        return true;
    }

private:
    /// Stores the inverse permutation, because it is what we will use.
    std::experimental::fixed_capacity_vector<Index, Capacity> _map;
};

/// Represents a bijection from a k-subset of the index set {0, ..., n-1}
/// to the index set {0, ..., k-1}.  In addition, the bijection is updated
/// automatically when an unmapped index is queried -- it is automatically
/// mapped to the least image.
template <std::integral Index, size_t Capacity>
class PartialPermutation
{
public:
    static constexpr Index not_mapped = Index(-1);

    /// Creates an empty bijection.
    constexpr explicit PartialPermutation(size_t n) noexcept
      : _map(n, not_mapped), _next_image(0)
    {
        assert(n >= 0 && n <= Capacity);
    }

    /// Returns `true` if every index in the range {0, ..., n-1} is
    /// explicitly mapped.
    constexpr bool is_full() const noexcept
    {
        return _next_image == _map.size();
    }

    /// Returns `true` if every index in the range {0, ..., n-1} is
    /// explicitly mapped to itself.
    constexpr bool is_identity() const noexcept
    {
        Index index = Index(0);
        for (Index image : _map)
        {
            if (image != index++)
                return false;
        }
        return true;
    }

    /// Returns the image of `index` under this bijection.
    ///
    /// If `index` is not explicitly mapped, it is automatically mapped
    /// to the least unmapped image and the bijection object is updated.
    /// That is why this member function is not `const`.
    constexpr Index map(Index index) noexcept
    {
        assert(index >= 0 && index < _map.size());
        if (_map[index] == not_mapped)
            _map[index] = _next_image++;
        return _map[index];
    }

private:
    /// `_map[index]` is the image under the bijection, or `not_mapped`
    /// if `index` is not explicitly mapped.
    std::experimental::fixed_capacity_vector<Index, Capacity> _map;

    /// The next image to map an unmapped index to.
    Index _next_image;
};

/// Bundles an inverse position permutation with a partial alphabet
/// permutation.
struct CodewordPermutation
{
    /// Represents a permutation of (position) indices {0, ..., m-1}.
    using PositionPermutation = Permutation<PositionIndex, MAX_CODEWORD_LENGTH>;

    /// Represents a bijection from a subset of k alphabet indices
    /// in the range {0, ..., n-1} to the first k alphabet indices
    /// {0, ..., k-1}.  The bijection grows automatically when an
    /// unmapped index is mapped.
    using AlphabetPermutation =
        PartialPermutation<AlphabetIndex, MAX_ALPHABET_SIZE>;

    PositionPermutation inverse_position_permutation;
    AlphabetPermutation partial_alphabet_permutation;

    explicit constexpr CodewordPermutation(const CodewordRules &rules) noexcept
      : inverse_position_permutation(rules.codeword_length()),
        partial_alphabet_permutation(rules.alphabet_size())
    {
    }
};

/// Represents the set of all automorphisms (permutation of letters and
/// positions) that map a guess sequence to itself.
class AutomorphismGroup
{
public:
    /// Creates the automorphism group for an empty guess sequence.
    explicit AutomorphismGroup(const CodewordRules &rules)
    {
        CodewordPermutation perm(rules);
        do
        {
            _perms.push_back(perm);
        }
        while (perm.inverse_position_permutation.next());
    }

    /// Returns `true` if this automorphism group contains exactly one
    /// automorphism, which is necessarily the identity morphism.
    constexpr bool is_singleton() const noexcept
    {
        assert(!_perms.empty());
        if (_perms.size() == 1)
        {
            const CodewordPermutation &perm = _perms.front();
            assert(perm.inverse_position_permutation.is_identity());
            if (perm.partial_alphabet_permutation.is_full())
            {
                assert(perm.partial_alphabet_permutation.is_identity());
                return true;
            }
        }
        return false;
    }

    /// If `guess` is canonical under the current automorphism, update
    /// the current automorphism group by appending `guess` to the
    /// sequence, and returns `true`.  Otherwise, returns `false`.
    ///
    /// `guess` is canonical if none of the permutations in the current
    /// automorphism group maps `guess` to a lexicographically smaller
    /// image.
    ///
    /// Note: The updated automorphism group may be a singleton.
    ///
    /// TODO: Some of the code branches might be optimized by
    /// TODO: studying the structure of automorphism.
    bool refine(const Codeword &guess)
    {
        // TODO: support Codeword::begin() and Codeword::end()
        const LetterSequence letters = guess.letters();

        // The guess is canonical if and only if for every automorphism
        // in the group, the mapped `guess` is greater than or equal to
        // `guess` in lexical order.

        std::vector<CodewordPermutation> perms;

        for (CodewordPermutation perm : _perms)
        {
            // To check whether h := perm(g) < g in lexical order,
            // we compare h[j] vs g[j] for each j = 0, ..., m-1.
            // Note that h[j] = letter_map[g[inverse_position_map[j]].

            bool is_automorphism = true;
            const PositionSize m = perm.inverse_position_permutation.size();
            for (PositionIndex j = 0; j < m && is_automorphism; j++)
            {
                AlphabetIndex index = letters[j];
                PositionIndex pos = perm.inverse_position_permutation.map(j);
                AlphabetIndex image = perm.partial_alphabet_permutation.map(letters[pos]);

                if (image < index) // guess is not canonical
                    return false;
                if (image > index) // not an automorphism
                    is_automorphism = false;
            }

            if (is_automorphism)
                perms.push_back(perm);
        }

        _guess_sequence.push_back(guess); // may throw
        std::swap(perms, _perms);
        return true;
    }

    constexpr const std::vector<Codeword> &guess_sequence() const
    {
        return _guess_sequence;
    }

private:
    /// The guess sequence whose automorphism group is represented.
    std::vector<Codeword> _guess_sequence;

    /// List of all permutations that map `_guess_sequence` to itself.
    std::vector<CodewordPermutation> _perms;
};

///// Outputs a codeword permutation to a stream.
//inline std::ostream& operator << (std::ostream& os, const CodewordPermutation &p)
//{
//	// Output peg permutation.
//	os << "(";
//	for (int i = 0; i < MM_MAX_PEGS; ++i)
//	{
//		if (i > 0)
//			os << ' ';
//		os << (size_t)p.peg[i];
//	}
//	os << ") o (";
//
//	// Output color permutation.
//	for (int i = 0; i < MM_MAX_COLORS; ++i)
//	{
//		if (i > 0)
//			os << ' ';
//		os << (size_t)p.color[i];
//	}
//	os << ")";
//	return os;
//}
//
//inline void get_canonical_guesses(const CodewordRules &rules)
//{
//    CodewordPopulation population(rules);
//    std::vector<Codeword> candidate_guesses(population.begin(), population.end());
//
//    AutomorphismGroup group(rules);
//    AutomorphismGroup current(group);
//#if 1
//    for (Codeword guess : candidate_guesses)
//    {
//        if (current.refine(guess))
//        {
//            std::cout << "Canonical: " << guess << std::endl;
//            current = group;
//        }
//    }
//#else
//    Codeword guess(rules);
//    std::istringstream("0111") >> guess;
//    if (current.refine(guess))
//    {
//        std::cout << "Canonical: " << guess << std::endl;
//        current = group;
//    }
//    else
//    {
//        std::cout << "Not canonical" << std::endl;
//    }
//#endif
//}

std::vector<AutomorphismGroup>
get_canonical_guesses(const AutomorphismGroup &group,
                      std::span<const Codeword> candidate_guesses);

} // namespace mastermind

