#pragma once

#include "codeword.hpp"

#include <array>
#include <cassert>
//#include <concepts>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <span>
#include <vector>

namespace mastermind {

/// Represents a mapping defined over the index set {0, ..., N-1} of type T
/// to itself.
template <class Index, size_t N>
class Mapping
{
public:
    /// Sentinel value returned by `map()` to indicate that an index
    /// is not mapped.
    constexpr static const Index not_mapped = Index(-1);

    /// Creates an empty mapping.
    constexpr Mapping() noexcept : _map{} //, _size{}
    {
        std::fill(_map.begin(), _map.end(), not_mapped);
    }

    /// Creates a mapping where the first `(end-begin)` indices are
    /// mapped to the images specified by `begin...end`.
    template <class Iter>
    constexpr Mapping(Iter begin, Iter end) noexcept : Mapping()
    {
        std::copy(begin, end, _map.begin());
//        _size = std::distance(begin, end);
    }

//    /// Returns the number of mapped indices.
//    constexpr size_t size() const noexcept { return _size; }

    /// Returns the image of `index`, or `not_mapped` if `index` is not mapped.
    constexpr Index map(const Index &index) const noexcept
    {
        assert(index >= 0 && index < _map.size());
        return _map[index];
    }

    /// Returns the image of `index` if it is mapped.  Otherwise, maps
    /// `index` to `image` and returns `image`.
    constexpr Index map_or_update(const Index &index, const Index &image) noexcept
    {
        assert(index >= 0 && index < _map.size());
        assert(image >= 0 && image < _map.size());
        if (_map[index] == not_mapped)
        {
            _map[index] = image;
//            _size++;
            return image;
        }
        else
            return _map[index];
    }

    /// Returns `true` if this mapping is the identity permutation of
    /// the first `n` indices, {0, ..., n-1}.
    constexpr bool is_identity(size_t n) const noexcept
    {
        assert(n <= _map.size());
        for (std::size_t j = 0; j < n; j++)
        {
            if (_map[j] != j)
                return false;
        }
        for (std::size_t j = n; j < _map.size(); j++)
        {
            if (_map[j] != not_mapped)
                return false;
        }
        return true;
    }

private:
    /// `_map[i]` is the image of `i` or `not_mapped`.
    std::array<Index, N> _map;

//    /// Number of mapped indices.
//    std::size_t _size;
};

/// Represents a permutation of (position) indices {0, ..., m-1}.
using PositionPermutation = Mapping<PositionIndex, MAX_CODEWORD_LENGTH>;

/// Represents a permutation of the first k letters of the alphabet.
using AlphabetPermutation = Mapping<AlphabetIndex, MAX_ALPHABET_SIZE>;

/// Bundles an inverse position permutation and a partial alphabet
/// permutation.
struct CodewordMorphism
{
    PositionPermutation inverse_position_permutation;
    AlphabetPermutation partial_alphabet_permutation;
};

/// Represents a canonical codeword sequence, i.e. one that is the
/// lexicographical minimum among all codeword sequences isomophic
/// to it.
class CanonicalCodewordSequence
{
public:
    /// Creates an empty sequence bound to the given rules.
    explicit CanonicalCodewordSequence(const CodewordRules &rules)
      : _rules(rules), _used_letters(0)
    {
        std::array<PositionIndex, MAX_CODEWORD_LENGTH> map;
        const auto begin = map.begin();
        const auto end = map.begin() + rules.codeword_length();
        std::iota(begin, end, PositionIndex(0));
        do
        {
            PositionPermutation inverse_position_perm(begin, end);
            AlphabetPermutation partial_alphabet_perm;
            CodewordMorphism morphism {
                inverse_position_perm, partial_alphabet_perm
            };
            _morphisms.push_back(morphism);
        }
        while (std::next_permutation(begin, end));
    }

//    /// Returns `true` if there is no other codeword sequence that is
//    /// isomorphic to this one under the given rules.
//    constexpr bool is_singleton() const noexcept
//    {
//        assert(!_morphisms.empty());
//        if (_morphisms.size() > 1)
//            return false;
//
//        assert(!_sequence.empty());
//        return false;
////
////        if (_perms.size() == 1)
////        {
////            const CodewordPermutation &perm = _perms.front();
////            assert(perm.inverse_position_permutation.is_identity());
////            if (perm.partial_alphabet_permutation.is_full())
////            {
////                assert(perm.partial_alphabet_permutation.is_identity());
////                return true;
////            }
////        }
////        return false;
//    }

    /// Appends `guess` to the sequence and returns `true` if the resulting
    /// sequence is canonical.  Otherwise, returns `false` and leaves the
    /// object unchanged.
    bool extend(const Codeword &guess)
    {
        std::array<AlphabetIndex, MAX_CODEWORD_LENGTH> letters;
        guess.copy_letters(letters.begin());
        size_t used_letters = _used_letters;

        // `guess` is canonical if and only if under every automorphism
        // that maps the current sequence to itself, the image of `guess`
        // is greater than or equal to `guess` in lexicographical order.

        std::vector<CodewordMorphism> morphisms;
        for (CodewordMorphism morphism : _morphisms)
        {
            AlphabetIndex next = _used_letters;

            // To compare h := morphism(g) against g in lexical order,
            // we compare h[j] vs g[j] for each j = 0, ..., m-1.
            // Note that h[j] = letter_map[g[inverse_position_map[j]].
            const PositionSize m = _rules.codeword_length();
            bool is_automorphism = true;
            for (PositionIndex j = 0; j < m; j++)
            {
                const AlphabetIndex index = letters[j];
                const AlphabetIndex image =
                    morphism.partial_alphabet_permutation.map_or_update(
                        letters[morphism.inverse_position_permutation.map(j)],
                        next);
                if (image < index) // not canonical
                    return false;
                if (image > index)
                {
                    is_automorphism = false;
                    break;
                }
                if (image == next)
                    next++;
            }

            if (is_automorphism)
            {
                morphisms.push_back(morphism);
                used_letters = next;
            }
        }

        assert(!morphisms.empty());
        _sequence.push_back(guess);
        std::swap(morphisms, _morphisms);
        _used_letters = used_letters;
        return true;
    }

    /*constexpr*/ std::vector<Codeword>::const_iterator begin() const noexcept
    {
        return _sequence.begin();
    }

    /*constexpr*/ std::vector<Codeword>::const_iterator end() const noexcept
    {
        return _sequence.end();
    }

    /*constexpr*/ std::size_t size() const noexcept { return _sequence.size(); }

    /*constexpr*/ bool contains(const Codeword &g) const noexcept
    {
        return std::find(begin(), end(), g) != end();
    }

    /*constexpr*/ Codeword back() const noexcept
    {
        return _sequence.back();
    }

//    constexpr const std::vector<Codeword> &sequence() const
//    {
//        return _sequence;
//    }

private:
    /// Rules that codewords in this sequence conform to.
    CodewordRules _rules;

    /// The canonical codeword sequence.
    std::vector<Codeword> _sequence;

    /// List of all codeword morphisms that map `_sequence` to itself.
    std::vector<CodewordMorphism> _morphisms;

    /// Number of letters that appear in the sequence.
    std::size_t _used_letters;
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

std::vector<CanonicalCodewordSequence>
get_canonical_guesses(const CanonicalCodewordSequence &sequence,
                      std::span<const Codeword> candidate_guesses);

} // namespace mastermind

