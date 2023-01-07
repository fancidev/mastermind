#pragma once

#include "codeword.hpp"

#include <bit>
#include <bitset>
#include <cassert>

#include "thirdparty/fixed_capacity_vector"

#include <iostream>
#include <cstdint>
#include <numeric>
#include <vector>

namespace mastermind {

/// Represents a permutation of (position) indices {0, ..., m-1}.
class PositionPermutation
{
public:
    /// Creates an identity position permutation.
    constexpr explicit PositionPermutation(const CodewordRules &rules) noexcept
      : _inverse(rules.codeword_length())
    {
        std::iota(_inverse.begin(), _inverse.end(), PositionIndex(0));
    }

    /// Changes to the next position permutation in lexical order.
    ///
    /// Returns `false` if the new permutation is the identity permutation;
    /// otherwise returns `true`.
    constexpr bool next() noexcept
    {
        return std::next_permutation(_inverse.begin(), _inverse.end());
    }

    constexpr PositionSize size() const noexcept
    {
        return _inverse.size();
    }

    /// Returns the inverse of j' under this permutation.  That is,
    /// if j == inv(j'), then j' == map(j).
    constexpr PositionIndex inv(PositionIndex j_prime) const noexcept
    {
        assert(j_prime >= 0 && j_prime < _inverse.size());
        return _inverse[j_prime];
    }

    /// Returns `true` if this is the identity permutation.
    constexpr bool is_identity() const noexcept
    {
        const PositionSize m = _inverse.size();
        for (PositionIndex j = 0; j < m; j++)
        {
            if (_inverse[j] != j)
                return false;
        }
        return true;
    }

private:
    /// Stores the inverse permutation, because it is what we will use.
    std::experimental::fixed_capacity_vector<PositionIndex, MAX_CODEWORD_LENGTH>
        _inverse;
};

/// Represents a bijection from a subset of k (alphabet) indices in the range
/// {0, ..., n-1} to the first k (alphabet) indices (0, ..., k-1).
class AlphabetPermutation // rename to AlphabetBijection, or just bijection
{
public:
    static constexpr AlphabetIndex NotMapped = AlphabetIndex(-1);

    /// Creates an empty bijection.
    constexpr explicit AlphabetPermutation(const CodewordRules &rules) noexcept
      : _map(rules.alphabet_size(), NotMapped), _num_mapped(0) { }

//    constexpr bool is_mapped(AlphabetIndex i) const noexcept
//    {
//        assert(i >= 0 && i < _mapping.size());
//        return (_mapped & (size_t(1) << i)) != size_t(0);
//    }

    /// Returns `true` if the entire alphabet index set {0, ..., n-1}
    /// is mapped.
    constexpr bool is_permutation() const noexcept
    {
        return _num_mapped == _map.size();
    }

    /// Returns `true` if the bijection is the identity permutation.
    constexpr bool is_identity() const noexcept
    {
        const AlphabetSize n = _map.size();
        for (AlphabetIndex i = 0; i < n; i++)
        {
            if (_map[i] != i)
                return false;
        }
        return true;
    }

    /// Returns the number of letters mapped by this bijection.
    constexpr AlphabetSize size() const noexcept
    {
        return _num_mapped;
    }

    constexpr AlphabetIndex map(AlphabetIndex i) const noexcept
    {
        assert(i >= 0 && i < _map.size());
        return _map[i];
    }

    /// Updates the bijection by mapping `index` to `image`.
    ///
    /// `index` must not be mapped, and `image` must be the smallest
    /// available image.
    constexpr void map(AlphabetIndex index, AlphabetIndex image) noexcept
    {
        assert(index >= 0 && index < _map.size());
        assert(_map[index] == NotMapped);
        assert(image == _num_mapped);

        _map[index] = image;
        ++_num_mapped;
    }

private:
    /// `_map[i]` is the image under the bijection, or the constant
    /// NotMapped if `i` is not mapped.
    std::experimental::fixed_capacity_vector<AlphabetIndex, MAX_ALPHABET_SIZE>
        _map;

    /// `_mapped[i]` is `true` if and only if `_mapping[i]` is valid.
    AlphabetSize _num_mapped;
};

///// Represents a set of permutations of the codeword, where each permutation
///// is the composition of:
///// - a fixed permutation of (all) the positions;
///// - a fixed permutation of a subset of the alphabet; and
///// - any one of the permutations of the other subset of the alphabet.
//class CodewordPermutation
//{
//	/// `j' = _position_map[j]` maps position j to j', `0 <= j, j' < m`.
//    PositionMap _position_map;
//
//    /// `j = _inverse_position_map[j']` maps position j' back to j.
//    PositionMap _inverse_position_map;
//
//    /// `i' = _letter_map[i]` maps letter i to i` for `0 <= i < n` where
//    /// `_letter_mapped[i]` is true.  If `_letter_mapped[i]` is false,
//    /// `_letter_map[i]` is unspecified and should not be used.
//    LetterMap _letter_map;
//
//    /// `_letter_mapped[i]` is `true` if the i-th letter is mapped.
//    std::bitset<MAX_ALPHABET_SIZE> _letter_mapped;
//
//public:
//
//	/// Creates an identity permutation.
//	constexpr CodewordPermutation(const CodewordRules &rules) noexcept
//      : _position_map(rules.codeword_length()),
//        _inverse_position_map(rules.codeword_length()),
//        _letter_map(rules.alphabet_size()),
//        _letter_mapped()
//	{
//        std::iota(_letter_map.begin(), _letter_map.end(), AlphabetIndex(0));
//        std::iota(_position_map.begin(), _position_map.end(), PositionIndex(0));
//        _inverse_position_map = _position_map;
//        _letter_mapped.set();
//	}
//
//    /// Returns the mapped value of position j.
//    constexpr PositionIndex map_position(PositionIndex j) const noexcept
//    {
//        assert(j >= 0 && j < _position_map.size());
//        return _position_map[j];
//    }
//
//    /// Returns the position that is mapped to j.
//    constexpr PositionIndex inverse_map_position(PositionIndex j) const noexcept
//    {
//        assert(j >= 0 && j < _inverse_position_map.size());
//        return _inverse_position_map[j];
//    }
//
//    AlphabetIndex map_letter(AlphabetIndex i)
//    {
//        if (_letter_mapped[i])
//        {
//            return _letter_map[i];
//        }
//        else
//        {
//            int least = std::countr_zero(_letter_mapped);
//            return _letter_map[i] = static_cast<AlphabetIndex>(least);
//        }
//    }
//
//#if 0
//	/// Returns the inverse of the permutation.
//	CodewordPermutation inverse() const
//	{
//		CodewordPermutation ret(_rules);
//		for (int i = 0; i < _rules.pegs(); ++i)
//		{
//			if (peg[i] >= 0)
//				ret.peg[peg[i]] = i;
//		}
//		for (int i = 0; i < _rules.colors(); ++i)
//		{
//			if (color[i] >= 0)
//				ret.color[color[i]] = i;
//		}
//		return ret;
//	}
//#endif

//	/// Permutes the pegs and colors in a codeword.
//	Codeword permute(const Codeword &w) const
//	{
//		Codeword ret;
//		for (int i = 0; i < MM_MAX_PEGS && w[i] != Codeword::EmptyColor; ++i)
//		{
//			ret.set(peg[i], color[w[i]]);
//		}
//		return ret;
//	}
//
//	/// Permutes the pegs in a codeword.
//	Codeword permute_pegs(const Codeword &w) const
//	{
//		Codeword ret;
//		for (int i = 0; i < MM_MAX_PEGS /* && w[i] != 0xFF */; ++i)
//			ret.set(peg[i], w[i]);
//		return ret;
//	}
//
//#if 0
//	Codeword permute_colors(const Codeword &w) const
//	{
//		Codeword ret;
//		for (int i = 0; i < _rules.pegs(); ++i)
//			ret.set(i, color[w[i]]);
//		return ret;
//	}
//#endif
//};

/// Represents the set of all automorphisms (permutation of letters and
/// positions) that map a guess sequence to itself.
class AutomorphismGroup
{
public:
    /// Creates the automorphism group for an empty guess sequence.
    explicit AutomorphismGroup(const CodewordRules &rules)
    {
        PositionPermutation position_perm(rules);
        AlphabetPermutation alphabet_perm(rules); // fully unrestricted
        do
        {
            _perms.push_back(std::make_pair(position_perm, alphabet_perm));
        }
        while (position_perm.next());
    }

    /// Returns `true` if this automorphism group contains exactly one
    /// automorphism, which is necessarily the identity morphism.
    constexpr bool is_singleton() const noexcept
    {
        assert(!_perms.empty());
        if (_perms.size() == 1)
        {
            assert(_perms[0].first.is_identity());
            if (_perms[0].second.is_permutation())
            {
                assert(_perms[0].second.is_identity());
                return true;
            }
        }
        return false;
    }

    /// If guess is canonical, refine the automorphism.
    /// TODO: Some of the code branches might be optimized by
    /// TODO: studying the structure of automorphism.
    bool refine(Codeword guess) noexcept
    {
        // TODO: support Codeword::begin() and Codeword::end()
        const LetterSequence letters = guess.letters();

        // The guess is canonical if and only if for every automorphism
        // in the group, the mapped `guess` is greater than or equal to
        // `guess` in lexical order.

        std::vector<std::pair<PositionPermutation, AlphabetPermutation>> perms;

        for (const std::pair<PositionPermutation, AlphabetPermutation> &perm : _perms)
        {
            PositionPermutation position_perm = perm.first;
            AlphabetPermutation alphabet_perm = perm.second;

            // To check whether h := perm(g) < g in lexical order,
            // we compare h[j] vs g[j] for each j = 0, ..., m-1.
            // Note that h[j] = letter_map[g[inverse_position_map[j]].

            bool is_automorphism = true;
            const PositionSize m = position_perm.size();
            for (PositionIndex j = 0; j < m && is_automorphism; j++)
            {
                AlphabetIndex index = letters[j];
                PositionIndex pos = position_perm.inv(j);
                AlphabetIndex image = alphabet_perm.map(letters[pos]);
                if (image == AlphabetPermutation::NotMapped)
                {
                    image = alphabet_perm.size();
                    alphabet_perm.map(letters[pos], image);
                }
                if (image < index) // guess is not canonical
                    return false;
                if (image > index) // not an automorphism
                    is_automorphism = false;
            }

            if (is_automorphism)
                perms.push_back(std::make_pair(position_perm, alphabet_perm));
        }

        std::swap(perms, _perms);
        return true;
    }

private:
    std::vector<std::pair<PositionPermutation, AlphabetPermutation>> _perms;
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

inline void get_canonical_guesses(const CodewordRules &rules)
{
    CodewordPopulation population(rules);
    std::vector<Codeword> candidate_guesses(population.begin(), population.end());

    AutomorphismGroup group(rules);
    AutomorphismGroup current(group);
#if 1
    for (Codeword guess : candidate_guesses)
    {
        if (current.refine(guess))
        {
            std::cout << "Canonical: " << guess << std::endl;
            current = group;
        }
    }
#else
    Codeword guess(rules);
    std::istringstream("0111") >> guess;
    if (current.refine(guess))
    {
        std::cout << "Canonical: " << guess << std::endl;
        current = group;
    }
    else
    {
        std::cout << "Not canonical" << std::endl;
    }
#endif
}

inline std::vector<std::pair<Codeword, AutomorphismGroup>>
get_canonical_guesses(const AutomorphismGroup &group,
                      const CodewordRules &rules)
{
    // TODO: add Codeword::enumerate()
    CodewordPopulation population(rules);
    std::vector<Codeword> candidate_guesses(population.begin(), population.end());

    std::vector<std::pair<Codeword, AutomorphismGroup>> result;

    AutomorphismGroup current(group);
    for (Codeword guess : candidate_guesses)
    {
        if (current.refine(guess))
        {
            result.emplace_back(guess, current);
            std::cout << "Canonical: " << guess << std::endl;
            current = group;
        }
    }
    return result;
}

} // namespace mastermind

