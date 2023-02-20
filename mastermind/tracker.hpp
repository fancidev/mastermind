#pragma once

#include "codeword.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <numeric>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mastermind {

/// Represents a bijection from a subset of the index set {0, ..., N-1}
/// of type Index to another subset of that index set.
template <class Index, size_t N>
class Bijection : public std::array<Index, N>
{
public:
    static_assert(std::is_integral_v<Index>, "Index must be integral");
    static_assert(N < static_cast<size_t>(Index(-1)), "no room for placeholder");

    /// A placeholder that indicates an unmapped index.
    static constexpr Index not_mapped = Index(-1);

    /// Creates an empty mapping.
    constexpr Bijection() noexcept
    {
        using base = std::array<Index, N>;
        std::fill(base::begin(), base::end(), not_mapped);
    }

    /// Creates an identity mapping of the first `n` indices.
    explicit constexpr Bijection(size_t n) noexcept : Bijection()
    {
        assert(n <= N);
        using base = std::array<Index, N>;
        std::iota(base::begin(), base::begin() + n, Index(0));
    }

    /// Returns a bijection defined on the full index set {0, ..., N-1}
    /// by extending the current bijection.
    constexpr Bijection complete() const noexcept
    {
        std::array<bool, N> used {};
        for (size_t i = 0; i < N; i++)
        {
            if ((*this)[i] != not_mapped)
                used[(*this)[i]] = true;
        }
        Bijection completed {*this};
        Index next{0};
        for (size_t i = 0; i < N; i++)
        {
            if (completed[i] == not_mapped)
            {
                while (used[next]) next++;
                completed[i] = next;
                used[next++] = true;
            }
        }
        return completed;
    }

    /// Returns the inverse mapping.
    constexpr Bijection inverse() const noexcept
    {
        Bijection inv;
        for (size_t i = 0; i < N; i++)
        {
            Index ii = (*this)[i];
            if (ii != not_mapped)
                inv[ii] = static_cast<Index>(i);
        }
        return inv;
    }

    /// Returns the bijection equal to this composed with `other`
    /// (`other` is applied first).
    constexpr Bijection composed_with(const Bijection &other) const noexcept
    {
        Bijection composed;
        for (size_t i = 0; i < N; i++)
        {
            Index ii = other[i];
            if (ii != not_mapped)
            {
                Index iii = (*this)[ii];
                assert(iii != not_mapped);
                composed[i] = iii;
            }
        }
        return composed;
    }

    /// Writes the bijection to an output stream in the form "(3*5)".
    template <class CharT, class Traits>
    friend std::basic_ostream<CharT, Traits> &
    operator <<(std::basic_ostream<CharT, Traits> &os, const Bijection &p)
    {
        size_t n = N;
        while (n > 0 && p[n - 1] == not_mapped)
            n--;
        os << CharT('(');
        for (size_t j = 0; j < n; j++)
        {
            if (p[j] == not_mapped)
                os << CharT('*');
            else
                os << p[j];
        }
        os << CharT(')');
        return os;
    }
};

/// Represents a bijection from a subset of codewords to another subset of
/// codewords.
struct CodewordMorphism
{
    /// Represents a permutation of positions.
    using PositionMap = Bijection<Position, MAX_CODEWORD_SIZE>;

    /// Represents a bijection between two subsets of letters.
    using LetterMap = Bijection<Letter, MAX_ALPHABET_SIZE>;

    PositionMap position_inv;
    LetterMap letter_map;

    /// Returns the inverse of this morphism.
    constexpr CodewordMorphism inverse() const noexcept
    {
        return CodewordMorphism {
            position_inv.inverse(), letter_map.inverse()
        };
    }

    /// Returns an morphism equal to this composed with `other` (`other`
    /// is applied first).
    constexpr CodewordMorphism composed_with(const CodewordMorphism &other)
    const noexcept
    {
        return CodewordMorphism {
            other.position_inv.composed_with(position_inv),
            letter_map.composed_with(other.letter_map)
        };
    }
};

/// Tracks possible secrets as constraints are added.
class CodewordSet
{
public:
    /// Copy constructor and copy assignment.
    explicit CodewordSet(const CodewordSet &other) = default;
    CodewordSet & operator=(const CodewordSet &other) = default;

    /// Move constructor and move assignment.
    explicit CodewordSet(CodewordSet &&other) = default;
    CodewordSet & operator=(CodewordSet &&other) = default;

    /// Creates a set of all codewords conforming to the given rules.
    explicit CodewordSet(const CodewordRules &rules);

    /// Adds a constraint.
    void push_constraint(const Constraint &constraint);

    /// Returns the rules the codewords conform to.
    constexpr const CodewordRules &rules() const noexcept { return _rules; }

    /// Returns the constraints the codewords conform to.
    constexpr const std::vector<Constraint> &constraints() const noexcept
    {
        return _constraints;
    }

    constexpr std::span<const Codeword> possible_secrets() const noexcept
    {
        return _list;
    }

    constexpr const std::vector<CodewordMorphism> &morphisms() const noexcept
    {
        return _morphisms;
    }

    std::vector<Codeword> get_canonical_guesses() const;

private:
    /// Rules that codewords in the set conform to.
    CodewordRules _rules;

    /// Constraints that codewords in the set passes.
    std::vector<Constraint> _constraints;

    /// `_used[i]` is true iff letter i appeared in any of the guesses.
    std::bitset<MAX_ALPHABET_SIZE> _used;

    /// List of codewords conforming to `_rules` and passing `_constraints`.
    std::vector<Codeword> _list;

    /// List of all codeword morphisms that map the guess sequence of
    /// `constraints` to a canonical sequence.  Only the `_used` letters
    /// of the constraints are mapped in the `letter_map` of these
    /// morphisms.
    std::vector<CodewordMorphism> _morphisms;
};

} // namespace mastermind
