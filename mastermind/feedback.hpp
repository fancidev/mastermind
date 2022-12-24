#pragma once

#include "rules.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <regex>
#include <string>

/// Generates the mapping from feedback ordinal to feedback outcome (a, b)
/// for codewords of length M.
template <size_t M>
static constexpr std::array<std::pair<size_t, size_t>, (M+1)*(M+2)/2>
generate_outcomes()
{
    std::array<std::pair<size_t, size_t>, (M+1)*(M+2)/2> outcomes;
    for (size_t nAB = 0; nAB <= M; ++nAB)
    {
        for (size_t nA = 0; nA <= nAB; ++nA)
        {
            size_t nB = nAB - nA;
            size_t k = (nAB+1)*nAB/2+nA;
            outcomes[k] = std::pair<size_t, size_t>(nA, nB);
        }
    }
    return outcomes;
}

namespace mastermind
{
    /**
     * Represents the feedback from comparing a *guess* to a *secret*.
     *
     * A feedback takes the form `xAyB`, where `x` represents the number
     * of colors in the right peg, and `y` represents the number of colors
     * in the wrong peg.
     *
     * Not all *formally* valid feedbacks are *actually* valid for a given
     * set of rules.  For example:
     *
     *   - "5A0B" is invalid in a game with 10 colors and 4 pegs;
     *   - "3A1B" is invalid in a game with 10 colors and 4 pegs;
     *   - "2A0B" is invalid in a game with 5 colors, 4 pegs and no repetition.
     *
     * In addition, as constraints are added to the game, the set of valid
     * feedbacks gets smaller.
     *
     * Therefore, this class only checks for *formal validity* but not for
     * *actual* validity.
     *
     * To improve memory locality, formally valid feedbacks are arranged
     * in a triangle, such that `x` and `y` along each diagonal sum to the
     * same value, like below:
     *
     *   0A0B  1A0B  2A0B  3A0B  4A0B
     *   0A1B  1A1B  2A1B  3A1B
     *   0A2B  1A2B  2A2B
     *   0A3B  1A3B
     *   0A4B
     *
     * The feedbacks are then mapped to ordinals in the following order:
     *
     *   0A0B -> 0A1B -> 1A0B -> 0A2B -> 1A1B -> 2A0B -> ...
     *
     * To convert a feedback `xAyB` to its ordinal position, use the formula
     *
     *   (x+y)*(x+y+1)/2+x
     *
     * The largest feedback ordinal for a rule set with `p` pegs is
     *
     *   p*(p+1)/2+p = p*(p+3)/2
     *
     * Therefore 8-bit storage is able to represent any feedback up to
     * 21 pegs (corresponding to maximum ordinal value 252).  This is
     * sufficient for our use.
     *
     * There is no simple formula to convert feedback ordinal to (x, y).
     * We build a lookup table for that purpose.
     */
    class Feedback
    {
    public:

        /// Type of feedback ordinal.
        typedef uint8_t ordinal_type;

        /// Maximum number of formally valid distinct feedbacks.
        static constexpr size_t MaxOutcomes = (MM_MAX_PEGS+1)*(MM_MAX_PEGS+2)/2;

    private:

        /// Ordinal of the feedback.
        ordinal_type _ordinal;

        static constexpr std::array<std::pair<size_t, size_t>, MaxOutcomes>
            _outcomes = generate_outcomes<MM_MAX_PEGS>();

    public:

//        /// Creates an empty feedback.
//        Feedback() : _value(-1) { }

        /// Creates a feedback from its ordinal.
        /// Does not validate argument.
//        explicit Feedback(ordinal_type ordinal) : _ordinal(ordinal) { }

        /// Creates a feedback with the given `nA` and `nB`.
        /// Throws std::invalid_argument if formally invalid.
        Feedback(size_t nA, size_t nB)
        {
            if (!(nA >= 0 && nA <= MM_MAX_PEGS))
                throw std::invalid_argument("nA out of range");
            if (!(nB >= 0 && nB <= MM_MAX_PEGS))
                throw std::invalid_argument("nB out of range");
            size_t nAB = nA + nB;
            if (!(nAB <= MM_MAX_PEGS))
                throw std::invalid_argument("(nA + nB) out of range");

            _ordinal = static_cast<ordinal_type>((nAB+1)*nAB/2+nA);
        }

//        /// Creates a feedback from a string of the form "1A2B".
//        /// Throws std::invalid_argument if the string is malformed or
//        /// if the feedback is not formally valid.
//        explicit Feedback(const char *s) : _value(-1)
//        {
//            if (s &&
//                (s[0] >= '0' && s[0] <= '9') &&
//                (s[1] == 'A' || s[1] == 'a') &&
//                (s[2] >= '0' && s[2] <= '9') &&
//                (s[3] == 'B' || s[3] == 'b') &&
//                (s[4] == '\0'))
//            {
//                *this = Feedback(s[0]-'0', s[2]-'0');
//            }
//        }

        /// Returns the ordinal of the feedback.
        ordinal_type ordinal() const { return _ordinal; }

        /// Returns the number of matching colors in matching pegs.
        size_t a() const { return _outcomes[_ordinal].first; }

        /// Returns the number of matching colors in unmatched pegs.
        size_t b() const { return _outcomes[_ordinal].second; }

//        /**
//         * Compact format of a feedback.
//         * In this format, bits 0-3 stores @c nB, and bits 4-7 stores
//         * @c nA. To convert a feedback into the compact format, use
//         * the formula <code>x = (nA << 4) | nB</code>.
//         * To restore a feedback from compact format, use the formula
//         * <code>nA = x >> 4</code> and <code>nB = x & 0x0F</code>.
//         */
//        typedef unsigned char compact_type;
//
//        /// Converts the feedback into compact form.
//        compact_type pack() const
//        {
//            std::pair<int,int> ab = outcome_table::lookup(_value);
//            return (compact_type)((ab.first << 4) | ab.second);
//        }
//
//        /// Restores a feedback from compact form.
//        /// @todo check malformed input
//        static Feedback unpack(compact_type ab)
//        {
//            int nA = (ab >> 4), nB = (ab & 0x0F);
//            return Feedback(nA, nB);
//        }

        /// Returns the feedback corresponding to a perfect match for a
        /// given rule set.
        /// TODO: this function is not well-defined if codewords are not of
        /// TODO: uniform length.
        static Feedback perfect_match(const Rules &rules)
        {
            return Feedback(rules.num_pegs(), 0);
        }

//        /// Returns the size of the set of distinct feedback values under
//        /// a given set of rules. The practically impossible feedback
//        /// <code>(p-1,1)</code> is included.
//        static size_t size(const Rules &rules)
//        {
//            return (size_t)perfectValue(rules)._value + 1;
//        }
    };

    /// Tests whether two feedbacks are equal.
    inline bool operator ==(const Feedback &a, const Feedback &b)
    {
        return a.ordinal() == b.ordinal();
    }

    /// Tests whether two feedbacks are unequal.
    inline bool operator !=(const Feedback &a, const Feedback &b)
    {
        return a.ordinal() != b.ordinal();
    }

    /// Writes a feedback to an output stream in the form "1A2B".
    template<class CharT, class Traits>
    std::basic_ostream<CharT, Traits> &
    operator <<(std::basic_ostream<CharT, Traits> &os, const Feedback &feedback)
    {
        return os << feedback.a() << 'A' << feedback.b() << 'B';
    }

    /// Reads a feedback from an input stream in the form "1A2B", followed
    /// by a whitespace or end-of-input.
    template<class CharT, class Traits>
    std::basic_istream<CharT, Traits> &
    operator >>(std::basic_istream<CharT, Traits> &is, Feedback &feedback)
    {
        std::string s;
        if (!(is >> s))
            return is;

        std::regex re("([0-9])A([0-9])B", std::regex::icase);
        std::smatch m;
        if (!std::regex_match(s, m, re))
            throw std::invalid_argument("malformed rules string");

        int nA = std::stoi(m[1].str());
        int nB = std::stoi(m[2].str());
        feedback = Feedback(nA, nB);
        // is.setstate(std::ios_base::failbit);
        return is;
    }

} // namespace mastermind
