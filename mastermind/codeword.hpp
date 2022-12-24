// codeword.hpp -- basic data types used in the Mastermind game

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <regex>
#include <stdexcept>
#include <string>
#include <sstream>

//#include "util/intrinsic.hpp"

// ============================================================================
// Constants
// ============================================================================

#ifndef MM_MAX_PEGS
/// The maximum number of pegs supported by the program.
#define MM_MAX_PEGS 6
#endif

#if MM_MAX_PEGS > 9
#error MM_MAX_PEGS must be less than or equal to 9.
#endif

#ifndef MM_MAX_COLORS
/// The maximum number of colors supported by the program.
#define MM_MAX_COLORS 10
#endif

#if MM_MAX_COLORS > 10
#error MM_MAX_COLORS must be less than or equal to 10.
#endif

#if (MM_MAX_PEGS + MM_MAX_COLORS) != 16
#error MM_MAX_PEGS and MM_MAX_COLORS must sum to 16.
#endif

namespace mastermind
{
    /// Defines a set of rules that codewords must conform to.
    class Rules
    {
        /// Number of pegs in the game.
        /// Must be within the range [1, MM_MAX_PEGS].
        size_t _num_pegs;

        /// Number of colors in the game.
        /// Must be within the range [1, MM_MAX_COLORS].
        size_t _num_colors;

        /// Whether any color may appear more than once in a codeword.
        bool _repeatable;

    public:
        /// Creates a rule set with standard Mastermind rules.
        Rules() : _num_pegs(4), _num_colors(6), _repeatable(true) { }

        /// Creates a rule set with the given parameters.
        /// Throws `std::invalid_argument` if any argument is out of range
        /// or if the rules admit no codeword.
        Rules(size_t num_pegs, size_t num_colors, bool repeatable)
          : _num_pegs(num_pegs),
            _num_colors(num_colors),
            _repeatable(repeatable)
        {
            if (!(num_pegs >= 1 && num_pegs <= MM_MAX_PEGS))
                throw std::invalid_argument("num_pegs out of range");
            if (!(num_colors >= 1 && num_colors <= MM_MAX_COLORS))
                throw std::invalid_argument("num_colors out of range");
            if (!repeatable && num_colors < num_pegs)
                throw std::invalid_argument(
                    "num_colors must be greater than or equal to num_pegs if "
                    "the colors are not allowed to repeat");
        }

        /// Returns the number of pegs in the game.
        size_t num_pegs() const { return _num_pegs; }

        /// Returns the number of colors in the game.
        size_t num_colors() const { return _num_colors; }

        /// Returns true if a color is allowed to appear more than once in a
        /// codeword.
        bool repeatable() const { return _repeatable; }
        
        /// Returns the number of (distinct) codewords that conform to this
        /// rule set.  This is equal to `num_colors ** num_pegs` for
        /// non-repeatable rules and `P(num_colors, num_pegs)` for
        /// repeatable rules, where `P(n, r) := n*(n-1)*...*(n-r+1)`.
        size_t population_size() const
        {
            size_t n = 1;
            size_t c = _num_colors;
            for (size_t i = 0; i < _num_pegs; i++)
            {
                n *= c;
                c += _repeatable ? 0 : -1;
            }
            return n;
        }

        /// Returns the string representation of the rule set in the form
        /// "p4c6r" or "p4c10n".
        std::string to_str() const
        {
            std::ostringstream os;
            os << 'p' << num_pegs() << 'c' << num_colors()
               << (repeatable()? 'r' : 'n');
            return os.str();
        }
        
        /// Creates a rule set from string input in the form "p4c6r" or
        /// "p4c10n".  Throws `std::invalid_argument` if the input string
        /// is malformed or if the rules are invalid.
        static Rules from_str(const std::string &s)
        {
            std::regex re("p([0-9])c([0-9][0-9]?)([nr])", std::regex::icase);
            std::smatch m;
            if (!std::regex_match(s, m, re))
                throw std::invalid_argument("malformed rules string");

            int num_pegs = std::stoi(m[1].str());
            int num_colors = std::stoi(m[2].str());
            bool repeatable = (m[3].str() == "r") || (m[3].str() == "R");
            return Rules(num_pegs, num_colors, repeatable);
        }
    };

    /// Generates the mapping from feedback ordinal to feedback outcome (a, b)
    /// for codewords of length M.  Used by the Feedback class.
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

        /// Creates a feedback from its ordinal.
        /// Does not validate argument.
        //  explicit Feedback(ordinal_type ordinal) : _ordinal(ordinal) { }

        /// Creates a feedback with the given `a` and `b`.
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

        /// Returns the ordinal of the feedback.
        ordinal_type ordinal() const { return _ordinal; }

        /// Returns the number of matching colors in matching pegs.
        size_t a() const { return _outcomes[_ordinal].first; }

        /// Returns the number of matching colors in unmatched pegs.
        size_t b() const { return _outcomes[_ordinal].second; }

        /// Returns the feedback corresponding to a perfect match for a
        /// given rule set.
        /// TODO: this function is not well-defined if codewords are not of
        /// TODO: uniform length.
        static Feedback perfect_match(const Rules &rules) noexcept
        {
            return Feedback(rules.num_pegs(), 0);
        }
        
        /// Returns the string representation of the feedback in the form
        /// "1A2B".
        std::string to_str() const
        {
            std::ostringstream os;
            os << a() << 'A' << b() << 'B';
            return os.str();
        }
        
        /// Creates a feedback from an input string in the form "1A2B".
        /// Throws std::invalid_argument if the string is malformed, or
        /// if the resulting feedback is formally invalid.
        static Feedback from_str(const std::string &s)
        {
            std::regex re("([0-9])A([0-9])B", std::regex::icase);
            std::smatch m;
            if (!std::regex_match(s, m, re))
                throw std::invalid_argument("malformed feedback string");

            int a = std::stoi(m[1].str());
            int b = std::stoi(m[2].str());
            return Feedback(a, b);
        }
        
        /// Tests whether two feedbacks are equal.
        bool operator ==(const Feedback &other) const noexcept
        {
            return ordinal() == other.ordinal();
        }

        /// Tests whether two feedbacks are unequal.
        bool operator !=(const Feedback &other) const noexcept
        {
            return ordinal() != other.ordinal();
        }

    private:
        /// Ordinal of the feedback.
        ordinal_type _ordinal;

        static constexpr std::array<std::pair<size_t, size_t>, MaxOutcomes>
            _outcomes = generate_outcomes<MM_MAX_PEGS>();
    };

    /// Represents a codeword (such as 2587).
    class alignas(16) Codeword
    {
        /// Number of occurrence of each color.
        /// The value must be between 0 and `MM_MAX_PEGS` inclusive.
        int8_t _counter[MM_MAX_COLORS]; // TODO: rename to freq

        /// The (zero-based) color on each peg.
        /// If a peg is empty, the corresponding value is (int8_t)(-1).
        int8_t _digit[MM_MAX_PEGS]; // TODO: rename to content

    public:

        /// Constant representing an 'empty' color.
        static const int EmptyColor = -1;
            
        /// <summary>Creates an empty codeword.</summary>
        Codeword()
        {
            std::memset(_counter, 0, sizeof(_counter));
            std::memset(_digit, -1, sizeof(_digit));
        }

        /// Return `true` if the codeword is empty.
        bool IsEmpty() const { return _digit[0] < 0; }

        /// <summary>Gets the color in the given peg.</summary>
        /// <param name="peg">Zero-based index of the peg.</param>
        /// <return>Color in the given peg.</return>
        int operator [](int peg) const
        {
            assert(peg >= 0 && peg < MM_MAX_PEGS);
            return _digit[peg];
        }

        /// Sets the color on a given peg.
        void set(int peg, int color)
        {
            assert(peg >= 0 && peg < MM_MAX_PEGS);
            assert(color == -1 || (color >= 0 && color < MM_MAX_COLORS));
            if (_digit[peg] >= 0)
                --_counter[(int)_digit[peg]];
            if ((_digit[peg] = (char)color) >= 0)
                ++_counter[color];
        }

        /// Returns the frequency (number of occurrence) of the given `color`.
        int count(int color) const
        {
            assert(color >= 0 && color < MM_MAX_COLORS);
            return _counter[color];
        }

        /// Returns `true` if any color appears more than once in the codeword.
        bool has_repetition() const
        {
#if 0
            int mask = util::simd::byte_mask(
                *(util::simd::simd_t<int8_t,16>*)this > (int8_t)1);
            return (mask & ((1 << MM_MAX_COLORS) - 1)) != 0;
#else
            return std::any_of(std::begin(_counter),
                               std::end(_counter),
                               [](int8_t n) { return n > 1; });
#endif
        }

        /// Return `true` if the codeword conforms to the given rule set.
//        bool conforms_to(const Rules &rules) const;
        
        /// Tests whether two codewords are equal.
        bool operator ==(const Codeword &other) const noexcept
        {
            return memcmp(this, &other, sizeof(Codeword)) == 0;
        }

        /// Tests whether two codewords are unequal.
        bool operator !=(const Codeword &other) const noexcept
        {
            return memcmp(this, &other, sizeof(Codeword)) != 0;
        }
    };

//    /// Writes a codeword to an output stream.
//    std::ostream& operator << (std::ostream &os, const Codeword &c);
//
//    /// Reads a codeword from an input stream.
//    /// @ingroup Codeword
//    std::istream& operator >> (std::istream &is, Codeword &c);

} // namespace mastermind
