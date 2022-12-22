#pragma once

#include <cstddef>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>

#ifndef MM_MAX_PEGS
/// The maximum number of pegs supported by the program.
/// This value must be smaller than or equal to @c 9.
#define MM_MAX_PEGS 6
#endif

#if MM_MAX_PEGS > 9
#error MM_MAX_PEGS must be smaller than or equal to 9.
#endif

#ifndef MM_MAX_COLORS
/// The maximum number of colors supported by the program.
/// This value must be smaller than or equal to @c 10.
#define MM_MAX_COLORS 10
#endif

#if MM_MAX_COLORS > 10
#error MM_MAX_COLORS must be smaller than or equal to 10.
#endif

#if (MM_MAX_PEGS + MM_MAX_COLORS) != 16
#error MM_MAX_PEGS and MM_MAX_COLORS must add to 16.
#endif

namespace mastermind
{
    /// Defines the rules of a game.
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

        /// Returns the number of (distinct) conforming codewords for this
        /// rule set.  This is equal to `c ** p` for non-repeatable rules
        /// and `P(c,p) := c*(c-1)*...*(c-p+1)` for repeatable rules (where
        /// c := num_colors and p := num_pegs).
        size_t num_admissible() const
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

    };

    /// Writes the string representation of rule set `r` to output, in the
    /// form "p4c6r" or "p4c10n".  No delimiter is output.
    template<class CharT, class Traits>
    std::basic_ostream<CharT, Traits> &
    operator <<(std::basic_ostream<CharT, Traits> &os, const Rules &r)
    {
        os << 'p' << r.num_pegs() << 'c' << r.num_colors()
           << (r.repeatable()? 'r' : 'n');
        return os;
    }

    /// Reads rule set `r` from string input, in the form "p4c6r" or "p4c10n",
    /// followed by a whitespace or end-of-input.
    /// Throws `std::invalid_argument` if the input is malformed or if the
    /// resulting rules are invalid.
    template<class CharT, class Traits>
    std::basic_istream<CharT, Traits> &
    operator >>(std::basic_istream<CharT, Traits> &is, Rules &r)
    {
        std::string s; // TODO: should use template parameter
        is >> s;

        std::regex re("p([0-9])c([0-9][0-9]?)([nr])", std::regex::icase);
        std::smatch m;
        if (!std::regex_match(s, m, re))
            throw std::invalid_argument("malformed rules string");

        int num_pegs = std::stoi(m[1].str());
        int num_colors = std::stoi(m[2].str());
        bool repeatable = (m[3].str() == "r") || (m[3].str() == "R");
        r = Rules(num_pegs, num_colors, repeatable);

        return is;
    }

//    /// @cond DETAILS
//    namespace details {
//
//    struct RulesFormatter
//    {
//        Rules rules;
//
//        RulesFormatter(const Rules &r) : rules(r) { }
//
//        /// Returns the index to a custom ios format field.
//        static int index()
//        {
//            static int i = std::ios_base::xalloc();
//            return i;
//        }
//    };
//
//    inline std::istream& operator >> (std::istream &s, const RulesFormatter &f)
//    {
//        s.iword(RulesFormatter::index()) = f.rules.pack();
//        return s;
//    }
//
//    } // namespace details
//    /// @endcond
//
//    /// Sets the rules to be associated with a stream.
//    /// @ingroup Rules
//    inline details::RulesFormatter setrules(const Rules &rules)
//    {
//        return details::RulesFormatter(rules);
//    }
//
//    /// Gets the rules associated with a stream.
//    /// @ingroup Rules
//    inline Rules getrules(std::ios_base &s)
//    {
//        return Rules::unpack(s.iword(details::RulesFormatter::index()));
//    }

} // namespace mastermind
