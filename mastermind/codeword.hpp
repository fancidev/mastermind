// codeword.hpp -- basic data types used in the Mastermind game

#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <regex>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

// ============================================================================
// Constants
// ============================================================================

/// Maximum supported alphabet size times codeword length.  This is a hard
/// limit imposed by the internal storage format of the Codeword class.
#define MAX_ALPHABET_SIZE_X_CODEWORD_LENGTH 64

/// Maximum number of letters supported in the alphabet.  This is a soft
/// limit that restricts the solver to practical problem scale.
#define MAX_ALPHABET_SIZE 10

/// Maximum number of letters suppoorted in a codeword.  This is a soft
/// limit that restricts solver to practical problem scale.
#define MAX_CODEWORD_LENGTH 6

/// Helper macro that converts its argument to a string.
#define STRINGIFY(x) _STRINGIFY(x)
#define _STRINGIFY(x) #x

namespace mastermind {

/// Typedef to represent the number of letters in an alphabet.
typedef size_t AlphabetSize;

/// Typedef to represent an index into a letter in an alphabet.
typedef size_t AlphabetIndex;

/// Typedef to represent the number of letters in a codeword.
typedef size_t PositionSize;

/// Typedef to represent an index into a letter in a codeword.
typedef size_t PositionIndex;

/// Defines any particular structure among the letters in a codeword.
enum class CodewordStructure : uint8_t
{
    /// No particular structure
    None,
    /// Each letter in the codeword appears exactly once in that codeword.
    Heterogram,
};

constexpr const char *to_string(const CodewordStructure &structure) noexcept
{
    switch (structure)
    {
        case CodewordStructure::None:
            return "None";
        case CodewordStructure::Heterogram:
            return "Heterogram";
        default:
            return "Invalid";
    }
}

/// Defines a set of rules that the secret codeword conforms to.
class CodewordRules // TODO: rename to CodewordFormat? CodewordForm?
{
    /// Number of letters in the alphabet.
    uint8_t _alphabet_size;

    /// Number of letters in the codeword.
    uint8_t _codeword_length;

    /// Structual restriction of the codeword.
    CodewordStructure _structure;

public:
    /// Creates a rule set with standard Mastermind rules.
    constexpr CodewordRules() noexcept
      : _alphabet_size(6),
        _codeword_length(4),
        _structure(CodewordStructure::None)
    {
    }

    /// Creates a rule set with the given parameters.
    /// Throws `std::invalid_argument` if any argument is out of the
    /// supported range, or if no codeword could satisfy the rules.
    constexpr CodewordRules(AlphabetSize alphabet_size,
                            PositionSize codeword_length,
                            CodewordStructure structure)
      : _alphabet_size(alphabet_size),
        _codeword_length(codeword_length),
        _structure(structure)
    {
        if (!(alphabet_size >= 1))
            throw std::invalid_argument("alphabet size must be positive");
        if (!(alphabet_size <= MAX_ALPHABET_SIZE))
            throw std::invalid_argument("alphabet size must not exceed "
                                        STRINGIFY(MAX_ALPHABET_SIZE));
        if (!(codeword_length >= 1))
            throw std::invalid_argument("codeword length must be positive");
        if (!(codeword_length <= MAX_CODEWORD_LENGTH))
            throw std::invalid_argument("codeword length must not exceed "
                                        STRINGIFY(MAX_CODEWORD_LENGTH));
    }

    /// Returns the number of letters in the alphabet.
    constexpr AlphabetSize alphabet_size() const noexcept
    {
        return _alphabet_size;
    }

    /// Returns the number of letters in the codeword.
    constexpr PositionSize codeword_length() const noexcept
    {
        return _codeword_length;
    }

    /// Returns the structure among the letters of the secret.
    constexpr CodewordStructure structure() const noexcept
    {
        return _structure;
    }

    // TODO: move to COdewordPopulation class
    /// Returns the number of codewords that conform to this set of rules.
    /// This is equal to `alphabet_size ** codeword_length` for general
    /// codewords, and `P(alphabet_size, codeword_length)` for heterograms,
    /// where `P(n, m) := n*(n-1)*...*(n-m+1)`.
    constexpr size_t population_size() const noexcept
    {
        size_t count = 1;
        size_t n = _alphabet_size;
        for (PositionIndex j = 0; j < _codeword_length; j++)
        {
            count *= n;
            if (_structure == CodewordStructure::Heterogram)
                --n;
        }
        return count;
    }
//
//    /// Returns the string representation of the rule set in the form
//    /// "p4c6r" or "p4c10n".
//    std::string to_str() const
//    {
//        std::ostringstream os;
//        os << 'p' << num_pegs() << 'c' << num_colors()
//           << (repeatable()? 'r' : 'n');
//        return os.str();
//    }

//    /// Creates a rule set from string input in the form "p4c6r" or
//    /// "p4c10n".  Throws `std::invalid_argument` if the input string
//    /// is malformed or if the rules are invalid.
//    static Rules from_str(const std::string &s)
//    {
//        std::regex re("p([0-9])c([0-9][0-9]?)([nr])", std::regex::icase);
//        std::smatch m;
//        if (!std::regex_match(s, m, re))
//            throw std::invalid_argument("malformed rules string");
//
//        int num_pegs = std::stoi(m[1].str());
//        int num_colors = std::stoi(m[2].str());
//        bool repeatable = (m[3].str() == "r") || (m[3].str() == "R");
//        return Rules(num_pegs, num_colors, repeatable);
//    }
};

/// Helper function used by the Feedback class to generate the mapping
/// from feedback ordinals to feedback outcomes (a, b) for codewords of
/// length M.
template <size_t M>
static constexpr std::array<std::pair<size_t, size_t>, (M+1)*(M+2)/2>
generate_outcomes() noexcept
{
    std::array<std::pair<size_t, size_t>, (M+1)*(M+2)/2> outcomes;
    for (size_t ab = 0; ab <= M; ++ab)
    {
        for (size_t a = 0; a <= ab; ++a)
        {
            size_t b = ab - a;
            size_t k = (ab + 1) * ab / 2 + a;
            outcomes[k] = std::pair<size_t, size_t>(a, b);
        }
    }
    return outcomes;
}

/**
 * Represents the feedback from comparing two codewords.
 *
 * Given codewords x = (x[1], ..., x[M]) and y = (y[1], ..., y[M]) drawn
 * from alphabet A = (A[1], ..., A[N]), feedback is represented by the pair
 * (a, b), where
 *
 *   a := #{ j | 1 <= j <= M, x[j] == y[j] }
 *   b := sum_{i=1}^N min(#{ j | x[j] == A[i] }, #{ j | y[j] == A[i] }) - a
 *
 * In words, `a` represents the number matching letters in matching positions,
 * and `b` represents the number of matching colors in different positions.
 *
 * Note that x and y are commutative.
 *
 * Not all feedback pairs are valid for a given set of rules.  For example:
 *
 *   - (5, 0) is invalid in a game with 10 letters over 4 positions;
 *   - (3, 1) is invalid in a game with 10 letters over 4 positions;
 *   - (2, 0) is invalid in a game with 5 letters over 4 positions where
 *     each letter appears at most once.
 *
 * In addition, as constraints are added to the game, the set of valid
 * feedbacks gets smaller.
 *
 * Therefore, this class does not check for the validity of the feedback.
 *
 * To improve memory locality, the set of all feedbacks are arranged
 * in a triangle, such that `a` and `b` along each diagonal sum to the
 * same value, like below:
 *
 *   0,0  0,1  0,2  0,3  0,4
 *   1,0  1,1  1,2  1,3
 *   2,0  2,1  2,2
 *   3,0  3,1
 *   4,0
 *
 * The feedbacks are then mapped to ordinals in increasing (a+b) followed
 * by increasing (a) order, like below:
 *
 *   0,0 -> 0,1 -> 1,0 -> 0,2 -> 1,1 -> 2,0 -> ...
 *
 * To convert a feedback (a, b) to its ordinal position, use the formula
 *
 *   (a+b)*(a+b+1)/2+a
 *
 * The largest feedback ordinal for a rule set M positions is
 *
 *   M*(M+1)/2+M
 *
 * and the number of (different) feedbacks is
 *
 *   M*(M+1)/2+M+1 = (M+1)*(M+2)/2
 *
 * Therefore 8-bit storage is able to represent any feedback up to
 * M = 21 (corresponding to maximum ordinal value 252).  This is
 * sufficient for our use.
 *
 * There is no simple formula to convert feedback ordinal to (a, b).
 * We build a lookup table for that purpose.
 */
class Feedback
{
public:
    /// Type of feedback ordinal.
    typedef uint8_t ordinal_type;

    /// Maximum number of formally valid distinct feedbacks.
    static constexpr size_t MaxOutcomes =
        (MAX_CODEWORD_LENGTH + 1) * (MAX_CODEWORD_LENGTH + 2) / 2;

    /// Creates a feedback from its ordinal.
    constexpr explicit Feedback(ordinal_type ordinal) noexcept
      : _ordinal(ordinal) {}

    /// Creates a feedback with the given a and b.
    /// Throws std::invalid_argument if any argument is out of the supported
    /// range.
    constexpr Feedback(size_t a, size_t b)
    {
        if (!(a >= 0 && a <= MAX_CODEWORD_LENGTH))
            throw std::invalid_argument("a must be between 0 and "
                                        STRINGIFY(MAX_CODEWORD_LENGTH));
        if (!(b >= 0 && b <= MAX_CODEWORD_LENGTH))
            throw std::invalid_argument("b must be between 0 and "
                                        STRINGIFY(MAX_CODEWORD_LENGTH));
        size_t ab = a + b;
        if (!(ab <= MAX_CODEWORD_LENGTH))
            throw std::invalid_argument("(a + b) must be between 0 and "
                                        STRINGIFY(MAX_CODEWORD_LENGTH));

        _ordinal = static_cast<ordinal_type>((ab+1)*ab/2+a);
    }

    /// Returns the ordinal of the feedback.
    constexpr ordinal_type ordinal() const noexcept { return _ordinal; }

    /// Returns the number of matching letters in matching positions.
    constexpr size_t a() const noexcept { return _outcomes[_ordinal].first; }

    /// Returns the number of matching letters in unmatched positions.
    constexpr size_t b() const noexcept { return _outcomes[_ordinal].second; }

    /// Returns the feedback corresponding to a perfect match for a
    /// given rule set.
    static constexpr Feedback perfect_match(const CodewordRules &rules) noexcept
    {
        return Feedback(rules.codeword_length(), 0);
    }

    /// Returns the string representation of the feedback in the form "1A2B".
    template <class CharT = char,
              class Traits = std::char_traits<CharT>,
              class Allocator = std::allocator<CharT>>
    std::basic_string<CharT, Traits, Allocator> to_string() const
    {
        std::basic_ostringstream<CharT, Traits, Allocator> os;
        os << a() << 'A' << b() << 'B';
        return os.str();
    }

//    /// Creates a feedback from an input string in the form "1A2B".
//    /// Throws std::invalid_argument if the string is malformed, or
//    /// if the resulting feedback is formally invalid.
//    static Feedback from_str(const std::string &s)
//    {
//        std::regex re("([0-9])A([0-9])B", std::regex::icase);
//        std::smatch m;
//        if (!std::regex_match(s, m, re))
//            throw std::invalid_argument("malformed feedback string");
//
//        int a = std::stoi(m[1].str());
//        int b = std::stoi(m[2].str());
//        return Feedback(a, b);
//    }

    /// Default comparison operators by member comparison.
    constexpr bool operator <=>(const Feedback &other) const noexcept = default;

private:
    /// Ordinal of the feedback.
    ordinal_type _ordinal;

    static constexpr std::array<std::pair<size_t, size_t>, MaxOutcomes>
        _outcomes = generate_outcomes<MAX_CODEWORD_LENGTH>();
};

/// Represents a codeword (such as 2587).
class Codeword
{
public:
      /// Creates a codeword with the given letters.
    constexpr Codeword(AlphabetSize n,
                       std::span<AlphabetIndex> letters) noexcept
      : _position_mask(0), _alphabet_mask(0)
    {
        PositionSize m = letters.size();
        assert(m >= 1 && m <= MAX_CODEWORD_LENGTH);
        assert(n >= 1 && n <= MAX_ALPHABET_SIZE);
        assert(m * n <= MAX_ALPHABET_SIZE_X_CODEWORD_LENGTH);

        uint64_t unity = 0;
        for (PositionIndex j = 0; j < m; j++)
        {
            AlphabetIndex i = letters[j];
            assert(i >= 0 && i < n);

            _position_mask |= uint64_t(1) << (j * n + i);
            unity <<= n;
            unity |= uint64_t(1);
        }

        for (AlphabetIndex i = 0; i < n; i++)
        {
            int times = std::popcount(_position_mask & (unity << i));
            uint64_t times_mask = (uint64_t(1) << (times * n)) - uint64_t(1);
            _alphabet_mask |= (unity & times_mask) << i;
        }
    }

    /// Returns the length of the codeword (in number of letters).
    constexpr PositionSize length() const noexcept
    {
        return std::popcount(_position_mask);
    }

    /// Gets the letter at the given position.
    constexpr AlphabetIndex get(PositionIndex j, AlphabetSize n) const noexcept
    {
        assert((_position_mask >> (j * n)) != 0);
        return std::countr_zero(_position_mask >> (j * n));
    }

    /// Returns the number of occurrence of the given letter.
    constexpr size_t get_frequency(AlphabetIndex i,
                                   AlphabetSize n) const noexcept
    {
        assert(i >= 0 && i < n);
        size_t frequency = 0;
        uint64_t mask = _alphabet_mask >> i;
        while (mask & uint64_t(1))
        {
            frequency++;
            mask >>= n;
        }
        return frequency;
    }

    /// Returns `true` if any letter in the codeword appears exactly once
    /// in the codeword.
    constexpr bool is_heterogram(AlphabetSize n) const noexcept
    {
        return (_alphabet_mask >> n) == 0;
    }

//    template <class CharT,
//              class Traits = std::char_traits<CharT>,
//              class Allocator = std::allocator<CharT>>
//    std::basic_string<CharT, Traits, Allocator>
//    to_string(std::basic_string_view<CharT, Traits> alphabet) const
//    {
//        const AlphabetSize n = alphabet.size();
//        const PositionSize m = length();
//        std::basic_string<CharT, Traits, Allocator> s(m, CharT(0));
//        for (PositionIndex j = 0; j < m; j++)
//            s[j] = alphabet[get(j, n)];
//        return s;
//    }

    std::string to_string(std::string_view alphabet) const
    {
        const AlphabetSize n = alphabet.size();
        const PositionSize m = length();
        std::string s(m, '\0');
        for (PositionIndex j = 0; j < m; j++)
            s[j] = alphabet[get(j, n)];
        return s;
    }

    static constexpr Codeword from_string(std::string_view s,
                                          std::string_view alphabet)
    {
        if (!(s.size() > 0))
            throw std::invalid_argument("codeword must not be empty");
        if (!(s.size() <= MAX_CODEWORD_LENGTH))
            throw std::invalid_argument("codeword length must not exceed "
                                        STRINGIFY(MAX_CODEWORD_LENGTH));
        const PositionSize m = s.size();

        if (!(alphabet.size() > 0))
            throw std::invalid_argument("alphabet must not be empty");
        if (!(alphabet.size() <= MAX_ALPHABET_SIZE))
            throw std::invalid_argument("alphabet size must not exceed "
                                        STRINGIFY(MAX_ALPHABET_SIZE));
        const AlphabetSize n = alphabet.size();

        if (!(n * m <= MAX_ALPHABET_SIZE_X_CODEWORD_LENGTH))
            throw std::invalid_argument(
                "codeword length times alphabet size must not exceed "
                STRINGIFY(MAX_ALPHABET_SIZE_X_CODEWORD_LENGTH));

        std::array<AlphabetIndex, MAX_CODEWORD_LENGTH> letters;
        for (PositionIndex j = 0; j < m; j++)
        {
            size_t i = alphabet.find(s[j]);
            if (i == std::string::npos)
                throw std::invalid_argument("codeword contains letters "
                                            "not in the alphabet");
            letters[j] = static_cast<AlphabetIndex>(i);
        }
        return Codeword(n, std::span<AlphabetIndex>(letters).first(m));
    }

    /// Default == and != operators by member comparison.
    constexpr bool operator ==(const Codeword &other) const noexcept = default;

    /// Compares two codewords and returns the feedback.
    friend constexpr Feedback compare(const Codeword &x,
                                      const Codeword &y) noexcept
    {
        int ab = std::popcount(x._alphabet_mask & y._alphabet_mask);
        int a = std::popcount(x._position_mask & y._position_mask);
        int ordinal = (ab+1)*ab/2+a;
        return Feedback(static_cast<Feedback::ordinal_type>(ordinal));
    }

private:
    /// Bitset consisting of M groups of N bits, where bit i in group j
    /// is set if the letter at position j is the i-th letter in the
    /// alphabet.  (N := alphabet_size, M := codeword_length)
    uint64_t _position_mask;

    /// Bitset consisting of M groups of N bits, where bit i in group j
    /// is set if the i-th letter in the alphabet occurs at least (j + 1)
    /// times in the codeword.  (N := alphabet_size, M := codeword_length)
    uint64_t _alphabet_mask;
};

} // namespace mastermind
