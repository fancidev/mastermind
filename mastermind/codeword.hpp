// codeword.hpp -- basic data types used in the Mastermind game

#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <numeric>
#include <regex>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

// ============================================================================
// Constants
// ============================================================================

/// Maximum supported alphabet size times codeword length.  This value
/// must be less than the number of bits in Codeword::mask_type.
#define MAX_ALPHABET_SIZE_X_CODEWORD_LENGTH 63

/// Maximum number of letters supported in the alphabet.  This is a soft
/// limit that restricts the solver to practical problem scale.
#define MAX_ALPHABET_SIZE 10

/// Maximum number of letters supported in a codeword.  This is a soft
/// limit that restricts solver to practical problem scale.
#define MAX_CODEWORD_LENGTH 6

/// The alphabet used to format codewords as strings.
#define ALPHABET "1234567890"

static_assert(sizeof(ALPHABET) / sizeof(ALPHABET[0]) - 1 == MAX_ALPHABET_SIZE,
              "ALPHABET must contain MAX_ALPHABET_SIZE characters");

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

/// Defines a set of rules that secrets and guesses must conform to.
class CodewordRules // TODO: rename to CodewordFormat? CodewordForm? CodewordAttributes?
{
public:
    /// Creates a rule set with the given parameters.
    ///
    /// Throws `std::invalid_argument` if any of the following requirements
    /// are not satisfied:
    ///
    /// - `1 <= alphabet_size <= MAX_ALPHABET_SIZE`
    /// - `1 <= codeword_length <= MAX_CODEWORD_LENGTH`
    /// - `alphabet_size*codeword_length <= MAX_ALPHABET_SIZE_X_CODEWORD_LENGTH`
    ///
    constexpr CodewordRules(AlphabetSize alphabet_size,
                            PositionSize codeword_length,
                            CodewordStructure structure)
      : _alphabet_size(alphabet_size),
        _codeword_length(codeword_length),
        _structure(structure)
    {
        if (!(alphabet_size >= 1 && alphabet_size <= MAX_ALPHABET_SIZE))
            throw std::invalid_argument("alphabet size must be between 1 and "
                                        STRINGIFY(MAX_ALPHABET_SIZE));
        if (!(codeword_length >= 1 && codeword_length <= MAX_CODEWORD_LENGTH))
            throw std::invalid_argument("codeword length must be between 1 and "
                                        STRINGIFY(MAX_CODEWORD_LENGTH));
        if (!(alphabet_size * codeword_length <=
              MAX_ALPHABET_SIZE_X_CODEWORD_LENGTH))
            throw std::invalid_argument(
                "codeword length times alphabet size must not exceed "
                STRINGIFY(MAX_ALPHABET_SIZE_X_CODEWORD_LENGTH));
    }

    /// Creates a rule set with standard Mastermind rules.
    constexpr CodewordRules() noexcept
      : CodewordRules(AlphabetSize(6), PositionSize(4), CodewordStructure::None)
    {
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

    /// Returns the alphabet from which letters in this codeword are drawn.
    constexpr std::string_view alphabet() const noexcept
    {
        return std::string_view(ALPHABET, _alphabet_size);
    }

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

private:
    /// Number of letters in the alphabet.
    AlphabetSize _alphabet_size;

    /// Number of letters in the codeword.
    PositionSize _codeword_length;

    /// Structual restriction of the codeword.
    CodewordStructure _structure;
};

/// Writes the given codeword rule set to an output stream, in the form
/// "n10m4d" or "n6m4".
template <class CharT, class Traits>
std::basic_ostream<CharT, Traits> &
operator <<(std::basic_ostream<CharT, Traits> &os, const CodewordRules &rules)
{
    os << CharT('n') << rules.alphabet_size();
    os << CharT('m') << rules.codeword_length();
    if (rules.structure() == CodewordStructure::Heterogram)
        os << CharT('d');
    return os;
}

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

/// Writes the given feedback to an output stream, in the form "1A2B".
template <class CharT, class Traits>
std::basic_ostream<CharT, Traits> &
operator <<(std::basic_ostream<CharT, Traits> &os, const Feedback &feedback)
{
    return os << feedback.a() << 'A' << feedback.b() << 'B';
}

/// Returns a bit pattern of type T
template <class T>
constexpr T cyclic_mask(size_t bits_in_group, size_t num_groups)
{
    T mask(0);
    for (size_t j = 0; j < num_groups; j++)
    {
        mask <<= bits_in_group;
        mask |= 1;
    }
    return mask;
}

/// Represents a codeword (such as 2587).
class Codeword
{
public:
//    /// Creates an empty codeword.
//    constexpr Codeword() noexcept
//      : Codeword(AlphabetSize(0), std::span<AlphabetIndex>()) {}

    /// Creates a codeword with the given letters under the given rules.
    ///
    /// Throws `std::invalid_argument` if letters do not conform to the rules.
    constexpr Codeword(std::span<AlphabetIndex> letters,
                       const CodewordRules &rules)
      : _position_mask(0), _alphabet_mask(0)
    {
        const AlphabetSize n = rules.alphabet_size();
        const PositionSize m = rules.codeword_length();
        if (letters.size() != m)
            throw std::invalid_argument("letters length does not conform");

        for (PositionIndex j = 0; j < m; j++)
        {
            AlphabetIndex i = letters[j];
            if (!(i >= 0 && i < n))
                throw std::invalid_argument("letter out of range");
            _position_mask |= mask_type(1) << (j * n + i);
        }

        mask_type unity = cyclic_mask<mask_type>(n, m);
        for (AlphabetIndex i = 0; i < n; i++)
        {
            int times = std::popcount(_position_mask & (unity << i));
            mask_type times_mask = (mask_type(1) << (times * n)) - mask_type(1);
            _alphabet_mask |= (unity & times_mask) << i;
        }

        _alphabet_mask |= mask_type(1) << (m * n);

        if (rules.structure() == CodewordStructure::Heterogram &&
            !is_heterogram())
            throw std::invalid_argument("heterogram expected");
    }

    /// Returns the length of the codeword (in number of letters).
    constexpr PositionSize length() const noexcept
    {
        return std::popcount(_position_mask);
    }

    /// Returns the size of the alphabet from which the codeword is drawn.
    constexpr AlphabetSize alphabet_size() const noexcept
    {
        int m = std::popcount(_position_mask);
        int mn = MaskBits - 1 - std::countl_zero(_alphabet_mask);
        return mn / m;
    }

    /// Gets the letter at the given position.
    constexpr AlphabetIndex get(PositionIndex j) const noexcept
    {
        assert(j >= 0 && j < length());
        const AlphabetSize n = alphabet_size();
        assert((_position_mask >> (j * n)) != 0);
        return std::countr_zero(_position_mask >> (j * n));
    }

    /// Returns the number of times the i-th letter of the alphabet appears
    /// in the codeword.
    constexpr size_t get_frequency(AlphabetIndex i) const noexcept
    {
        const PositionSize m = length();
        const AlphabetSize n = alphabet_size();
        assert(i >= 0 && i < n);

        const mask_type mask = cyclic_mask<mask_type>(n, m);
        return std::popcount(_alphabet_mask & (mask << i));
    }

    /// Returns `true` if every letter in the codeword appears exactly once
    /// in the codeword.  (Note: returns `true` for an empty codeword.)
    constexpr bool is_heterogram() const noexcept
    {
        return std::has_single_bit(_alphabet_mask >> alphabet_size());
    }

    /// Returns true if this codeword conforms to the given rules.
    constexpr bool conforms_to(const CodewordRules &rules) const noexcept
    {
        return (length() == rules.codeword_length())
            && (alphabet_size() == rules.alphabet_size())
            && (is_heterogram() || rules.structure() == CodewordStructure::None);
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

    template <class CharT = char,
              class Traits = std::char_traits<CharT>,
              class Allocator = std::allocator<CharT>>
    std::basic_string<CharT, Traits, Allocator> to_string() const
    {
        PositionSize m = length();
        std::basic_string<CharT, Traits, Allocator> s(m, CharT(0));
        for (PositionIndex j = 0; j < m; j++)
            s[j] = ALPHABET[get(j)];
        return s;
    }

    static constexpr Codeword from_string(std::string_view s,
                                          const CodewordRules &rules)
    {
        const PositionSize m = rules.codeword_length();
        if (s.size() != m)
            throw std::invalid_argument("codeword length does not conform");

        std::array<AlphabetIndex, MAX_CODEWORD_LENGTH> letters;
        for (PositionIndex j = 0; j < m; j++)
        {
            size_t i = rules.alphabet().find(s[j]);
            if (i == std::string::npos)
                throw std::invalid_argument("codeword contains letters "
                                            "not in the alphabet");
            letters[j] = static_cast<AlphabetIndex>(i);
        }
        return Codeword(std::span(letters).first(m), rules);
    }

    /// Default == and != operators by member comparison.
    constexpr bool operator ==(const Codeword &other) const noexcept = default;

    /// Compares this codeword to another and returns the feedback.
    friend constexpr Feedback compare(const Codeword &x,
                                      const Codeword &y) noexcept
    {
        assert(x.length() == y.length());
        assert(x.alphabet_size() == y.alphabet_size());
        int ab = std::popcount(x._alphabet_mask & y._alphabet_mask) - 1;
        int a = std::popcount(x._position_mask & y._position_mask);
        int ordinal = (ab+1)*ab/2+a;
        return Feedback(static_cast<Feedback::ordinal_type>(ordinal));
    }

private:
    /// Integral type used to store the codeword's bit-mask.
    typedef uint64_t mask_type;

    /// Number of bits mask_type;
    static constexpr size_t MaskBits = sizeof(mask_type) * 8;

    static_assert(MAX_ALPHABET_SIZE_X_CODEWORD_LENGTH < MaskBits,
                  "mask_type does not have enough bits");

    /// Bitset consisting of m groups of n bits, where bit i in group j
    /// is set if the letter at position j is the i-th letter of the
    /// alphabet.  (n := alphabet_size, m := codeword_length)
    /// The most significant (MaskBits-m*n) bits are always 0.
    mask_type _position_mask;

    /// Bitset consisting of m groups of n bits, where bit i in group j
    /// is set if the i-th letter of the alphabet occurs at least (j + 1)
    /// times in the codeword.  (n := alphabet_size, m := codeword_length)
    ///
    /// In addition, bit `m*n` is set to 1 to mark the alphabet size.
    /// The remaining most significant bits are set to 0.
    mask_type _alphabet_mask;
};

/// Writes the given codeword to an output stream, in the form "1224".
template <class CharT, class Traits>
std::basic_ostream<CharT, Traits> &
operator <<(std::basic_ostream<CharT, Traits> &os, const Codeword &codeword)
{
    PositionSize m = codeword.length();
    for (PositionIndex j = 0; j < m; j++)
        os << ALPHABET[codeword.get(j)];
    return os;
}

/// Returns a vector v with (m + 1) elements (m := codeword length), where
/// v[j] := the size of the sub-population where the first j letters in the
/// codeword are fixed.
template <size_t M>
constexpr std::array<size_t, M + 1>
sub_population_sizes(const CodewordRules &rules) noexcept
{
    std::array<size_t, M + 1> sizes {};
    size_t count = 1;
    const size_t n = rules.alphabet_size();
    const size_t m = rules.codeword_length();
    for (PositionIndex j = m; j > 0; --j)
    {
        sizes[j] = count;
        if (rules.structure() == CodewordStructure::Heterogram)
            count *= (n - j + 1);
        else
            count *= n;
    }
    sizes[0] = count;
    return sizes;
}

/// Represents the collection of all codewords conforming to a given rule set.
class CodewordPopulation
{
public:
    /// Creates the population from the given rules.
    constexpr CodewordPopulation(const CodewordRules &rules) noexcept
      : _rules(rules),
        _sub_population_sizes(sub_population_sizes<MAX_CODEWORD_LENGTH>(rules))
    {
    }

    /// Returns the rules that the population of codewords conform to.
    constexpr const CodewordRules &rules() const noexcept { return _rules; }

    /// Returns the number of (distinct) codewords in the population.
    /// (Note: The empty codeword is the unique codeword of length zero.)
    constexpr size_t size() const noexcept { return _sub_population_sizes[0]; }

    /// Returns the codeword at the given index in the population in
    /// lexicographical order.
    constexpr Codeword get(size_t index) const noexcept
    {
        assert(index >= 0 && index < size());
        const size_t n = _rules.alphabet_size();
        const size_t m = _rules.codeword_length();

        std::array<AlphabetIndex, MAX_CODEWORD_LENGTH> letters {};
        if (_rules.structure() == CodewordStructure::Heterogram)
        {
            std::array<AlphabetIndex, MAX_ALPHABET_SIZE> alphabet {};
            std::iota(alphabet.begin(), alphabet.end(), AlphabetIndex(0));
            for (PositionIndex j = 0; j < m; j++)
            {
                size_t i = index / _sub_population_sizes[j + 1];
                letters[j] = alphabet[i];
                std::copy(alphabet.cbegin() + i + 1,
                          alphabet.cbegin() + n - j,
                          alphabet.begin() + i);
                index %= _sub_population_sizes[j + 1];
            }
        }
        else
        {
            for (PositionIndex j = 0; j < m; j++)
            {
                size_t i = index / _sub_population_sizes[j + 1];
                letters[j] = static_cast<AlphabetIndex>(i);
                index %= _sub_population_sizes[j + 1];
            }
        }
        return Codeword(std::span(letters).first(m), _rules);
    }

private:
    /// Rule set that all codewords in the population conform to.
    CodewordRules _rules;

    /// Vector with (m + 1) elements (m := codeword length), where element j
    /// denotes the size of the sub-population where the first j letters of
    /// a codeword are fixed.
    std::array<size_t, MAX_CODEWORD_LENGTH + 1> _sub_population_sizes;
};

} // namespace mastermind
