// codeword.hpp -- basic data types used in the Mastermind game

#pragma once

#include <algorithm>
#include <array>
#include <bit>
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

// ============================================================================
// Constants
// ============================================================================

/// Maximum number of letters allowed in an alphabet.
#define MAX_ALPHABET_SIZE 10

/// Maximum number of letters allowed in a codeword.
#define MAX_CODEWORD_SIZE 6

/// The alphabet used to format codewords as strings.
#define ALPHABET "0123456789"

static_assert(sizeof(ALPHABET) / sizeof(ALPHABET[0]) - 1 == MAX_ALPHABET_SIZE,
              "ALPHABET must contain MAX_ALPHABET_SIZE letters");

/// Helper macro that converts its argument to a string.
#define STRINGIFY(x) _STRINGIFY(x)
#define _STRINGIFY(x) #x

namespace mastermind {

/// Typedef to represent the number of letters in an alphabet.
typedef size_t AlphabetSize;

/// Typedef to represent an index into the alphabet.
typedef size_t Letter;

/// Typedef to represent the number of letters in a codeword.
typedef size_t CodewordSize;

/// Typedef to represent an index into a codeword.
typedef size_t Position;

/// Defines the attributes that secrets and guesses must adhere to.
class CodewordRules // TODO: rename to CodewordAttributes?
{
public:
    /// Creates a rule set with the given parameters.
    ///
    /// Throws `std::invalid_argument` if any of the following requirements
    /// are not satisfied:
    ///
    /// - `1 <= alphabet_size <= MAX_ALPHABET_SIZE
    /// - `1 <= codeword_size <= MAX_CODEWORD_SIZE
    /// - `alphabet_size >= codeword_size` if `is_heterogram` is true
    ///
    constexpr CodewordRules(AlphabetSize alphabet_size,
                            CodewordSize codeword_size,
                            bool is_heterogram)
      : _alphabet_size(alphabet_size),
        _codeword_size(codeword_size),
        _is_heterogram(is_heterogram)
    {
        if (!(alphabet_size >= 1 && alphabet_size <= MAX_ALPHABET_SIZE))
            throw std::invalid_argument("alphabet_size must be between 1 and "
                                        STRINGIFY(MAX_ALPHABET_SIZE));
        if (!(codeword_size >= 1 && codeword_size <= MAX_CODEWORD_SIZE))
            throw std::invalid_argument("codeword_size must be between 1 and "
                                        STRINGIFY(MAX_CODEWORD_SIZE));
        if (is_heterogram && alphabet_size < codeword_size)
            throw std::invalid_argument("alphabet_size must be greater than "
                                        "or equal to codeword_size for "
                                        "heterograms");
    }

    /// Creates a rule set with standard Mastermind rules.
    constexpr CodewordRules() noexcept
      : CodewordRules(AlphabetSize(6), CodewordSize(4), false)
    {
    }

    /// Returns the number of letters in the alphabet.
    constexpr AlphabetSize alphabet_size() const noexcept
    {
        return _alphabet_size;
    }

    /// Returns the number of letters in the codeword.
    constexpr CodewordSize codeword_size() const noexcept
    {
        return _codeword_size;
    }

    /// Returns `true` every letter in the codeword may appear only once.
    constexpr bool is_heterogram() const noexcept { return _is_heterogram; }

    /// Default equality operators by member comparison.
    constexpr bool operator==(const CodewordRules &) const noexcept = default;

    /// Writes codeword attributes to an output stream in the form
    /// "6x4" or "6x4h".
    template <class CharT, class Traits>
    friend std::basic_ostream<CharT, Traits> &
    operator <<(std::basic_ostream<CharT, Traits> &os,
                const CodewordRules &rules)
    {
        os << rules.alphabet_size() << CharT('x') << rules.codeword_size();
        if (rules.is_heterogram())
            os << CharT('h');
        return os;
    }

    /// Reads dictionary definition from an input stream, in the form
    /// "6x4" or "10p4".
    template <class CharT, class Traits>
    friend std::basic_istream<CharT, Traits> &
    operator >>(std::basic_istream<CharT, Traits> &is, CodewordRules &rules)
    {
        size_t n, m;
        CharT sp;
        if (!(is >> n >> sp >> m))
            return is;

        bool is_heterogram;
        if (sp == CharT('x') || sp == CharT('X'))
            is_heterogram = false;
        else if (sp == CharT('p') || sp == CharT('P'))
            is_heterogram = true;
        else
        {
            is.setstate(is.failbit);
            return is;
        }

        try
        {
            rules = CodewordRules(n, m, is_heterogram);
        }
        catch (const std::invalid_argument &)
        {
            is.setstate(is.failbit);
        }
        return is;
    }

private:
    /// Number of letters in the alphabet.
    AlphabetSize _alphabet_size;

    /// Number of letters in the codeword.
    CodewordSize _codeword_size;

    /// `true` if every letter in the codeword may appear only once.
    bool _is_heterogram;
};

/// Helper function used by the Feedback class to generate the mapping
/// from feedback ordinals to feedback outcomes (a, b) for codewords of
/// length M.
template <size_t M>
constexpr std::array<std::pair<size_t, size_t>, (M+1)*(M+2)/2>
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
        (MAX_CODEWORD_SIZE + 1) * (MAX_CODEWORD_SIZE + 2) / 2;

    static_assert(MaxOutcomes - 1 <= std::numeric_limits<ordinal_type>::max(),
                  "ordinal_type is not large enough for MAX_CODEWORD_LENGTH");

    /// Creates the feedback "0A0B".
    constexpr Feedback() noexcept : _ordinal(0) {}

    /// Creates a feedback from its ordinal.
    constexpr explicit Feedback(ordinal_type ordinal) noexcept
      : _ordinal(ordinal) {}

    /// Creates a feedback with the given a and b.
    /// Throws std::invalid_argument if any argument is out of the supported
    /// range.
    constexpr Feedback(size_t a, size_t b)
    {
        if (!(a <= MAX_CODEWORD_SIZE))
            throw std::invalid_argument("a must be between 0 and "
                                        STRINGIFY(MAX_CODEWORD_SIZE));
        if (!(b <= MAX_CODEWORD_SIZE))
            throw std::invalid_argument("b must be between 0 and "
                                        STRINGIFY(MAX_CODEWORD_SIZE));
        size_t ab = a + b;
        if (!(ab <= MAX_CODEWORD_SIZE))
            throw std::invalid_argument("(a + b) must be between 0 and "
                                        STRINGIFY(MAX_CODEWORD_SIZE));

        _ordinal = static_cast<ordinal_type>((ab+1)*ab/2+a);
    }

    /// Returns the ordinal of the feedback.
    constexpr ordinal_type ordinal() const noexcept { return _ordinal; }

    /// Returns the number of matching letters in matching positions.
    constexpr size_t a() const noexcept { return _outcomes[_ordinal].first; }

    /// Returns the number of matching letters in unmatched positions.
    constexpr size_t b() const noexcept { return _outcomes[_ordinal].second; }

    /// Default comparison operators by member comparison.
    constexpr auto operator <=>(const Feedback &other) const noexcept = default;

    /// Returns the feedback corresponding to a perfect match for a
    /// given rule set.
    static constexpr Feedback perfect_match(const CodewordRules &rules) noexcept
    {
        return Feedback(rules.codeword_size(), 0);
    }

    /// Writes a feedback to an output stream in the form "1A2B".
    template <class CharT, class Traits>
    friend std::basic_ostream<CharT, Traits> &
    operator <<(std::basic_ostream<CharT, Traits> &os, const Feedback &feedback)
    {
        return os << feedback.a() << 'A' << feedback.b() << 'B';
    }

    /// Reads a feedback from an input stream in the form "1A2B".
    template <class CharT, class Traits>
    friend std::basic_istream<CharT, Traits> &
    operator >>(std::basic_istream<CharT, Traits> &is, Feedback &feedback)
    {
        size_t a, b;
        CharT A, B;
        if (!(is >> a >> A >> b >> B))
            return is;
        if (!(a <= MAX_CODEWORD_SIZE
              && b <= MAX_CODEWORD_SIZE
              && a + b <= MAX_CODEWORD_SIZE
              && (A == CharT('A') || A == CharT('a'))
              && (B == CharT('B') || B == CharT('b'))))
        {
            is.setstate(is.failbit);
            return is;
        }
        feedback = Feedback(a, b);
        return is;
    }

private:
    /// Ordinal of the feedback.
    ordinal_type _ordinal;

    static constexpr std::array<std::pair<size_t, size_t>, MaxOutcomes>
        _outcomes = generate_outcomes<MAX_CODEWORD_SIZE>();
};

/// Returns a bit pattern of type T
template <class T>
constexpr T cyclic_mask(size_t bits_in_group, size_t num_groups)
{
    T mask(0);
    for (size_t j = 0; j < num_groups; j++)
    {
        mask <<= bits_in_group;
        mask |= T(1);
    }
    return mask;
}

/// Represents a codeword (such as 2587).
class Codeword
{
public:
    /// Creates an empty codeword.
    constexpr Codeword() noexcept : _position_mask(0), _alphabet_mask(0) {}

    /// Creates a codeword with the given letters.
    ///
    /// Throws `std::invalid_argument` if the codeword cannot be represented.
    template <class Iter> // forward iterator
    constexpr Codeword(Iter begin, Iter end)
      : _position_mask(0), _alphabet_mask(0)
    {
        size_t m = 0;
        for (Iter it = begin; it != end; ++it, ++m)
        {
            if (m >= MAX_CODEWORD_SIZE)
                throw std::invalid_argument("codeword too long");

            Letter i = *it;
            if (!(static_cast<size_t>(i) < MAX_ALPHABET_SIZE))
                throw std::invalid_argument("letter out of range");

            _position_mask |= (mask_type(1) << (m * MAX_ALPHABET_SIZE + i));

            const mask_type one = 1;
            const mask_type mask = (one << MAX_CODEWORD_SIZE) - one;
            int k = static_cast<int>(i) * MAX_CODEWORD_SIZE;
            _alphabet_mask |= (_alphabet_mask << 1) & (mask << k);
            _alphabet_mask |= one << k;
        }
    }

    /// Returns the number of letters in the codeword.
    constexpr CodewordSize codeword_size() const noexcept
    {
        return std::popcount(_position_mask);
    }

    /// Returns the minimum size of the alphabet needed to represent
    /// the codeword.  It is equal to the maximum letter plus one,
    /// or zero if the codeword is empty.
    constexpr AlphabetSize alphabet_size() const noexcept
    {
        constexpr int M = MAX_CODEWORD_SIZE;
        static_assert(M > 0, "MAX_CODEWORD_SIZE must be greater than zero");
        return (std::bit_width(_alphabet_mask) + M - 1) / M;
    }

    /// Returns `true` if no letter in the codeword appears more than once.
    constexpr bool is_heterogram() const noexcept
    {
        return !(_alphabet_mask & ~_alphabet_probe);
    }

    /// Returns `true` if this codeword conforms to the given rules.
    constexpr bool conforms_to(const CodewordRules &rules) const noexcept
    {
        return (codeword_size() == rules.codeword_size())
            && (alphabet_size() <= rules.alphabet_size())
            && (!rules.is_heterogram() || is_heterogram());
    }

    /// Default == and != operators by member comparison.
    constexpr bool operator ==(const Codeword &other) const noexcept = default;

    /// Integral type used to store the codeword's bit-mask.
    typedef uint64_t mask_type;

    static_assert(std::is_unsigned_v<mask_type>, "mask_type must be unsigned");

    class LetterIterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = Letter;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using reference = Letter;

        constexpr LetterIterator(Codeword::mask_type position_mask) noexcept
          : _position_mask(position_mask) { }

        constexpr reference operator *() const noexcept
        {
            assert(_position_mask != mask_type(0));
            return std::countr_zero(_position_mask);
        }

        constexpr LetterIterator& operator ++() noexcept
        {
            assert(_position_mask != mask_type(0));
            _position_mask >>= MAX_ALPHABET_SIZE;
            return *this;
        }

        constexpr bool operator ==(const LetterIterator &) const noexcept = default;

    private:
        Codeword::mask_type _position_mask;
    };

    constexpr LetterIterator begin() const noexcept
    {
        return LetterIterator(_position_mask);
    }

    constexpr LetterIterator end() const noexcept
    {
        return LetterIterator(0);
    }

    /// Compares this codeword to another and returns the feedback.
    friend constexpr Feedback
    compare(const Codeword &x, const Codeword &y) noexcept
    {
        int ab = std::popcount(x._alphabet_mask & y._alphabet_mask);
        int a = std::popcount(x._position_mask & y._position_mask);
        int ordinal = (ab+1)*ab/2+a;
        return Feedback(static_cast<Feedback::ordinal_type>(ordinal));
    }

    /// Writes a codeword to an output stream in the form "1224".
    template <class CharT, class Traits>
    friend std::basic_ostream<CharT, Traits> &
    operator <<(std::basic_ostream<CharT, Traits> &os, const Codeword &codeword)
    {
        for (Letter i : codeword)
            os << ALPHABET[i];
        return os;
    }

    /// Reads a codeword from an input stream in the form "1224".
    template <class CharT, class Traits>
    friend std::basic_istream<CharT, Traits> &
    operator >>(std::basic_istream<CharT, Traits> &is, Codeword &codeword)
    {
        std::array<Letter, MAX_CODEWORD_SIZE> letters;

        size_t m = 0;
        for (int cc; (cc = is.peek()) != Traits::eof(); )
        {
            char c = static_cast<char>(cc);
            if (static_cast<int>(c) != cc)
                break;

            size_t i = std::string_view(ALPHABET).find(c);
            if (i == std::string::npos)
                break;

            if (m >= MAX_CODEWORD_SIZE)
            {
                is.setstate(is.failbit);
                return is;
            }

            letters[m++] = static_cast<Letter>(i);
            is.get();
        }

        codeword = Codeword(letters.data(), letters.data() + m);
        return is;
    }

private:
    /// Bit mask where the least significant `M*N` bits are divided
    /// into M groups of N bits, where bit i in group j is set if
    /// the j-th position of the codeword contains the letter i.
    ///
    /// The remaining (most significant) bits are set to zero.
    ///
    /// (`N := MAX_ALPHABET_SIZE`, `M := MAX_CODEWORD_SIZE`.)
    ///
    mask_type _position_mask;

    /// Bit mask where the least significant `M*N` bits are divided
    /// into N groups of M bits, where bit j in group i is set if
    /// letter i appears at least (j+1) times in the codeword.
    ///
    /// The remaining (most significant) bits are set to zero.
    ///
    /// (`N := MAX_ALPHABET_SIZE`, `M := MAX_CODEWORD_SIZE`.)
    ///
    mask_type _alphabet_mask;

    /// Position mask corresponding to codeword `[0] * MAX_CODEWORD_SIZE`.
    static constexpr mask_type _position_probe =
        cyclic_mask<mask_type>(MAX_ALPHABET_SIZE, MAX_CODEWORD_SIZE);

    /// Alphabet mask corresponding to codeword `[0, ..., MAX_ALPHABET_SIZE-1]`.
    static constexpr mask_type _alphabet_probe =
        cyclic_mask<mask_type>(MAX_CODEWORD_SIZE, MAX_ALPHABET_SIZE);
};

struct Constraint
{
    Codeword guess;    // challenge
    Feedback feedback; // response

    constexpr bool operator()(const Codeword &secret) const noexcept
    {
        return compare(secret, guess) == feedback;
    }

    /// Writes a constraint to an output stream in the form "1234:1A1B".
    template <class CharT, class Traits>
    friend std::basic_ostream<CharT, Traits> &
    operator <<(std::basic_ostream<CharT, Traits> &os, const Constraint &c)
    {
        return os << c.guess << CharT(':') << c.feedback;
    }

    /// Reads a constraint from an input stream in the form "1234:1A1B".
    template <class CharT, class Traits>
    friend std::basic_istream<CharT, Traits> &
    operator >>(std::basic_istream<CharT, Traits> &is, Constraint &c)
    {
        Codeword guess;
        Feedback feedback;
        char sp;

        if (!(is >> guess >> sp))
            return is;
        if (sp != CharT(':'))
        {
            is.setstate(is.failbit);
            return is;
        }
        if (!(is >> feedback))
            return is;

        c = Constraint{guess, feedback};
        return is;
    }
};

/// Represents a set of codewords that conform to the given rules and that
/// satisfy the given constraints.
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

    /// Creates a subset of a base set by filtering by a constraint.
    CodewordSet(const CodewordSet &base, const Constraint &constraint);

    /// Returns the number of codewords in the set.
    constexpr size_t size() const noexcept { return _list.size(); }

    /// Returns the codeword at the given index.
    constexpr Codeword operator[](size_t index) const noexcept
    {
        return _list[index];
    }

    /// Returns an iterator to the first codeword in the set.
    std::vector<Codeword>::const_iterator begin() const noexcept
    {
        return _list.begin();
    }

    /// Returns an iterator past the last codeword in the set.
    std::vector<Codeword>::const_iterator end() const noexcept
    {
        return _list.end();
    }

private:
    /// List of codewords.
    std::vector<Codeword> _list;
};

} // namespace mastermind
