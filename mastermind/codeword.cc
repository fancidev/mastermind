#include "codeword.hpp"

#include <algorithm>
#include <iterator>
#include <random>

namespace mastermind {

/// Returns a vector v with (m + 1) elements (m := codeword length), where
/// v[j] := the size of the sub-population where the first j letters in the
/// codeword are fixed.
template <size_t M>
static constexpr std::array<size_t, M + 1>
get_sub_population_sizes(const CodewordRules &rules) noexcept
{
    std::array<size_t, M + 1> sizes {};
    size_t count = 1;
    const size_t n = rules.alphabet_size();
    const size_t m = rules.codeword_size();
    for (Position j = m; j > 0; --j)
    {
        sizes[j] = count;
        if (rules.is_heterogram())
            count *= (n - j + 1);
        else
            count *= n;
    }
    sizes[0] = count;
    return sizes;
}

template <size_t M>
static constexpr Codeword
get_codeword_at(size_t index,
                const CodewordRules &rules,
                const std::array<size_t, M + 1> &sub_population_sizes) noexcept
{
    const size_t n = rules.alphabet_size();
    const size_t m = rules.codeword_size();

    std::array<Letter, MAX_CODEWORD_SIZE> letters {};
    if (rules.is_heterogram())
    {
        std::array<Letter, MAX_ALPHABET_SIZE> alphabet {};
        std::iota(alphabet.begin(), alphabet.end(), Letter(0));
        for (Position j = 0; j < m; j++)
        {
            size_t i = index / sub_population_sizes[j + 1];
            letters[j] = alphabet[i];
            std::copy(alphabet.cbegin() + i + 1,
                      alphabet.cbegin() + n - j,
                      alphabet.begin() + i);
            index %= sub_population_sizes[j + 1];
        }
    }
    else
    {
        for (Position j = 0; j < m; j++)
        {
            size_t i = index / sub_population_sizes[j + 1];
            letters[j] = static_cast<Letter>(i);
            index %= sub_population_sizes[j + 1];
        }
    }
    return Codeword(letters.data(), letters.data() + m);
}

CodewordSet::CodewordSet(const CodewordRules &rules)
{
    std::array<size_t, MAX_CODEWORD_SIZE + 1> sub_population_sizes(
        get_sub_population_sizes<MAX_CODEWORD_SIZE>(rules));

    size_t n = sub_population_sizes[0];
    _list.reserve(n);
    for (size_t index = 0; index < n; ++index)
        _list.push_back(get_codeword_at<MAX_CODEWORD_SIZE>(index, rules, sub_population_sizes));
}

CodewordSet::CodewordSet(const CodewordSet &base, const Constraint &constraint)
{
    std::copy_if(base.begin(),
                 base.end(),
                 std::back_inserter(_list),
                 constraint);
}

} // namespace mastermind
