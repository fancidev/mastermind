
#include "codebreaker.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <vector>

namespace mastermind {

template <class H>
concept Heuristic = requires(std::span<size_t> partition_sizes)
{
    /// `H::score_type` shall be the return type of `H::evaluate()`.
    typename H::score_type;

    /// `H::name()` shall return the name of the heuristic.  The returned
    /// pointer must have static lifetime.
    { H::name() } -> std::convertible_to<const char *>;

    /// `H::evaluate()` shall return a value of `H::score_type` for the
    /// given partitioning of potential secrets.
    ///
    /// `H::evaluate()` is guaranteed to be called with the size of the
    /// partitions corresponding to feedback ordinals 0A0B through mA0B,
    /// including empty ones as well as the impossible feedback (m-1)A1B.
    { H::evaluate(partition_sizes) } ->
        std::convertible_to<typename H::score_type>;
};

template <Heuristic H>
class HeuristicCodeBreaker : public CodeBreaker
{
public:
    HeuristicCodeBreaker(const CodewordRules &rules)
      : _partition_count(Feedback::perfect_match(rules).ordinal() + 1)
    {
        CodewordPopulation population(rules);
        size_t count = population.size();
        _population.reserve(count);
        for (size_t i = 0; i < count; i++)
            _population.push_back(population.get(i));
        _admissible = std::span(_population);
    }

    virtual const char *name() const override { return H::name(); }

    virtual Codeword make_guess() override
    {
        if (_admissible.empty())
            throw std::runtime_error("no admissible secret");

        // Note: admissible candidates are placed before inadmissible
        // candidates, which automatically favors them if equal score.
        std::span<Codeword> candidate_guesses(_population);

        Codeword chosen_guess;
        typename H::score_type chosen_score{};

        for (size_t index = 0; index < candidate_guesses.size(); index++)
        {
            Codeword guess = candidate_guesses[index];

            std::array<size_t, Feedback::MaxOutcomes> freq{};
            for (Codeword secret : _admissible)
            {
                Feedback feedback = compare(guess, secret);
                ++freq[feedback.ordinal()];
            }

            typename H::score_type score = H::evaluate(
                std::span(freq).first(_partition_count));
            if (index == 0 || score < chosen_score)
            {
                chosen_guess = guess;
                chosen_score = score;
            }
        }
        return chosen_guess;
    }

    virtual void step(const Codeword &guess, Feedback response) override
    {
        auto filter = [=](const Codeword &secret)
        {
            return compare(guess, secret) == response;
        };
#if 0
        auto it = std::partition(_admissible.begin(),
                                 _admissible.end(),
                                 filter);
#else
        auto it = std::stable_partition(_admissible.begin(),
                                        _admissible.end(),
                                        filter);
#endif
        _admissible = _admissible.first(it - _admissible.begin());
    }

private:
    std::vector<Codeword> _population;
    std::span<Codeword> _admissible;
    size_t _partition_count;
};

namespace heuristics {

/// Scores a guess by the worst-case number of posterior potential secrets
/// (Knuth, 1976).  If two guesses produce the same worst-case number of
/// posterior potential secrets, the second-to-worst number is compared,
/// and so on.
///
/// ??? If a correction is applied (which is the default), a cell corresponding
/// ??? to a perfect match is assumed to be empty, because it will not "remain"
/// ??? after this guess.
struct MinimizeWorstCase
{
    using score_type = std::array<size_t, Feedback::MaxOutcomes>;

    static constexpr const char *name() noexcept { return "minmax"; }

    /// Returns a copy of the partition sizes sorted in descending order.
    static constexpr score_type evaluate(std::span<size_t> partition_sizes) noexcept
    {
        score_type score{};
        std::copy(partition_sizes.begin(),
                  partition_sizes.end(),
                  score.begin());
        std::sort(score.begin(), score.end(), std::greater());
        return score;
    }
};

/// Similar to MinimizeWorstCase, except that the size of a perfect partition
/// (if any) is treated as being zero.  The rationale is that the perfect
/// partition requires no more guess.
struct MinimizeWorstCase2 : MinimizeWorstCase
{
    static constexpr const char *name() noexcept { return "minmax2"; }

    static constexpr score_type evaluate(std::span<size_t> partition_sizes) noexcept
    {
        partition_sizes[partition_sizes.size() - 1] = 0;
        return MinimizeWorstCase::evaluate(partition_sizes);
    }
};

/// Scores a guess by the posterior expected number of potential secrets,
/// assuming equal prior probability of each potential secret (Irving, 1979).
///
/// Given P partitions of size N[1], ..., N[P] that sum to N, the objective
/// function to minimize is
///
///   `N[1]*(N[1]/N) + ... + N[P]*(N[P]/N)`
///
/// which is equivalent to minimizing the score
///
///   `N[1]**2 + ... + N[P]**2`
///
struct MinimizeAverage
{
    using score_type = uint64_t;

    static constexpr const char *name() noexcept { return "minavg"; }

    static constexpr score_type evaluate(std::span<size_t> partition_sizes) noexcept
    {
        auto op = [](score_type score, size_t count) -> score_type {
            return score + count * count;
        };
        return std::accumulate(partition_sizes.begin(),
                               partition_sizes.end(),
                               score_type(0),
                               op);
    }
};

/// Similar to MinimizeAverage, except that the partition of perfect match
/// (if any) is excluded from the score.  The rationale is that the perfect
/// partition requires no more guess.
struct MinimizeAverage2 : MinimizeAverage
{
    static constexpr const char *name() noexcept { return "minavg2"; }

    static constexpr score_type evaluate(std::span<size_t> partition_sizes) noexcept
    {
        return MinimizeAverage::evaluate(
            partition_sizes.first(partition_sizes.size() - 2));
    }
};

} // namespace heuristics

std::unique_ptr<CodeBreaker>
create_heuristic_breaker(const CodewordRules &rules, std::string_view name)
{
    if (name == "minavg")
        return std::make_unique<HeuristicCodeBreaker<heuristics::MinimizeAverage>>(rules);
    else if (name == "minavg2")
        return std::make_unique<HeuristicCodeBreaker<heuristics::MinimizeAverage2>>(rules);
    else if (name == "minmax")
        return std::make_unique<HeuristicCodeBreaker<heuristics::MinimizeWorstCase>>(rules);
    else
        throw std::invalid_argument("invalid heuristic name");
}

} // namespace mastermind
