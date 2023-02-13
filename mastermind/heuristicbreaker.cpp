
#include "codebreaker.hpp"

#include <algorithm>
#include <array>
#include <cmath>
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
concept Heuristic = requires(std::span<const size_t> partition_sizes)
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
        CodewordSet population(rules);
        _population.assign(population.begin(), population.end());
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
/// If `AdjustPerfectPartition` is true, the size of the perfect partition
/// is always treated as zero, essentially favoring a non-empty perfect
/// partition.  The rationale is that the perfect partition requires no
/// more guesses.
template <bool AdjustPerfectPartition = false>
struct MinimizeWorstCase
{
    using score_type = std::array<size_t, Feedback::MaxOutcomes>;

    static constexpr const char *name() noexcept
    {
        return AdjustPerfectPartition ? "minmax~" : "minmax";
    }

    /// Returns a copy of the partition sizes sorted in descending order.
    static constexpr score_type evaluate(
        std::span<const size_t> partition_sizes) noexcept
    {
        score_type score{};
        std::copy(partition_sizes.begin(),
                  partition_sizes.end(),
                  score.begin());

        if constexpr (AdjustPerfectPartition)
            score[partition_sizes.size() - 1] = 0;

        std::sort(score.begin(),
                  score.begin() + partition_sizes.size(),
                  std::greater());
        return score;
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
/// If `AdjustPerfectPartition` is true, the size of the perfect partition
/// is excluded from the summation above, essentially favoring guesses that
/// are admissible.  The rationale is that the perfect partition requires
/// no more guesses.
template <bool AdjustPerfectPartition = false>
struct MinimizeAverage
{
    using score_type = uint64_t;

    static constexpr const char *name() noexcept
    {
        return AdjustPerfectPartition ? "minavg~" : "minavg";
    }

    static constexpr score_type evaluate(
        std::span<const size_t> partition_sizes) noexcept
    {
        auto op = [](score_type score, size_t count) -> score_type {
            return score + count * count;
        };

        if constexpr (AdjustPerfectPartition)
            partition_sizes = partition_sizes.first(partition_sizes.size() - 2);

        return std::accumulate(partition_sizes.begin(),
                               partition_sizes.end(),
                               score_type(0),
                               op);
    }
};

/// Scores a guess by (a proxy for) the expected number of further guesses
/// needed to reveal the secret (Neuwirth, 1982).
///
/// Specifically, given P partitions of size N[1], ..., N[P] that sum to N,
/// the objective is to maximize the posterior entropy, defined as
///
///   `-[ (N[1]/N)*log(N[1]/N) + ... + (N[P]/N)*log(N[P]/N) ]`
///
/// which is equivalent to minimizing
///
///   `N[1]*log(N[1]) + ... + N[P]*log(N[P])`
///
/// Assuming `log(N[i])` to be proportional to the number of further guesses
/// needed to reveal a secret in the i-th partition, and assumining secrets
/// to be uniformly distributed, the objective function is then proportional
/// to the expected number of further guesses needed.
///
/// If `AdjustPerfectPartition` is true, the perfect partition is excluded
/// from the summation above.  The rationale is that the perfect partition
/// requires no more guesses.
///
template <bool AdjustPerfectPartition = false>
struct MaximizeEntropy
{
    using score_type = double;
//    /// Type of the score (wrapped double to avoid numerical instability
//    /// during comparison).
//    typedef util::wrapped_float<double, 100> score_t;

    static constexpr const char *name() noexcept
    {
        return AdjustPerfectPartition ? "entropy~" : "entropy";
    }

    static constexpr score_type evaluate(
        std::span<const size_t> partition_sizes) noexcept
    {
        auto op = [](score_type score, size_t count) -> score_type
        {
            if (count >= 2)
                return score + std::log(static_cast<double>(count)) * count;
            else
                return score;
        };

        if constexpr (AdjustPerfectPartition)
            partition_sizes = partition_sizes.first(partition_sizes.size() - 2);

        // TODO: check the below
//        if (apply_correction && freq[freq.size()-1]) // 4A0B
//        {
//            s -= 2.0 * std::log(2.0);
//        }

        return std::accumulate(partition_sizes.begin(),
                               partition_sizes.end(),
                               0.0,
                               op);
    }
};

/// Scores a guess by the number of (non-empty) partitions it produces.
/// The rationale is that more partitions means more information which
/// leads to fewer further steps.
///
/// If `AdjustPerfectPartition` is true, the perfect partition counts
/// as 1.5 partitions, making it slightly more favorable.  This makes
/// limited sense and is provided for completeness only.
template <bool AdjustPerfectPartition = false>
struct MaximizePartitions
{
    using score_type = int;

    static constexpr const char *name() noexcept
    {
        return AdjustPerfectPartition ? "maxparts~" : "maxparts";
    }

    /// Returns the negative number of non-empty partitions (to minimize).
    static constexpr int evaluate(
        std::span<const size_t> partition_sizes) noexcept
    {
        size_t count = std::count_if(partition_sizes.begin(),
                                     partition_sizes.end(),
                                     [](size_t k) { return k > 0; });

        if constexpr (AdjustPerfectPartition)
            count = 2 * count + partition_sizes.back();

        return -static_cast<int>(count); // negate for minimization
    }
};

} // namespace heuristics

std::unique_ptr<CodeBreaker>
create_heuristic_breaker(const CodewordRules &rules, std::string_view name)
{
    if (name == "minavg")
    {
        return std::make_unique<HeuristicCodeBreaker<
            heuristics::MinimizeAverage<false>>>(rules);
    }
    else if (name == "minavg~")
    {
        return std::make_unique<HeuristicCodeBreaker<
            heuristics::MinimizeAverage<true>>>(rules);
    }
    else if (name == "minmax")
    {
        return std::make_unique<HeuristicCodeBreaker<
            heuristics::MinimizeWorstCase<false>>>(rules);
    }
    else if (name == "minmax~")
    {
        return std::make_unique<HeuristicCodeBreaker<
            heuristics::MinimizeWorstCase<true>>>(rules);
    }
    else if (name == "entropy")
    {
        return std::make_unique<HeuristicCodeBreaker<
            heuristics::MaximizeEntropy<false>>>(rules);
    }
    else if (name == "entropy~")
    {
        return std::make_unique<HeuristicCodeBreaker<
            heuristics::MaximizeEntropy<true>>>(rules);
    }
    else if (name == "maxparts")
    {
        return std::make_unique<HeuristicCodeBreaker<
            heuristics::MaximizePartitions<false>>>(rules);
    }
    else if (name == "maxparts~")
    {
        return std::make_unique<HeuristicCodeBreaker<
            heuristics::MaximizePartitions<true>>>(rules);
    }
    else
    {
        throw std::invalid_argument("invalid heuristic name");
    }
}

} // namespace mastermind
