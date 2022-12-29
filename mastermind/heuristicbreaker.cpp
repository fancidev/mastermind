
#include "codebreaker.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <vector>

namespace mastermind {

class Heuristic
{
public:
    typedef uint64_t score_type;
    virtual ~Heuristic() = default;

    /// Gets the name of the heuristic.  The returned pointer is valid
    /// during the lifetime of the heuristic object.
    virtual const char *name() const = 0;

    /// Computes a score for the given partitioning of potential secrets.
    /// The number of partitions (including empty ones) MUST cover each
    /// feedback ordinal from 0A0B to mA0B, inclusive, including the
    /// impossible feedback ordinal (m-1)A1B.
    virtual score_type score(std::span<const size_t> partition_sizes) const = 0;

protected:
    Heuristic() = default;
};

class HeuristicCodeBreaker : public CodeBreaker
{
public:
    HeuristicCodeBreaker(const CodewordRules &rules,
                         std::unique_ptr<Heuristic> heuristic)
      : _heuristic(std::move(heuristic)),
        _partition_count(Feedback::perfect_match(rules).ordinal() + 1)
    {
        CodewordPopulation population(rules);
        size_t count = population.size();
        _population.reserve(count);
        for (size_t i = 0; i < count; i++)
            _population.push_back(population.get(i));
        _admissible = std::span(_population);
    }

    virtual const char *name() const override
    {
        return _heuristic->name();
    }

    virtual Codeword make_guess() override
    {
        if (_admissible.empty())
            throw std::runtime_error("no admissible secret");

        // Note: admissible candidates are placed before inadmissible
        // candidates, which automatically favors them if equal score.
        std::span<Codeword> candidate_guesses(_population);

        Codeword chosen_guess;
        Heuristic::score_type chosen_score =
            std::numeric_limits<Heuristic::score_type>::max();

        for (Codeword guess : candidate_guesses)
        {
            std::array<size_t, 256> freq{};
            for (Codeword secret : _admissible)
            {
                Feedback feedback = compare(guess, secret);
                ++freq[feedback.ordinal()];
            }

            Heuristic::score_type score = _heuristic->score(
                std::span(freq).first(_partition_count));
            if (score < chosen_score)
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
    std::unique_ptr<Heuristic> _heuristic;
    size_t _partition_count;
};

namespace heuristics {

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
struct MinimizeAverage : public Heuristic
{
    const char *name() const override { return "minavg"; }

    score_type score(std::span<const size_t> partition_sizes) const override
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
struct MinimizeAverage2 : public MinimizeAverage
{
    const char *name() const override { return "minavg2"; }

    score_type score(std::span<const size_t> partition_sizes) const override
    {
        return MinimizeAverage::score(
            partition_sizes.first(partition_sizes.size() - 2));
    }
};

} // namespace heuristics

std::unique_ptr<CodeBreaker>
create_heuristic_breaker(const CodewordRules &rules, std::string_view name)
{
    std::unique_ptr<Heuristic> heuristic;
    if (name == "minavg")
        heuristic = std::make_unique<heuristics::MinimizeAverage>();
    else if (name == "minavg2")
        heuristic = std::make_unique<heuristics::MinimizeAverage2>();
    else
        throw std::invalid_argument("invalid heuristic name");

    return std::make_unique<HeuristicCodeBreaker>(rules, std::move(heuristic));
}

} // namespace mastermind
