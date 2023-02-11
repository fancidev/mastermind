#include "codemaker.hpp"
#include <cassert>
#include <random>

namespace mastermind {

class StaticCodemaker : public Codemaker
{
public:
    constexpr StaticCodemaker(const Codeword &secret) noexcept
      : _secret(secret) {}

    virtual Feedback respond(const Codeword &guess) override
    {
        return compare(_secret, guess);
    }

private:
    Codeword _secret;
};

std::unique_ptr<Codemaker> create_static_codemaker(const Codeword &secret)
{
    return std::make_unique<StaticCodemaker>(secret);
}

/// Returns a random integer in the range [0, count).
static size_t get_random_index(size_t count)
{
    std::random_device device;
    std::default_random_engine engine(device());
    std::uniform_int_distribution<size_t> distribution(0, count - 1);
    return distribution(engine);
}

/// Creates a standard code maker.
Codeword sample(const CodewordRules &rules)
{
    CodewordPopulation population(rules);
    size_t index = get_random_index(population.size());
    return population.get(index);
}

}
