#include "codemaker.hpp"
#include "tracker.hpp"

#include <cassert>
#include <random>
#include <vector>

namespace mastermind {

class StaticCodemaker : public Codemaker
{
public:
    constexpr StaticCodemaker(const Codeword &secret) noexcept
      : _secret(secret) {}

    virtual Feedback respond(const Codeword &guess) override
    {
        return Correlation(_secret, guess);
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
    CodewordSet population(rules);
    return population.possible_secrets()[get_random_index(population.possible_secrets().size())];
}

//class DynamicCodemaker : public Codemaker
//{
//public:
//    DynamicCodemaker(const CodewordRules &rules)
//      : _rules(rules)
//    {
//        CodewordPopulation p(rules);
//        _secrets.assign(p.begin(), p.end());
//    }
//
//    virtual Feedback respond(const Codeword &guess) override
//    {
//        assert(guess.conforms_to(_rules));
//
//        // Partition
//        std::vector<Feedback> feedbacks;
//        feedbacks.reserve(_secrets.size());
//        for (Codeword secret : _secrets)
//            feedbacks.push_back(compare(secret, guess));
//
//        // partition it
//
//
//        return compare(_secret, guess);
//    }
//
//private:
//    CodewordRules _rules;
//    std::vector<Codeword> _secrets;
//};

}
