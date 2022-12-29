#include "codemaker.hpp"
#include <cassert>
#include <random>

namespace mastermind {

class StandardCodeMaker : public CodeMaker
{
public:
    constexpr StandardCodeMaker(const CodewordRules &rules,
                                const Codeword &secret) noexcept
      : _rules(rules), _secret(secret)
    {
        assert(secret.conforms_to(rules));
    }

    virtual Feedback respond(const Codeword &guess) override
    {
        return compare(_secret, guess);
    }

    virtual std::optional<Codeword> secret() const override
    {
        return _secret;
    }

private:
    CodewordRules _rules;
    Codeword _secret;
};

/// Returns a random integer in the range [0, count).
static size_t get_random_index(size_t count)
{
    std::random_device device;
    std::default_random_engine engine(device());
    std::uniform_int_distribution<size_t> distribution(0, count - 1);
    return distribution(engine);
}

/// Creates a standard code maker.
std::unique_ptr<CodeMaker> create_code_maker(
    const CodewordRules &rules, std::optional<Codeword> secret)
{
    if (secret.has_value())
    {
        if (!secret.value().conforms_to(rules))
            throw std::invalid_argument("secret does not conform to rules");
    }
    else
    {
        CodewordPopulation population(rules);
        size_t index = get_random_index(population.size());
        secret = population.get(index);
    }
    return std::make_unique<StandardCodeMaker>(rules, secret.value());
}

}
