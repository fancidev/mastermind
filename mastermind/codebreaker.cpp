#include "codebreaker.hpp"

#include <algorithm>
#include <vector>

namespace mastermind {

class SimpleCodeBreaker : public CodeBreaker
{
public:
    SimpleCodeBreaker(const CodewordRules &rules)
    {
        CodewordPopulation population(rules);
        size_t count = population.size();
        _potential_secrets.reserve(count);
        for (size_t i = 0; i < count; i++)
            _potential_secrets.emplace_back(population.get(i));
    }

    virtual const char *name() const override
    {
        return "simple";
    }

    virtual Codeword make_guess() override
    {
        if (_potential_secrets.empty())
            throw std::runtime_error("no potential secret");
        return _potential_secrets[0];
    }

    virtual void step(const Codeword &guess, Feedback response) override
    {
        auto filter = [=](const Codeword &secret)
        {
            return compare(guess, secret) == response;
        };
#if 0
        auto it = std::partition(_potential_secrets.begin(),
                                 _potential_secrets.end(),
                                 filter);
#else
        auto it = std::stable_partition(_potential_secrets.begin(),
                                        _potential_secrets.end(),
                                        filter);
#endif
        _potential_secrets.resize(it - _potential_secrets.begin());
    }

private:
    std::vector<Codeword> _potential_secrets;
};

std::unique_ptr<CodeBreaker>
create_simple_code_breaker(const CodewordRules &rules)
{
    return std::make_unique<SimpleCodeBreaker>(rules);
}

} // namespace mastermind
