#include "codebreaker.hpp"

#include <algorithm>
#include <vector>

namespace mastermind {

class SimpleCodeBreaker : public CodeBreaker
{
public:
    SimpleCodeBreaker(const CodewordRules &rules)
    {
        CodewordSet population(rules);
        _potential_secrets.assign(population.begin(), population.end());
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
#if 0
        auto it = std::partition(_potential_secrets.begin(),
                                 _potential_secrets.end(),
                                 Constraint{guess, response});
#else
        auto it = std::stable_partition(_potential_secrets.begin(),
                                        _potential_secrets.end(),
                                        Constraint{guess, response});
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
