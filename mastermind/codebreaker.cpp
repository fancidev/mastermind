#include "codebreaker.hpp"

#include <algorithm>
#include <vector>

namespace mastermind {

class SimpleCodeBreaker : public CodeBreaker
{
public:
    SimpleCodeBreaker(const CodewordRules &rules) : _possible_secrets(rules)
    {
    }

    virtual const char *name() const override
    {
        return "simple";
    }

    virtual Codeword make_guess() override
    {
        if (_possible_secrets.empty())
            throw std::runtime_error("no possible secret");
        return _possible_secrets[0];
    }

    virtual void step(const Constraint &constraint) override
    {
        _possible_secrets.push_constraint(constraint);
    }

private:
    CodewordSet _possible_secrets;
};

std::unique_ptr<CodeBreaker>
create_simple_code_breaker(const CodewordRules &rules)
{
    return std::make_unique<SimpleCodeBreaker>(rules);
}

} // namespace mastermind
