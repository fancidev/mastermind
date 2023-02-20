#include "codebreaker.hpp"
#include "tracker.hpp"

#include <algorithm>
#include <vector>

namespace mastermind {

class SimpleCodeBreaker : public CodeBreaker
{
public:
    SimpleCodeBreaker(const CodewordRules &rules) : _tracker(rules)
    {
    }

    virtual const char *name() const override
    {
        return "simple";
    }

    virtual Codeword make_guess() override
    {
        if (_tracker.possible_secrets().empty())
            throw std::runtime_error("no possible secret");
        return _tracker.possible_secrets().front();
    }

    virtual void step(const Constraint &constraint) override
    {
        _tracker.push_constraint(constraint);
    }

private:
    CodewordSet _tracker;
};

std::unique_ptr<CodeBreaker>
create_simple_code_breaker(const CodewordRules &rules)
{
    return std::make_unique<SimpleCodeBreaker>(rules);
}

} // namespace mastermind
