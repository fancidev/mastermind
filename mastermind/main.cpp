#include "codeword.hpp"
#include "codemaker.hpp"
#include "codebreaker.hpp"
#include "canonical.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace mastermind;

int usage(const char *error)
{
    if (error)
        std::cerr << "error: " << error << std::endl;
    std::cerr <<
        "usage: mastermind [options]\n"
        "Options:\n"
        "    -d      require any letter in codeword to appear only once\n"
        "    -h      display this help text and exit\n"
        "    -m M    number of letters in codeword (default: 4)\n"
        "    -n N    number of letters in alphabet (default: 6)\n"
        ;
    return error ? 1 : 0;
}

void play(const mastermind::CodewordRules &rules)
{
    std::unique_ptr<CodeMaker> code_maker(create_code_maker(rules));

    for (int round = 1; ; ++round)
    {
        std::cout << "Guess #" << round << "> ";
        std::cout.flush();

        std::string s;
        std::cin >> s;
        if (s == "?")
        {
            std::optional<Codeword> secret = code_maker->secret();
            if (secret.has_value())
                std::cout << "Secret: " << secret.value() << std::endl;
            else
                std::cout << "Secret is not available." << std::endl;
            --round;
            continue;
        }
        if (s == "quit")
            break;

        Codeword guess(rules);
        try
        {
            if (!(std::istringstream(s) >> guess))
                throw std::invalid_argument("cannot read guess");
        }
        catch (const std::invalid_argument &ex)
        {
            std::cerr << "error: " << ex.what() << std::endl;
            --round;
            continue;
        }

        const Feedback response = code_maker->respond(guess);
        std::cout << "Score #" << round << ": " << response << std::endl;
        if (response == Feedback::perfect_match(rules))
            break;
    }
}

void self_play(const mastermind::CodewordRules &rules)
{
    using namespace mastermind;

    std::unique_ptr<CodeMaker> maker(create_code_maker(rules));
    //std::unique_ptr<CodeBreaker> breaker(create_simple_code_breaker(rules));
    const char *name = "minavg";
    std::unique_ptr<CodeBreaker> breaker(create_heuristic_breaker(rules, name));

    for (int round = 1; round < 10000; ++round)
    {
        Codeword guess = breaker->make_guess();
        Feedback response = maker->respond(guess);

        std::cout << "#" << round << "  ";
        std::cout << guess << "  " << response << std::endl;
        if (response == Feedback::perfect_match(rules))
            break;
        breaker->step(guess, response);
    }
}

void test_breaker(const mastermind::CodewordRules &rules,
                  const char *name)
{
    using namespace mastermind;

    CodewordPopulation population(rules);
    size_t total_steps = 0;
    size_t worst_steps = 0;
    for (size_t index = 0; index < population.size(); ++index)
    {
        std::unique_ptr<CodeMaker> maker(create_code_maker(rules, population.get(index)));
        //std::unique_ptr<CodeBreaker> breaker(create_simple_code_breaker(rules));
        std::unique_ptr<CodeBreaker> breaker(create_heuristic_breaker(rules, name));

        size_t round = 0;
        while (++round < 10000)
        {
            Codeword guess = breaker->make_guess();
            Feedback response = maker->respond(guess);
            if (response == Feedback::perfect_match(rules))
                break;
            breaker->step(guess, response);
        }
        total_steps += round;
        worst_steps = std::max(worst_steps, round);
        if (index > 10000)
            break;
    }
    std::cout << "**" << name << std::endl;
    std::cout << "  Avg steps: " << total_steps << "/" << population.size()
        << " = " << static_cast<double>(total_steps) / population.size()
        << std::endl;
    std::cout << "  Max steps: " << worst_steps << std::endl;
}

static void _display_canonical_guesses(const CanonicalCodewordSequence &group,
                                       std::span<const Codeword> candidate_guesses,
                                       size_t &counter)
{
    if (group.sequence().size() >= 2)
        return;

    std::vector<CanonicalCodewordSequence> guesses =
        get_canonical_guesses(group, candidate_guesses);
    for (const CanonicalCodewordSequence &g : guesses)
    {
        ++counter;
        std::cout << counter << ">";
        for (auto cw : g.sequence())
            std::cout << " " << cw;
        std::cout << std::endl;

//        if (g.is_singleton())
//        {
//            for (const Codeword &guess : g.sequence())
//                std::cout << guess << " ";
//            std::cout << "-> all" << std::endl;
//        }
//        else
//        {
//            std::cout << g.sequence().back() << std::endl;
//            _display_canonical_guesses(g, candidate_guesses, counter);
//        }
        _display_canonical_guesses(g, candidate_guesses, counter);
    }
}

void display_canonical_guesses(const CodewordRules &rules)
{
    // TODO: add Codeword::enumerate()
    CodewordPopulation population(rules);
    std::vector<Codeword> candidate_guesses(population.begin(), population.end());
    CanonicalCodewordSequence sequence(rules);

    size_t counter = 0;
    _display_canonical_guesses(sequence, candidate_guesses, counter);
    std::cout << "*** Counter = " << counter << std::endl;
}

int main(int argc, const char *argv[])
{
    using namespace mastermind;

    CodewordRules rules;

#if 1
    display_canonical_guesses(rules);
    return 0;
#endif

#if 0
    AlphabetSize n = rules.alphabet_size();
    PositionSize m = rules.codeword_length();
    bool heterogram = rules.heterogram();

//    self_play(rules);
#if 1
    const char *heuristics[] = {
        "maxparts", // "maxparts~", "minavg", "minavg~", "minmax", "minmax~", "entropy", "entropy~",
    };
    for (const char *heuristic : heuristics)
    {
        test_breaker(rules, heuristic);
    }
    return 0;
#endif

    const std::string options_requring_argument = "mn";

    // Parse command line arguments.
    for (int i = 1; i < argc; i++)
    {
        const char *opt = argv[i];
        if (opt[0] != '-')
            return usage("missing option");
        if (opt[1] == '\0')
            return usage("missing option");

        const char *val = nullptr;
        if (options_requring_argument.find(opt[1]) != std::string::npos)
        {
            if (opt[2])
                val = &opt[2];
            else if (i + 1 == argc)
                return usage("option requires an argument");
            else
                val = argv[++i];
        }

        switch (opt[1])
        {
            case 'd':
                heterogram = true;
                break;
            case 'h':
                return usage(nullptr);
            case 'm':
                m = std::atoi(val);
                break;
            case 'n':
                n = std::atoi(val);
                break;
            default:
                return usage("unknown option");
        }
    }

    try
    {
        rules = mastermind::CodewordRules(n, m, heterogram);
    }
    catch (const std::invalid_argument &ex)
    {
        return usage(ex.what());
    }

    std::cout << "Codeword rules: " << rules << std::endl;
    std::cout << "  Alphabet size: " << rules.alphabet_size() << std::endl;
    std::cout << "  Codeword length: " << rules.codeword_length() << std::endl;
    std::cout << "  Requires hetero: " << rules.heterogram() <<  std::endl;
    std::cout << "Perfect match is " <<
        Feedback::perfect_match(rules) << std::endl;

    CodewordPopulation population(rules);
    std::cout << "Population size: " << population.size() << std::endl;

    std::vector<Codeword> all(population.begin(), population.end());

    std::cout << "First 5:";
    for (size_t index = 0; index < 5 && index < all.size(); index++)
    {
        std::cout << " " << all[index];
    }
    std::cout << std::endl;

    std::cout << "Last  5:";
    for (size_t index = std::max<size_t>(5, all.size()) - 5; index < all.size(); index++)
    {
        std::cout << " " << all[index];
    }
    std::cout << std::endl;

//    Codeword guess = Codeword::from_string("1357", rules);
//    Codeword secret = Codeword::from_string("2337", rules);
//    std::cout << "compare(" << guess.to_string("ABCDEFGHIJ") << ", "
//        << secret.to_string("ABCDEFGHIJ") << ") = "
//        << compare(guess, secret).to_string() << std::endl;

    return 0;
#endif
}
