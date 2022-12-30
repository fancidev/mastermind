#include "codeword.hpp"
#include "codemaker.hpp"
#include "codebreaker.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

static int usage(const char *error)
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

static void play(const mastermind::CodewordRules &rules)
{
    using namespace mastermind;

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

        Codeword guess = code_maker->secret().value(); // TODO: FIXME
        try
        {
            guess = Codeword::from_string(s, rules);
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

static void self_play(const mastermind::CodewordRules &rules)
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

static void test_breaker(const mastermind::CodewordRules &rules,
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

int main(int argc, const char *argv[])
{
    using namespace mastermind;

    CodewordRules rules;

    std::string alphabet(rules.alphabet());
    PositionSize m = rules.codeword_length();
    CodewordStructure structure = rules.structure();

//    self_play(rules);
    const char *heuristics[] = {
        "maxparts", "maxparts~", "minavg", "minavg~", "minmax", "minmax~", "entropy", "entropy~",
    };
    for (const char *heuristic : heuristics)
    {
        test_breaker(rules, heuristic);
    }
    return 0;

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
                structure = CodewordStructure::Heterogram;
                break;
            case 'h':
                return usage(nullptr);
            case 'm':
                m = std::atoi(val);
                break;
            case 'n':
            {
                int n = std::atoi(val);
                if (!(n >= 0 && n <= MAX_ALPHABET_SIZE))
                    return usage("alphabet size out of range");
                static_assert(MAX_ALPHABET_SIZE <= 10, "not supported");
                alphabet = std::string("1234567890").substr(0, n);
                break;
            }
            default:
                return usage("unknown option");
        }
    }

//    try
//    {
//        rules = mastermind::CodewordRules(alphabet, m, structure);
//    }
//    catch (const std::invalid_argument &ex)
//    {
//        return usage(ex.what());
//    }

    std::cout << "Codeword rules:" << std::endl;
    std::cout << "  Alphabet size: " << rules.alphabet_size() << std::endl;
    std::cout << "  Codeword length: " << rules.codeword_length() << std::endl;
    std::cout << "  Structure: " << to_string(rules.structure()) <<  std::endl;
    std::cout << "Perfect match is " <<
        Feedback::perfect_match(rules) << std::endl;

    CodewordPopulation population(rules);
    std::cout << "Population size: " << population.size() << std::endl;

    std::cout << "First 5:";
    for (size_t index = 0; index < 5 && index < population.size(); index++)
    {
        std::cout << " " << population.get(index);
    }
    std::cout << std::endl;

    std::cout << "Last  5:";
    for (size_t index = 0; index < 5; index++)
    {
        //if (index + )
        std::cout << " " << population.get(population.size() - (5 - index));
    }
    std::cout << std::endl;

//    Codeword guess = Codeword::from_string("1357", rules);
//    Codeword secret = Codeword::from_string("2337", rules);
//    std::cout << "compare(" << guess.to_string("ABCDEFGHIJ") << ", "
//        << secret.to_string("ABCDEFGHIJ") << ") = "
//        << compare(guess, secret).to_string() << std::endl;

    return 0;
}
