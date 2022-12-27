#include "codeword.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

static int usage(const char *error)
{
    if (error)
        std::cerr << "error: " << error << std::endl;
    std::cerr
        << "usage: mastermind [options]\n"
        << "Options:\n"
        << "    -n N   Set number of letters in the alphabet (default: 6)\n"
        << "    -m M   Set number of letters in the secret (default: 4)\n"
        << "    -d     Require any letter in the secret to appear only once\n"
        << "    -h     display usage and exit\n"
        ;
    return error ? 1 : 0;
}

static size_t get_random_index(size_t count)
{
    std::random_device device;
    std::default_random_engine engine(device());
    std::uniform_int_distribution<size_t> distribution(0, count - 1);
    return distribution(engine);
}

int main(int argc, const char *argv[])
{
    using namespace mastermind;

    CodewordRules rules;
    AlphabetSize n = rules.alphabet_size();
    PositionSize m = rules.codeword_length();
    CodewordStructure structure = rules.structure();

    const std::string options_requring_argument = "nm";

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
                n = std::atoi(val);
                break;
            default:
                return usage("unknown option");
        }
    }

    try
    {
        rules = mastermind::CodewordRules(n, m, structure);
    }
    catch (const std::invalid_argument &ex)
    {
        return usage(ex.what());
    }

    std::cout << "Codeword rules:" << std::endl;
    std::cout << "  Alphabet size: " << rules.alphabet_size() << std::endl;
    std::cout << "  Codeword length: " << rules.codeword_length() << std::endl;
    std::cout << "  Structure: " << to_string(rules.structure()) <<  std::endl;
    std::cout << "Perfect match is " <<
        Feedback::perfect_match(rules).to_string() << std::endl;

    CodewordPopulation population(rules);
    std::string alpha("ABCDEFGHIJKLMN");
    alpha = alpha.substr(0, rules.alphabet_size());
    std::cout << "Population size: " << population.size() << std::endl;
    std::cout << "Random choice: "
        << population.get(get_random_index(population.size())).to_string(alpha)
        << std::endl;

    std::cout << "First 5:";
    for (size_t index = 0; index < 5 && index < population.size(); index++)
    {
        std::cout << " " << population.get(index).to_string(alpha);
    }
    std::cout << std::endl;

    std::cout << "Last  5:";
    for (size_t index = 0; index < 5; index++)
    {
        //if (index + )
        std::cout << " " << population.get(population.size() - (5 - index)).to_string(alpha);
    }
    std::cout << std::endl;

    const char *alphabet = "1234567890";
    Codeword guess = Codeword::from_string("1357", alphabet);
    Codeword secret = Codeword::from_string("2337", alphabet);
    std::cout << "compare(" << guess.to_string("ABCDEFGHIJ") << ", "
        << secret.to_string("ABCDEFGHIJ") << ") = "
        << compare(guess, secret).to_string() << std::endl;

    return 0;
}
