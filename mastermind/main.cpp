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

static size_t get_random_index(size_t count)
{
    std::random_device device;
    std::default_random_engine engine(device());
    std::uniform_int_distribution<size_t> distribution(0, count - 1);
    return distribution(engine);
}

static void play(const mastermind::CodewordRules &rules)
{
    using namespace mastermind;

    CodewordPopulation population(rules);
    // TODO: check rules with no admissible codewords

    const Codeword secret = population.get(get_random_index(population.size()));

    for (int round = 1; ; ++round)
    {
        std::cout << "Guess #" << round << "> ";
        std::cout.flush();

        std::string s;
        std::cin >> s;
        Codeword guess = population.get(0);
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

        const Feedback response = compare(guess, secret);
        std::cout << "Score #" << round << ": " << response << std::endl;
        if (response == Feedback::perfect_match(rules))
            break;
    }
}

int main(int argc, const char *argv[])
{
    using namespace mastermind;

    CodewordRules rules;

    std::string alphabet(rules.alphabet());
    PositionSize m = rules.codeword_length();
    CodewordStructure structure = rules.structure();

    //play(rules);
    //return 0;

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
    std::cout << "Random choice: "
        << population.get(get_random_index(population.size())) << std::endl;

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
