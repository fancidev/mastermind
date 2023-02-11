#include "codeword.hpp"
#include "codemaker.hpp"
#include "codebreaker.hpp"
#include "canonical.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>

using namespace mastermind;

int usage(const char *error)
{
    if (error)
        std::cerr << "error: " << error << std::endl;
    std::cerr <<
        "usage: mastermind [options] action\n"
        "Options:\n"
        "    -b breaker  Specify the codebreaker; must be one of:\n"
        "                simple,minavg(default),minmax,entropy,maxparts\n"
        "    -c guess:feedback\n"
        "                Restrict codewords to those that compare as feedback\n"
        "                to guess.  This option may be specified multiple times.\n"
        "    -h          Display this help screen and exit\n"
        "    -l level    (For make mode) Set the difficulty level to one of:\n"
        "                novice, easy, normal (default), hard, expert\n"
        "    -r rule     Set codeword rule, such as 6x4 (default), 10p4\n"
        "    -s secret   (For make mode) Set the secret\n"
        "action:\n"
        "    count       Display number of possible secrets\n"
        "    canonical   Display canonical guesses in lexicographical order\n"
        "    list        Display possible secrets in lexicographical order\n"
        "    make        Accept guesses from stdin and print feedbacks to stdout\n"
        "                until the secret is revealed.\n"
        "    test        Run the codebreaker against the codemaker.\n"
        ;
    return error ? 1 : 0;
}

void play(const mastermind::CodewordRules &rules)
{
    std::unique_ptr<Codemaker> code_maker(create_static_codemaker(sample(rules)));

    for (int round = 1; ; ++round)
    {
        std::cout << "Guess #" << round << "> ";
        std::cout.flush();

        std::string s;
        std::cin >> s;
        if (s == "quit")
            break;

        Codeword guess;
        try
        {
            if (!(std::istringstream(s) >> guess))
                throw std::invalid_argument("cannot read guess");
            if (!guess.conforms_to(rules))
                throw std::invalid_argument("guess does not conform to rules");
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

    std::unique_ptr<Codemaker> maker(create_static_codemaker(sample(rules)));
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
    CodewordPopulation population(rules);
    size_t total_steps = 0;
    size_t worst_steps = 0;
    for (size_t index = 0; index < population.size(); ++index)
    {
        std::unique_ptr<Codemaker> maker(create_static_codemaker(population.get(index)));
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

static void _display_canonical_guesses(const CanonicalCodewordSequence &sequence,
                                       std::span<const Codeword> candidate_guesses,
                                       size_t &counter)
{
    if (sequence.size() >= 2)
        return;

    std::vector<CanonicalCodewordSequence> extended_sequences =
        get_canonical_guesses(sequence, candidate_guesses);
    for (const CanonicalCodewordSequence &extended_sequence : extended_sequences)
    {
        ++counter;
        std::cout << counter << ">";
        for (auto cw : extended_sequence)
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
        _display_canonical_guesses(extended_sequence, candidate_guesses, counter);
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

int main(int argc, const char **argv)
{
    using namespace std::literals;

    CodewordRules rules;

#if 0
    display_canonical_guesses(rules);
    //return 0;
#endif

#if 0
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

    const char *action = nullptr;
    int level = 0;
    const char *heuristic = "minavg";
    std::optional<Codeword> secret;
    std::vector<Constraint> constraints;

    // Parse command line arguments.
    for (int i = 1; i < argc; i++)
    {
        if (action != nullptr)
            return usage("unexpected arguments after action");

        const char *opt = argv[i];
        if (opt[0] == '-')
        {
            if (opt[1] == '\0')
                return usage("missing option");
            if (opt[2] != '\0')
                return usage("invalid option");
            if (opt[1] == 'h')
                return usage(nullptr);
            else if (opt[1] == 'r')
            {
                if (++i >= argc)
                    return usage("missing argument");
                if (!(std::istringstream(argv[i]) >> rules))
                    return usage("invalid rules supplied to -r option");
            }
            else if (opt[1] == 's')
            {
                if (++i >= argc)
                    return usage("missing argument");
                Codeword w;
                if (!(std::istringstream(argv[i]) >> w))
                    return usage("invalid secret supplied to -s option");
                secret = w;
            }
            else if (opt[1] == 'l')
            {
                if (++i >= argc)
                    return usage("missing argument");
                const char *val = argv[i];
                if (val == "novice"sv)
                    level = -2;
                else if (val == "easy"sv)
                    level = -1;
                else if (val == "normal"sv)
                    level = 0;
                else if (val == "hard"sv)
                    level = 1;
                else if (val == "expert"sv)
                    level = 2;
                else
                    return usage("invalid level supplied to -l option");
            }
            else if (opt[1] == 'c')
            {
                if (++i >= argc)
                    return usage("missing argument");
                Constraint constraint;
                if (!(std::istringstream(argv[i]) >> constraint))
                    return usage("invalid constraint supplied to -c option");
                constraints.push_back(constraint);
            }
            else
                return usage("unknown option");
        }
        else
        {
            action = opt;
        }
    }

    (void)level;

    if (action == nullptr)
        return usage("missing action");

    CodewordPopulation population(rules);
    std::vector<Codeword> all(population.begin(), population.end());
    for (const Constraint &constraint : constraints)
    {
        if (!constraint.guess.conforms_to(rules))
            return usage("invalid constraint");
        auto it = std::stable_partition(all.begin(), all.end(), constraint);
        all.resize(it - all.begin());
    }

    if (action == "count"sv)
    {
        std::cout << all.size() << std::endl;
    }
    else if (action == "list"sv)
    {
        for (Codeword codeword : all)
        {
            std::cout << codeword << std::endl;
        }
    }
    else if (action == "test"sv)
    {
        std::unique_ptr<Codemaker> codemaker;
        if (secret.has_value())
            codemaker = create_static_codemaker(secret.value());
        else if (level == 0)
            codemaker = create_static_codemaker(sample(rules));
//        else if (level == -1)
        else
            return usage("level not supported");

        std::unique_ptr<CodeBreaker> codebreaker;
        codebreaker = create_heuristic_breaker(rules, heuristic);

        while (true)
        {
            Codeword guess = codebreaker->make_guess();
            Feedback feedback = codemaker->respond(guess);
            std::cout << Constraint{guess, feedback} << std::endl;
            if (feedback == Feedback::perfect_match(rules))
                break;
            codebreaker->step(guess, feedback);
        }
    }
    else
        return usage("invalid action");

//    std::cout << "Codeword rules: " << rules << std::endl;
//    std::cout << "  Alphabet size: " << rules.alphabet_size() << std::endl;
//    std::cout << "  Codeword size: " << rules.codeword_size() << std::endl;
//    std::cout << "  Is heterogram: " << std::boolalpha << rules.is_heterogram() << std::endl;
//    std::cout << "  Perfect match: " <<
//        Feedback::perfect_match(rules) << std::endl;

    return 0;
}
