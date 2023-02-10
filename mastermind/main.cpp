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
        "    -c guess:feedback\n"
        "                Restrict codewords to those that compare as feedback\n"
        "                to guess.  This option may be specified multiple times.\n"
        "    -h          Display this help screen and exit\n"
        "    -l level    (For play mode) Set the difficulty level to one of:\n"
        "                novice, easy, normal (default), hard, expert\n"
        "    -r rule     Set codeword rule, such as 6x4 (default), 10p4\n"
        "action:\n"
        "    count       Display the number of codewords\n"
        "    canonical   Display canonical guesses (in lexical order)\n"
        "    list        Display all codewords (in lexical order)\n"
        "    play        Make guesses from stdin and get feedbacks from stdout\n"
        "                until the secret is revealed.  See also the -l option.\n"
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
    std::vector<std::pair<Codeword, Feedback>> constraints;

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
                Codeword guess;
                Feedback feedback;
                char sp;
                if (!(std::istringstream(argv[i]) >> guess >> sp >> feedback))
                    return usage("invalid constraint supplied to -c option");
                if (sp != ':')
                    return usage("invalid constraint supplied to -c option");
                constraints.push_back({guess, feedback});
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
    for (const std::pair<Codeword, Feedback> &constraint : constraints)
    {
        if (!constraint.first.conforms_to(rules))
            return usage("invalid constraint");
        auto it = std::stable_partition(
            all.begin(), all.end(), [&constraint](const Codeword &secret) {
            return compare(secret, constraint.first) == constraint.second;
        });
        // TODO: make constraint a callable object
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
    else
        return usage("invalid action");

//    std::cout << "Codeword rules: " << rules << std::endl;
//    std::cout << "  Alphabet size: " << rules.alphabet_size() << std::endl;
//    std::cout << "  Codeword size: " << rules.codeword_size() << std::endl;
//    std::cout << "  Is heterogram: " << std::boolalpha << rules.is_heterogram() << std::endl;
//    std::cout << "  Perfect match: " <<
//        Feedback::perfect_match(rules) << std::endl;

//    std::cout << "First 5:";
//    for (size_t index = 0; index < 5 && index < all.size(); index++)
//    {
//        std::cout << " " << all[index];
//    }
//    std::cout << std::endl;
//
//    std::cout << "Last  5:";
//    for (size_t index = std::max<size_t>(5, all.size()) - 5; index < all.size(); index++)
//    {
//        std::cout << " " << all[index];
//    }
//    std::cout << std::endl;

//    Codeword guess = Codeword::from_string("1357", rules);
//    Codeword secret = Codeword::from_string("2337", rules);
//    std::cout << "compare(" << guess.to_string("ABCDEFGHIJ") << ", "
//        << secret.to_string("ABCDEFGHIJ") << ") = "
//        << compare(guess, secret).to_string() << std::endl;

    return 0;
}
