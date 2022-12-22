#include "rules.hpp"

#include <iostream>
#include <sstream>
#include <string>

static int usage(const char *error)
{
    if (error)
        std::cerr << "error: " << error << std::endl;
    std::cerr << "usage: mastermind -r rules" << std::endl;
    return error ? 1 : 0;
}

int main(int argc, const char *argv[])
{
    using namespace std::string_literals;

    mastermind::Rules rules;

    // Parse command line arguments.
    for (int i = 1; i < argc; i++)
    {
        const char *opt = argv[i];
        const char *val = (i + 1 < argc)? argv[i + 1] : nullptr;
        if (opt == "-r"s)
        {
            if (!val)
                return usage("option -r requires an argument");
            try
            {
                std::istringstream(val) >> rules;
            }
            catch (const std::invalid_argument &ex)
            {
                return usage(ex.what());
            }
        }
    }

    std::cout << "Using rules " << rules << " with "
        << rules.num_admissible() << " admissible codewords" << std::endl;
    return 0;
}
