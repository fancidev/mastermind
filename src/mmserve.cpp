/* mmserve.cpp - Mastermind codemaker */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "Mastermind.hpp"
using namespace Mastermind;

struct Constraint
{
	Codeword guess;
	Feedback response;

	Constraint() { }
	Constraint(const Codeword &_guess, const Feedback &_response)
		: guess(_guess), response(_response) { }
};

class Analyst
{
	Engine e;

	// List of all codewords.
	CodewordList all;

	// Stack of remaining possibilities corresponding to each constraint.
	std::vector<CodewordRange> _secrets;

	// Stack of constraints.
	std::vector<Constraint> _constraints;

public:

	explicit Analyst(const Rules &rules) 
		: e(rules), all(e.generateCodewords())
	{
		_secrets.push_back(CodewordRange(all));
	}

#if 0
	FeedbackFrequencyTable 
	void partition(const Codeword &guess)
	{
		FeedbackFrequencyTable freq = e.partition(_remaining, guess);

	}
#endif

	Engine& engine() { return e; }

	void push_constraint(const Codeword &guess, const Feedback &response)
	{
		// Partition the remaining possibilities.
		CodewordRange remaining = _secrets.back();
		CodewordPartition cells = e.partition(remaining, guess);
		remaining = cells[response.value()];

		// Update internal state.
		_constraints.push_back(Constraint(guess, response));
		_secrets.push_back(remaining);
	}

	void pop_constraint()
	{
		assert(!_constraints.empty());
		_constraints.pop_back();
		_secrets.pop_back();
	}

	/// Returns a list of remaining possibilities.
	CodewordConstRange possibilities() const 
	{
		return _secrets.back();
	}

	/// Returns a list of constraints.
	util::range<std::vector<Constraint>::const_iterator> constraints() const
	{
		return util::range<std::vector<Constraint>::const_iterator>
			(_constraints.begin(), _constraints.end());
	}
};

/// Displays help screen for player mode.
static void help()
{
	std::cout << 
		"Input your guess (e.g. 1234) or type one of the following commands:\n"
		"  !,cheat       show the secret\n"
		"  h,help        display this help screen\n"
		"  l,list        list remaining possibilities\n"
		"  q,quit        quit the program\n"
		"  r,recap       display guesses and responses so far\n"
		"";
}

/// Displays remaining possibilities.
static void list(CodewordConstRange secrets)
{
	for (CodewordConstIterator it = secrets.begin(); it != secrets.end(); ++it)
	{
		std::cout << *it << ' ';
	}
	std::cout << std::endl;
}

/// Displays current constraints.
static void recap(const Analyst &game)
{
	if (game.constraints().empty())
	{
		std::cout << "You haven't made any guess yet!\n";
	}
	else
	{
		for (size_t i = 0; i < game.constraints().size(); ++i)
		{
			const Constraint &c = game.constraints()[i];
			std::cout << c.guess << " " << c.response << std::endl;
		}
	}
}

/// Enters interactive player mode. This starts a prompt for user commands,
/// until either the user enter "quit" or reveals the secret.
/// Important output are put to STDOUT. 
/// Informational messages are put to STDOUT.
static int serve(const Engine *e, bool verbose, const Codeword &given_secret)
{
	Analyst game(e->rules());

	// Generate all codewords.
	CodewordList all = e->generateCodewords();
	if (verbose)
	{
		std::cout << "There are " << all.size() << " codewords. "
			<< "Please make guesses or type help for help." << std::endl;
	}

	// Generate a random secret if one is not specified.
	Codeword secret = given_secret;
	if (secret.IsEmpty())
	{
		srand((unsigned int)time(NULL));
		rand();
		int i = rand() % all.size();
		secret = all[i];
	}

	for (;;)
	{
		// Display prompt.
		if (verbose)
		{
			std::cout << "> ";
			std::flush(std::cout);
		}

		// Read a line of command.
		std::string line;
		std::getline(std::cin, line);
		std::istringstream ss(line);

		// Read command.
		std::string cmd;
		if (!(ss >> cmd))
			break;

		// If we're in interactive mode, we would accept a command.
		if (verbose)
		{
			if (cmd == "!" || cmd == "cheat")
			{
				if (verbose)
					std::cout << "Secret is ";
				std::cout << secret << std::endl;
				continue;
			}
			if (cmd == "h" || cmd == "help")
			{
				help();
				continue;
			}
			if (cmd == "l" || cmd == "list")
			{
				list(game.possibilities());
				continue;
			}
			if (cmd == "q" || cmd == "quit")
			{
				break;
			}
			if (cmd == "r" || cmd == "recap")
			{
				recap(game);
				continue;
			}
		}

		// Now we expect a guess from the input.
		Codeword guess;
		std::istringstream ss2(cmd);
		if (!(ss2 >> setrules(e->rules()) >> guess))
		{
			std::cout << "Invalid command or guess: " << cmd << std::endl;
			continue;
		}

		// Compute the response and display it.
		Feedback response = e->compare(guess, secret);
		if (verbose)
			std::cout << guess << " ";
		std::cout << response << std::endl;
		game.push_constraint(guess, response);
		if (response == Feedback::perfectValue(e->rules()))
			break;
	}
	return 0;
}

/// Displays command line usage information.
static void usage()
{
	std::cerr <<
		"Usage: mmserve [-r rules] [options]\n"
		"Serve as a codemaker for a Mastermind game.\n"
		"Rules: 'p' pegs 'c' colors 'r'|'n'\n"
		"    mm,p4c6r    [default] Mastermind (4 pegs, 6 colors, with repetition)\n"
		"    bc,p4c10n   Bulls and Cows (4 pegs, 10 colors, no repetition)\n"
		"    lg,p5c8r    Logik (5 pegs, 8 colors, with repetition)\n"
		"Options:\n"
		"    -h          display this help screen and exit\n"
		"    -i          interactive mode; display instructions\n"
		"    -u secret   use the given secret instead of generating a random one\n"
		"    -v          displays version and exit\n"
		"";
}

/// Displays version information.
static void version()
{
	std::cout << 
		"Mastermind Strategies Version " << MM_VERSION_MAJOR << "."
		<< MM_VERSION_MINOR << "." << MM_VERSION_TWEAK << std::endl
		<< "Configured with max " << MM_MAX_PEGS << " pegs and "
		<< MM_MAX_COLORS << " colors.\n"
		"Visit http://code.google.com/p/mastermind-strategy/ for updates.\n"
		"";
}

#define USAGE_ERROR(msg) do { \
		std::cerr << "Error: " << msg << ". Type -h for help." << std::endl; \
		return 1; \
	} while (0)

#define USAGE_REQUIRE(cond,msg) do { \
		if (!(cond)) USAGE_ERROR(msg); \
	} while (0)

int main(int argc, char* argv[])
{
	Rules rules(4, 6, true);
	bool verbose = false;
	Codeword secret;

	// Parse command line arguments.
	for (int i = 1; i < argc; i++)
	{
		std::string s = argv[i];
		if (s == "-h")
		{
			usage();
			return 0;
		}
		else if (s == "-i")
		{
			verbose = true;
		}
		else if (s == "-u")
		{
			USAGE_REQUIRE(++i < argc, "missing argument for option -u");
			std::string arg(argv[i]);
			std::istringstream ss(arg);
			USAGE_REQUIRE(ss >> setrules(rules) >> secret, "expecting secret after -u");
		}
		else if (s == "-v")
		{
			version();
			return 0;
		}
		else
		{
			USAGE_ERROR("unknown option: " << s);
		}
	}

	// Create an algorithm engine and begin the game.
	Engine engine(rules);
	serve(&engine, verbose, secret);
	return 0;
}
