#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <stack>
#include <vector>
#include <cassert>
#include <time.h>
#include <stdlib.h>

#include "Codeword.hpp"
#include "Algorithm.hpp"
#include "Engine.hpp"

using namespace Mastermind;

// List codewords.
static void listCodewords(CodewordConstRange list)
{
	for (CodewordConstIterator it = list.begin(); it != list.end(); ++it)
	{
		std::cout << *it << ' '; // std::endl;
	}
	std::cout << std::endl;
}

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

static void displayInfo(const Analyst &game, bool verbose = false)
{
	std::cout << game.constraints().size() << " constraints, "
		<< game.possibilities().size() << " remaining possibilities." 
		<< std::endl;
	if (verbose && !game.constraints().empty())
	{
		std::cout << "Constraints are:" << std::endl;
		for (size_t i = 0; i < game.constraints().size(); ++i)
		{
			const Constraint &c = game.constraints()[i];
			std::cout << c.guess << " " << c.response << std::endl;
		}
	}
}

// sort_order: 0 = no sort, 2 = lexical sort by output
static void displayPartitions(Analyst &game, CodewordConstRange guesses, 
	int sort_order = 2)
{
	if (guesses.empty())
		return;

	Engine &e = game.engine();

	// Create an array of feedbacks in the following order:
	// a04,a03,a02,a01,a00;a13,a12,a11,a10;a22,a21,a20;a30;a40
	int p = e.rules().pegs();
	std::vector<Feedback> header;
	for (int a = 0; a <= p; a++)
	{
		for (int b = p - a; b >= 0; b--)
		{
			if (a == p - 1 && b == 1)
				continue;
			header.push_back(Feedback(a, b));
		}
	}

	// Compute the list of partitions.
	size_t n = guesses.size();
	std::vector<FeedbackFrequencyTable> freqs(n);
	std::vector<size_t> index(n);
	for (size_t i = 0; i < n; ++i)
	{
		Codeword guess = guesses[i];
		freqs[i] = e.compare(guess, game.possibilities());
		index[i] = i;
	}

	// Sort the output.
	if (sort_order == 2)
	{
		auto comparer = [&header,&freqs](size_t i, size_t j) -> bool
		{
			const FeedbackFrequencyTable &t1 = freqs[i], &t2 = freqs[j];
			for (size_t k = 0; k < header.size(); ++k)
			{
				unsigned int f1 = t1[header[k].value()];
				unsigned int f2 = t2[header[k].value()];
				if (f1 < f2)
					return true;
				if (f2 < f1)
					return false;
			}
			return (i < j);
		};
		std::sort(index.begin(), index.end(), comparer);
	}

	// Display header.
	std::cout << "Guess " << ' ';
	for (size_t j = 0; j < header.size(); ++j)
	{
		std::cout << std::setw(5) << header[j];
	}
	std::cout << std::endl;

	// Display the partition of each guess.
	for (size_t i = 0; i < n; ++i)
	{
		Codeword guess = guesses[index[i]];
		const FeedbackFrequencyTable &freq = freqs[index[i]];

		std::cout << std::left << std::setw(6) << guess << ' ' << std::right;
		for (size_t j = 0; j < header.size(); ++j)
		{
			size_t f = freq[header[j].value()];
			if (f == 0)
				std::cout << std::setw(5) << '-';
			else
				std::cout << std::setw(5) << f;
		}
		std::cout << std::endl;
	}
}

/// Displays help screen for player mode.
static void help_player()
{
	std::cout << 
		"Input your guess (e.g. 1234) or type one of the following commands:\n"
		"  !,cheat       show the secret\n"
		"  h,help        display this help screen\n"
		"  i,info        display information\n"
		"  l,list        list remaining possibilities\n"
		"  q,quit,exit   quit the program\n"
		//"  r,recommend   let the computer recommend a guess\n"
		"";
}

/// Enters interactive player mode. This starts a prompt for user commands,
/// until either the user enter "quit" or reveals the secret.
/// Important output are put to STDOUT. 
/// Informational messages are put to STDOUT.
int interactive_player(const Engine *e, int verbose, const Codeword &_secret)
{
	Analyst game(e->rules());

	// Generate all codewords.
	CodewordList all = e->generateCodewords();
	if (verbose)
	{
		std::cout << "There are " << all.size() << " codewords. "
			<< "Please make guesses or type help for help." << std::endl;
	}

	// Generate a secret if one is not specified.
	Codeword secret = _secret;
	if (secret.empty())
	{
		srand((unsigned int)time(NULL));
		rand();
		int i = rand() % all.size();
		secret = all[i];
	}

	for (;;)
	{
		// Display prompt.
		std::cout << "> ";
		std::flush(std::cout);

		// Read a line of command.
		std::string line;
		std::getline(std::cin, line);
		std::istringstream ss(line);

		// Read command.
		std::string cmd;
		if (!(ss >> cmd))
			break;

		// Dispatch command.
		if (cmd == "q" || cmd == "quit" || cmd == "exit")
		{
			break;
		}
		else if (cmd == "!" || cmd == "cheat")
		{
			if (verbose)
				std::cout << "Secret is ";
			std::cout << secret << std::endl;
		}
		else if (cmd == "h" || cmd == "help")
		{
			help_player();
		}
		else if (cmd == "i" || cmd == "info")
		{
			displayInfo(game, true);
		}
		else if (cmd == "l" || cmd == "list")
		{
			listCodewords(game.possibilities());
		}
		else // guess
		{
			Codeword guess;
			std::istringstream ss(cmd);
			if (!(ss >> setrules(e->rules()) >> guess))
			{
				std::cout << "Invalid command or guess: " << cmd << std::endl;
				continue;
			}
#if 0
			if (!guess.conforming(e.rules()))
			{
				std::cout << "Invalid guess: " << guess << std::endl;
				continue;
			}
#endif
			Feedback response = e->compare(guess, secret);
			if (verbose)
				std::cout << guess << " ";
			std::cout << response << std::endl;
			game.push_constraint(guess, response);
			if (response == Feedback::perfectValue(e->rules()))
				break;
		}
	}
	return 0;
}

/// Display help screen for diagnostic mode.
static void help_analyst()
{
	std::cout << 
		"Commands:\n"
		"  +,push 1234 0a2b  push a constraint to the constraint stack\n"
		"  -,pop             pop the last constraint from the stack\n"
		"  c,canonical       list canonical guesses\n"
		"  e,eval [guess...] evaluate guesses by heuristic score of their partitions\n"
		"  h,help            display this help screen\n"
		"  i,info            display a summary of the current situation\n"
		"  l,list            list remaining possibilities\n"
		"  p,part [guess...] partition the remaining possibilities by each guess\n"
		"  q,quit,exit       quit the program\n"
		"";
}

int interactive_analyst(const Engine *e, int /* verbose */)
{
	Analyst game(e->rules());

	// Display available commands.
	help_analyst();

	// Generate all codewords.
	std::cout << "Generating all codewords..." << std::endl;
	CodewordList all = e->generateCodewords();
	std::cout << "Done. There are " << all.size() << " codewords." << std::endl;

	for (;;)
	{
		// Display prompt.
		std::cout << "> ";
		std::flush(std::cout);

		// Read a line of command.
		std::string line;
		std::getline(std::cin, line);
		std::istringstream ss(line);

		// Read command.
		std::string cmd;
		if (!(ss >> cmd) || cmd == "quit" || cmd == "exit" || cmd == "q")
			break;

		// Dispatch command.
		if (cmd == "push")
		{
			Codeword guess;
			Feedback response;
			if (!(ss >> setrules(e->rules()) >> guess >> response))
			{
				std::cout << "Expecting: push guess response\n";
				continue;
			}
#if 0
			if (!guess.valid(e.rules()))
			{
				std::cout << "Invalid guess: " << guess << std::endl;
				continue;
			}
			if (!response.valid(e.rules()))
			{
				std::cout << "Invalid response: " << response << std::endl;
				continue;
			}
#endif
			game.push_constraint(guess, response);
			displayInfo(game);
		}
		else if (cmd == "pop")
		{
			if (game.constraints().empty())
			{
				std::cout << "There are no constraints." << std::endl;
				continue;
			}
			game.pop_constraint();
			displayInfo(game);
		}
		else if (cmd == "list")
		{
			listCodewords(game.possibilities());
		}
		else if (cmd == "part")
		{
			CodewordList guesses;

			// If no guesses are given, display all canonical guesses.
			char c;
			if (!(ss >> c))
			{
				guesses = all;
			}
			else
			{
				ss.unget();
				for (Codeword guess; ss >> setrules(e->rules()) >> guess; )
				{
#if 0
					if (!guess.valid(e.rules()))
					{
						std::cout << "Invalid guess: " << guess << std::endl;
						continue;
					}
#endif
					guesses.push_back(guess);
				}
				if (ss.fail())
				{
					std::cout << "Expecting guess." << std::endl;
				}
			}
			displayPartitions(game, guesses);
		}
		else
		{
			std::cout << "Unknown command: " << cmd << std::endl;
		}
	}
	return 0;
}
