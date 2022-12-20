#include <iostream>
#include <fstream>
#include <stdio.h>
#include <malloc.h>
#include <iomanip>
#include <memory>

#include "Rules.hpp"
#include "Codeword.hpp"
#include "Feedback.hpp"
#include "Algorithm.hpp"
#include "Engine.hpp"
#include "CodeBreaker.hpp"
#include "SimpleStrategy.hpp"
#include "HeuristicStrategy.hpp"
#include "OptimalStrategy.hpp"
#include "Equivalence.hpp"
#include "Heuristics.hpp"
#include "Benchmark.hpp"

#include "util/call_counter.hpp"
#include "util/hr_timer.hpp"
#include "util/io_format.hpp"

using namespace Mastermind;

void pause_output()
{
#ifdef _WIN32
	system("PAUSE");
#endif
}

#if 0
static void simulate_guessing(
	Engine &e, Strategy* strats[], size_t n,
	EquivalenceFilter *filter,
	const CodeBreakerOptions &options)
{
	CodewordList all = e.generateCodewords();
	Rules rules = e.rules();

	std::cout
		<< "Game Settings" << std::endl
		<< "---------------" << std::endl
		<< "Number of pegs:      " << rules.pegs() << std::endl
		<< "Number of colors:    " << rules.colors() << std::endl
		<< "Color repeatable:    " << std::boolalpha << rules.repeatable() << std::endl
		<< "Number of codewords: " << rules.size() << std::endl;

	// Pick a secret "randomly".
	Codeword secret = all[all.size()/4*3];
	std::cout << std::endl;
	std::cout << "Secret: " << secret << std::endl;

	// Use an array to store the status of each codebreaker.
	std::vector<bool> finished(n);

	// Create an array of codebreakers and output strategy names.
	// BUG: when changing from CodeBreaker* to CodeBreaker, there is
	// a memory error. Find out why.
	std::vector<CodeBreaker*> breakers;
	std::cout << std::left << " # ";
	for (size_t i = 0; i < n; ++i)
	{
		std::string name = strats[i]->name();
		std::cout << std::setw(10) << name;
		breakers.push_back(new CodeBreaker(e,
			std::unique_ptr<Strategy>(strats[i]),
			std::unique_ptr<EquivalenceFilter>(filter->clone()),
			options));
	}
	std::cout << std::right << std::endl;

	// Output horizontal line.
	std::cout << "---";
	for (size_t i = 0; i < n; ++i)
	{
		std::cout << "----------";
	}
	std::cout << std::endl;

	// Step-by-step guessing.
	int step = 0;
	for (size_t finished_count = 0; finished_count < n; )
	{
		std::cout << std::setw(2) << (++step);
		std::cout.flush();

		// Test each code breaker in turn.
		for (size_t i = 0; i < n; ++i)
		{
			if (finished[i])
			{
				std::cout << "          ";
				continue;
			}

			CodeBreaker &breaker = *breakers[i];
			Codeword guess = breaker.MakeGuess();
			if (guess.empty())
			{
				std::cout << " FAIL";
				finished[i] = true;
				++finished_count;
			}
			else
			{
				Feedback feedback = e.compare(secret, guess);
				std::cout << " " << guess << ":" << feedback;
				std::cout.flush();
				if (feedback == Feedback::perfectValue(e.rules()))
				{
					finished[i] = true;
					++finished_count;
				}
				breaker.AddConstraint(guess, feedback);
			}
		}
		std::cout << std::endl;
	}
}
#endif

#if 1
static void test_strategy_tree(
	const Engine *e,
	Strategy *strategies[],
	size_t n,
	EquivalenceFilter *filter,
	const CodeBreakerOptions &options)
	// const Codeword& first_guess)
{
	//CodewordList all = e.generateCodewords();
	Rules rules = e->rules();
	//Feedback target = Feedback::perfectValue(rules);
	
	std::cout
		<< "Game Settings" << std::endl
		<< "---------------" << std::endl
		<< "Number of pegs:      " << rules.pegs() << std::endl
		<< "Number of colors:    " << rules.colors() << std::endl
		<< "Color repeatable:    " << std::boolalpha << rules.repeatable() << std::endl
		<< "Number of codewords: " << rules.size() << std::endl;

	std::cout << std::endl
		<< "Options" << std::endl
		<< "---------" << std::endl
		<< "Optimize obvious guess: " << std::boolalpha
			<< options.optimize_obvious << std::endl
		<< "Guess possibility only: " << std::boolalpha
			<< options.possibility_only << std::endl;


	std::cout << std::endl << util::header;

	for (size_t i = 0; i < n; ++i)
	{
		Strategy *strat = strategies[i];

		// Build a strategy tree of this code breaker
		EquivalenceFilter *copy = filter->clone();
		StrategyTree tree = BuildStrategyTree(e, strat, copy, options);
		delete copy;

		StrategyTreeInfo info(strat->name(), tree, tree.root());
		std::cout << info;
	}
}
#endif

#if 0
static void display_canonical_guesses(
	Engine &e,
	const EquivalenceFilter *filter,
	int max_level,
	int level = 0)
{
	CodewordConstRange candidates = e.universe();
	CodewordList canonical = filter->get_canonical_guesses(candidates);

	// Display each canonical guess, and expand one more level is needed.
	if (level >= max_level)
	{
		std::cout << "[" << level << ":" << canonical.size() << "]";
#if 1
		if (canonical.size() > 20)
		{
			std::cout << " ... " << std::endl;
		}
		else
#endif
		{
			for (size_t i = 0; i < canonical.size(); ++i)
			{
				Codeword guess = canonical[i];
				std::cout << " " << guess;
			}
			std::cout << std::endl;
		}
	}
	else
	{
		for (size_t i = 0; i < canonical.size(); ++i)
		{
			Codeword guess = canonical[i];
			std::cout << "[" << level << ":" << i << "] " << guess << std::endl;

			EquivalenceFilter *child = filter->clone();
			child->add_constraint(guess, Feedback(), candidates);
			display_canonical_guesses(e, child, max_level, level+1);
			delete child;
			//std::cout << "[" << level << "] Total: " << canonical.size() << std::endl;
		}
	}
}

static void test_initial_guesses_in_equivalence_filter(
	Engine &e, const EquivalenceFilter *filter)
{
	display_canonical_guesses(e, filter, 2);
	delete filter;
}

static void display_partition_size(
	Engine &e, CodewordRange secrets, const Codeword &guess)
{
	CodewordPartition parts = e.partition(secrets, guess);
	for (size_t i = 0; i < parts.size(); ++i)
	{
		if (i > 0)
			std::cout << ",";
		std::cout << parts[i].size();
	}
}

static void test_partition_size(Engine &e, const char *name)
{
	auto filter = std::unique_ptr<EquivalenceFilter>(
		RoutineRegistry<CreateEquivalenceFilterRoutine>::get(name)(e));

	// Get an initial canonical guess.
	CodewordList all = e.generateCodewords();
	CodewordList initials = filter->get_canonical_guesses(e.universe());
	Codeword g1 = initials[0];
	// std::cout << g1 << std::endl;

	// Partition the codeword set by g1.
	CodewordPartition parts = e.partition(all, g1);
	for (size_t i = 0; i < parts.size(); ++i)
	{
		CodewordRange part = parts[i];
		if (part.empty())
			continue;

		// Add constraint.
		auto child = std::unique_ptr<EquivalenceFilter>(filter->clone());
		child->add_constraint(g1, Feedback(i), part);
		CodewordList canonical = child->get_canonical_guesses(e.universe());
		std::cout << "Feedback: " << Feedback(i) << ", #remaining = "
			<< part.size() <<  ", #{g2} = " << canonical.size() << std::endl;

		if (1)
		{
			// Display how each guess partitions the remaining secrets.
			for (size_t j = 0; j < canonical.size() && j < 5; ++j)
			{
				std::cout << canonical[j] << ": ";
				display_partition_size(e, part, canonical[j]);
				std::cout << std::endl;
			}
		}
	}
}

#endif

static void test_serialization()
{
	Rules rules(4, 6, true);
	StrategyTree tree(rules);

	if (1)
	{
		std::cout << "File 1" << std::endl;
		std::ifstream fs("strats/p4c6r-koyama-1.txt");
		fs >> tree;
	}
	if (1)
	{
		std::cout << "File 2" << std::endl;
		std::ifstream fs("strats/p4c6r-koyama-2.txt");
		fs >> tree;
	}
	if (1)
	{
		std::cout << "File 3" << std::endl;
		std::ifstream fs("strats/p4c6r-knuth.txt");
		fs >> tree;
	}
	if (1)
	{
		std::cout << "File 4" << std::endl;
		std::ifstream fs("strats/p4c6r-neuwirth.txt");
		fs >> tree;
	}
	if (1)
	{
		std::cout << "File 5" << std::endl;
		std::ifstream fs("strats/p4c6r-irving.txt");
		fs >> tree;
	}
}

/// Runs regression and benchmark tests.
int test(const Rules &rules, bool /* verbose */)
{
#ifdef NDEBUG
#define LOOP_FLAG 1
#else
#define LOOP_FLAG 0
#endif

#if 0
	extern int interactive(const Rules &rules);
	return interactive(rules);
#endif

	// Set up the standard engine.
	Engine engine(rules);
	const Engine *e = &engine;

#if 1
	compareRoutines<ComparisonRoutine>(e, "generic", "norepeat", 100000*LOOP_FLAG);
	//compareRoutines<ComparisonRoutine>(e, "generic", "test", 100000*LOOP_FLAG);

	//compareRoutines<GenerationRoutine>(e, "generic", "generic", 100*LOOP_FLAG);
	//compareRoutines<MaskRoutine>(e, "generic", "unrolled", 100000*LOOP_FLAG);
	//testSumSquares(rules, "generic", "generic", 10000000*LOOP_FLAG);
	//return TestFrequencyCounting(rules, 250000*LOOP_FLAG);
	//return TestEquivalenceFilter(rules, 10000*LOOP_FLAG);
	system("PAUSE");
	return 0;
#endif

	// Choose an equivalence filter.
	std::unique_ptr<EquivalenceFilter> filter(
#if 0
		//RoutineRegistry<CreateEquivalenceFilterRoutine>::get("Dummy")(e)
		//RoutineRegistry<CreateEquivalenceFilterRoutine>::get("Constraint")(e)
		RoutineRegistry<CreateEquivalenceFilterRoutine>::get("Color")(e)
#else
		new CompositeEquivalenceFilter(
			RoutineRegistry<CreateEquivalenceFilterRoutine>::get("Color")(e),
			RoutineRegistry<CreateEquivalenceFilterRoutine>::get("Constraint")(e))
#endif
		);

#if 0
	test_initial_guesses_in_equivalence_filter(e, filter.get());
	system("PAUSE");
	return 0;
#endif

#if 0
	test_partition_size(e, filter_name);
	system("PAUSE");
	return 0;
#endif

#if 0
	extern void test_optimal_strategy(Engine &);
	test_optimal_strategy(e);

	std::cout << std::endl;
	std::cout << "Call statistics for OptimalRecursion:" << std::endl;
	std::cout << util::call_counter::get("OptimalRecursion") << std::endl;

#if 0
	std::cout << "Call statistics for Comparison:" << std::endl;
	std::cout << util::call_counter::get("Comparison") << std::endl;
#endif

	pause_output();
	return 0;
#endif

#if 1
	test_serialization();
	pause_output();
	return 0;
#endif

	// @todo output strategy tree size info.
	CodeBreakerOptions options;
	options.optimize_obvious = true;
	options.possibility_only = false;

	//Codeword first_guess = Codeword::emptyValue();
	//Codeword first_guess = Codeword::Parse("0011", rules);

	using namespace Mastermind::Heuristics;
	Strategy* strats[] = {
		new SimpleStrategy(),
		new HeuristicStrategy<MinimizeWorstCase>(e),
#if 0
		new HeuristicStrategy<MinimizeAverage>(e),
		new HeuristicStrategy<MaximizeEntropy<false>>(e),
		new HeuristicStrategy<MaximizeEntropy<true>>(e),
		new HeuristicStrategy<MaximizePartitions>(e),
		new HeuristicStrategy<MinimizeLowerBound>(e, MinimizeLowerBound(e))
#endif
		//new OptimalCodeBreaker(rules),
	};
	int nstrat = sizeof(strats)/sizeof(strats[0]);

	// Note: Only one of the following can be invoked because
	// a code breaker makes a unique pointer of the strategy.
	//simulate_guessing(e, strats, nstrat, filter.get(), options);
	test_strategy_tree(e, strats, nstrat, filter.get(), options);

#if 0
	std::cout << "Call statistics for ConstraintEquivalence:" << std::endl;
	std::cout << util::call_counter::get("ConstraintEquivalence") << std::endl;

	std::cout << "Call statistics for ConstraintEquivalencePermutation:" << std::endl;
	std::cout << util::call_counter::get("ConstraintEquivalencePermutation") << std::endl;

	std::cout << "Call statistics for ConstraintEquivalenceCrossout:" << std::endl;
	std::cout << util::call_counter::get("ConstraintEquivalenceCrossout") << std::endl;
#endif

#if 0
	// Display some statistics.
	std::cout << "Call statistics for Comparison:" << std::endl;
	std::cout << util::call_counter::get("Comparison") << std::endl;
#endif

#if 0
	if (0) {
		printf("\nRun again:\n");
		CountFrequenciesImpl->SelectRoutine("c_p8_il_os");
		TestGuessingByTree(rules, breakers, sizeof(breakers)/sizeof(breakers[0]), first_guess);
		printf("\n");
	}
#endif

	pause_output();
	return 0;
}

#if 0
static int RegressionTestUnit()
{
	Rules rules;
	rules.length = 4;
	rules.ndigits = 10;
	rules.allow_repetition = false;

	// Codeword creation
	Codeword cw;
	if (!cw.IsEmpty()) {
		printf("FAILED: Codeword::IsEmpty()\n");
		return -1;
	}

	// Enumeration
	CodewordList list;
	list = CodewordList::Enumerate(rules);
	if (!(list.GetCount() == 5040)) {
		printf("FAILED: Wrong enumeration count\n");
		return -1;
	}
	if (!(list[357].ToString() == "0741")) {
		printf("FAILED: Wrong enumeration result\n");
		return -1;
	}

	// Comparison
	FeedbackList fbl(list[0], list);
	if (!(fbl.GetCount() == 5040)) {
		printf("FAILED: Codeword::CompareTo() returns wrong length of feedback\n");
		return -1;
	}
	if (!(fbl[3] == Feedback(3, 0))) {
		printf("FAILED: Incorrect comparison results\n");
		return -1;
	}

	// Count frequency
	return 0;
}

static int RegressionTest()
{
	printf("Running regression test...\n");
	if (RegressionTestUnit() == 0) {
		printf("Passed\n");
	} else {
		printf("Failed\n");
	}

	system("PAUSE");
	return 0;
}
#endif
