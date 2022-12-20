/**
 * @defgroup type Data Types
 *
 * @defgroup algo Basic Operations
 *
 * @defgroup equiv Equivalence Filters
 *
 * @defgroup strat Strategies
 *
 * @defgroup prog Program
 *
 * @defgroup util Utility Routines
 */

/// @defgroup Rules Rules
/// @ingroup type

/// @defgroup Codeword Codeword
/// @ingroup type

/// @defgroup Feedback Feedback
/// @ingroup type

/// @defgroup Obvious Obvious Strategy
/// A strategy that makes an obviously-optimal guess when one exists.
/// @ingroup strat

/// @defgroup Simple Simple Strategy
/// A strategy that makes an arbitrary guess from the set of remaining
/// secrets.
/// @ingroup strat

/// @defgroup Heuristic Heuristic Strategy
/// A strategy that makes the guess that produces the best heuristic
/// score.
/// @ingroup strat

/// @defgroup Optimal Optimal Strategy
/// A strategy that makes the optimal guess through exhaustive search.
/// @ingroup strat

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include "Rules.hpp"
#include "Codeword.hpp"
#include "Equivalence.hpp"
#include "SimpleStrategy.hpp"
#include "HeuristicStrategy.hpp"
#include "OptimalStrategy.hpp"
#include "CodeBreaker.hpp"
#include "Heuristics.hpp"
#include "util/io_format.hpp"

using namespace Mastermind;

#if 0

// c.f. http://www.javaworld.com.tw/jute/post/view?bid=35&id=138372&sty=1&tpg=1&ppg=1&age=0#138372

static int GetMaxBreakableWithin(
	int nsteps,
	int f[],
	CodewordList possibilities,
	CodewordList all,
	unsigned short unguessed_mask,
	int expand_levels)
{
	assert(nsteps >= 0);
	assert(expand_levels >= 0);

	if (nsteps == 0) {
		return 0;
	} else if (nsteps == 1) {
		return 1;
	}
	//else if (nsteps == 2) {
	//	int npegs = possibilities.GetRules().length;
	//	return std::min(possibilities.GetCount(), npegs*(npegs+3)/2);
	//}

	unsigned short impossible_mask = ((1<<all.GetRules().ndigits)-1) & ~possibilities.GetDigitMask();
	CodewordList candidates = all.FilterByEquivalence(unguessed_mask, impossible_mask);
	int best = 0;
	for (int i = 0; i < candidates.GetCount(); i++) {
		Codeword guess = candidates[i];

		int nguessable = 0;
		if (expand_levels == 0) {
			FeedbackFrequencyTable freq(FeedbackList(guess, possibilities));
			for (int fbv = 0; fbv <= freq.GetMaxFeedbackValue(); fbv++) {
				int partition_size = freq[Feedback(fbv)];
				nguessable += std::min(f[nsteps-1], partition_size);
			}
		} else {
			FeedbackFrequencyTable freq;
			possibilities.Partition(guess, freq);
			int partition_start = 0;
			for (int fbv = 0; fbv <= freq.GetMaxFeedbackValue(); fbv++) {
				int partition_size = freq[Feedback(fbv)];
				if (partition_size == 0)
					continue;

				CodewordList filtered(possibilities, partition_start, partition_size);
				int tmpcount = GetMaxBreakableWithin(nsteps-1, f, filtered, all,
					unguessed_mask & ~guess.GetDigitMask(),
					expand_levels-1);
				nguessable += std::min(tmpcount, partition_size);
				partition_start += partition_size;
			}
		}

		best = std::max(best, nguessable);
	}
	return best;
}

static int TestBound(Rules rules)
{
	// Find out the maximum size of any partition
	CodewordList all = CodewordList::Enumerate(rules);
	FeedbackFrequencyTable freq(FeedbackList(all[0], all));

	int f[100];
	f[0] = 0;
	f[1] = 1;
	int min_rounds = 0;
	int expand_levels = 1;
	for (int n = 2; n < 100; n++) {
		unsigned short unguessed_mask = ((1<<all.GetRules().ndigits)-1);
		f[n] = GetMaxBreakableWithin(n, f, all.Copy(), all, unguessed_mask, expand_levels);
		if (f[n] >= all.GetCount()) {
			min_rounds = n;
			break;
		}
	}
	for (int n = 1; n <= min_rounds; n++) {
		printf("Maximum number of secrets that can be guessed in %2d rounds: <= %d\n",
			n, f[n]);
	}
	printf("**** Minimum number of rounds necessary to guess all secrets: >= %d\n",
		min_rounds);

	int min_steps = all.GetCount()*min_rounds;
	for (int n = 1; n < min_rounds; n++) {
		min_steps -= f[n];
	}
	printf("**** Minimum total steps necessary to guess all secrets: >= %d\n",
		min_steps);
	printf("**** Minimum average steps necessary to guess all secrets: >= %5.3f\n",
		(double)min_steps/all.GetCount());

	int *S = Mastermind::Heuristics::MinimizeSteps::partition_score;
	S[0] = 0;
	for (int i = 1; i <= all.GetCount(); i++) {
		int nr = 0;
		for (nr = 1; nr <= 100; nr++) {
			if (f[nr] >= i)
				break;
		}
		int score = i*nr;
		for (int j = 1; j < nr; j++) {
			score -= f[j];
		}
		S[i] = score;
		//if (i % 40==0)
		//printf("S[%d] = %d\n", i, score);
	}

#if 0
	int m=14;

	const int nmax=20;
	int S[nmax+1];
	S[0]=0;
	int n;
	for (n = 1; n <= nmax; n++) {
		if (n <= m) {
			S[n] = 1+(n-1)*2;
		} else {
			S[n]=n;
			int nper = (n-1)/(m-1);
			S[n] += S[nper]*(m-2);
			S[n] += S[(n-1)-(nper*(m-2))];
		}
		printf("S[%4d] = %d\n", n, S[n]);
	}
#endif

	//system("PAUSE");
	return 0;
}
#endif

static void usage()
{
	std::cerr <<
		"Usage: mmstrat [-r rules] -s strategy [options]\n"
		"Build the specified strategy for the given rules.\n"
		"Rules: 'p' pegs 'c' colors 'r'|'n'\n"
		"    mm,p4c6r    [default] Mastermind (4 pegs, 6 colors, with repetition)\n"
		"    bc,p4c10n   Bulls and Cows (4 pegs, 10 colors, no repetition)\n"
		"    lg,p5c8r    Logik (5 pegs, 8 colors, with repetition)\n"
#ifndef NDEBUG
		"Modes:\n"
		"    -d          interactive diagnostics\n"
		"    -p [secret] interactive player (optionally using the given secret)\n"
		"    -s strat    build strategy 'strat' and output strategy tree\n"
		"    -t          run tests\n"
#endif
		// @todo descriptions for heuristic strategies 
		"Strategies:\n"
#ifndef NDEBUG
		"    replay path replay the strategy in 'path'; use - for STDIN\n"
#endif
		"    simple      simple strategy\n"
		"    minmax      min-max heuristic strategy\n"
		"    minavg      min-average heuristic strategy\n"
		"    entropy     max-entropy heuristic strategy\n"
		"    parts       max-parts heuristic strategy\n"
#ifndef NDEBUG
		"    minlb       min-lowerbound heuristic strategy\n"
#endif
		"    optimal     optimal strategy\n"
		"General Options:\n"
		"    -h          display this help screen and exit\n"
#ifdef _OPENMP
		"    -mt [n]     enable parallel execution with n threads [default="
		<< omp_get_max_threads() << "]\n"
#endif
		"    -po         make guess from remaining possibilities only\n"
#if ENABLE_CALL_COUNTER
		"    -prof       collect and display profiling details before exit\n"
#endif
		"    -q          quiet mode; display minimal information\n"
		"    -S          output strategy summary instead of strategy tree\n"
		"    -v          displays version and exit\n"
		"Options for Heuristic Strategies:\n"
		"    -e filter   specify the equivalence filter to use, which is one of:\n"
		"                default     composite filter (color + constraint)\n"
		"                color       filter by color equivalence\n"
		"                constraint  filter by constraint equivalence\n"
		"                none        do not apply any filter\n"
		"    -nc         do not apply a correction to the heuristic score\n"
		"                which favors guesses from remaining possibilities.\n" 
		"    -no         Do not attempt to make an obvious guess before applying\n"
		"                the heuristic function. This option is useful for debugging\n"
		"                purpose if the heuristic function may yield a guess that\n"
		"                is different than an obvious guess when one exists.\n"
		"Options for Optimal Strategies:\n"
#ifndef NDEBUG
		"    -md depth   limit the maximum number of guesses to reveal any secret\n"
#endif
		"    -O level    specify the level of optimization, which is one of:\n"
		"                1 - (default) minimize steps\n"
#ifndef NDEBUG
		"                2 - minimize steps, then depth\n"
		"                3 - minimize steps, then depth, then worst count\n"
#endif
		"";
}

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

// TODO: Add progress display to OptimalCodeBreaker
// TODO: Output strategy tree after finishing a run

extern int interactive_player(const Engine *e, int verbose, const Codeword &secret);
extern int interactive_analyst(const Engine *e, int verbose);
extern int test(const Rules &rules, bool verbose);

extern void pause_output();

#define USAGE_ERROR(msg) do { \
		std::cerr << "Error: " << msg << ". Type -h for help." << std::endl; \
		return 1; \
	} while (0)

#define USAGE_REQUIRE(cond,msg) do { \
		if (!(cond)) USAGE_ERROR(msg); \
	} while (0)

static int build_heuristic_strategy_tree(
	const Engine *e, const EquivalenceFilter *filter, int /* verbose */,
	const std::string &name, StrategyConstraints constraints,
	bool no_correction, StrategyTree &tree)
{
	using namespace Mastermind::Heuristics;

	bool ac = !no_correction; // apply correction
	Strategy *strat = NULL;
	if (name == "simple")
		strat = new SimpleStrategy();
	else if (name == "minmax")
		strat = new HeuristicStrategy<MinimizeWorstCase>(e, MinimizeWorstCase(ac));
	else if (name == "minavg")
		strat = new HeuristicStrategy<MinimizeAverage>(e, MinimizeAverage(ac));
	else if (name == "entropy")
		strat = new HeuristicStrategy<MaximizeEntropy>(e, MaximizeEntropy(ac));
	else if (name == "parts")
		strat = new HeuristicStrategy<MaximizePartitions>(e, MaximizePartitions(ac));
	else if (name == "minlb")
		strat = new HeuristicStrategy<MinimizeLowerBound>(e, MinimizeLowerBound(e));
	else
		USAGE_ERROR("unknown strategy: " << name);

	CodeBreakerOptions options;
	options.optimize_obvious = (name == "simple")? false : constraints.use_obvious;
	options.possibility_only = constraints.pos_only;
	std::unique_ptr<EquivalenceFilter> copy(filter->clone());
	tree = BuildStrategyTree(e, strat, copy.get(), options);
	return 0;
}

extern StrategyTree build_optimal_strategy_tree(
	const Engine *e, StrategyObjective obj, StrategyConstraints constraints);

// verbose: 0 = quiet, 1 = verbose, 2 = very verbose
static int build_strategy(
	const Engine *e, const EquivalenceFilter *filter, int verbose,
	const std::string &name, const std::string & /* file */,
	StrategyConstraints constraints, bool no_correction,
	StrategyObjective obj, bool summary)
{
	using namespace Mastermind::Heuristics;

	StrategyTree tree(e->rules());

	if (name == "file")
	{
		USAGE_ERROR("Not implemented");
	}
	else if (name == "optimal")
	{
		tree = build_optimal_strategy_tree(e, obj, constraints);
	}
	else
	{
		int ret = build_heuristic_strategy_tree(e, filter, verbose, name,
			constraints, no_correction, tree);
		if (ret)
			return ret;
	}

	// Output result.
	if (summary)
	{
		StrategyTreeInfo info(name, tree, tree.root());
		if (verbose)
		{
			std::cout << util::header;
			std::cout << info;
		}
		else
		{
			// @todo The following output should be output directly from
			// a StrategyCost object.
			std::cout << info.total_depth() << ':' << info.max_depth() << ':' 
				<< info.count_depth(info.max_depth()) << std::endl;
		}
	}
	else
	{
		WriteStrategy_TextFormat(std::cout, tree);
	}
	
	return 0;
}

int main(int argc, char* argv[])
{
	using namespace Mastermind;

#if 0
	Rules rules(4, 3, true);
#elif 0
	Rules rules(4, 10, false);
#else
	Rules rules(4, 6, true);
#endif

	int verbose = 1;
	enum Mode
	{
		DefaultMode = 0,
		StrategyMode = 1,
		PlayerMode = 2,
		DebugMode = 3,
		TestMode = 4
	} mode = DefaultMode;
	std::string strat_name, strat_file, filter_name;
	Codeword secret;
#ifdef _OPENMP
	int mt = 1;
#endif
	StrategyConstraints constraints;
	StrategyObjective obj = MinSteps;
	bool prof = false; // whether to enable profiling (call counting)
	bool no_correction = false;
	bool summary = false;

	// Parse command line arguments.
	for (int i = 1; i < argc; i++)
	{
		std::string s = argv[i];
		if (s == "-d")
		{
			USAGE_REQUIRE(mode == DefaultMode, "only one mode may be specified");
			mode = DebugMode;
		}
		else if (s == "-p")
		{
			USAGE_REQUIRE(mode == DefaultMode, "only one mode may be specified");
			mode = PlayerMode;
			if (i+1 < argc && argv[i+1][0] != '-')
			{
				std::string arg(argv[++i]);
				std::istringstream ss(arg);
				USAGE_REQUIRE(ss >> setrules(rules) >> secret,
					"expecting secret after -p");
			}
		}
		else if (s == "-O")
		{
			USAGE_REQUIRE(++i < argc, "missing argument for option -O");
			std::string level = argv[i];
			if (level == "1")
				obj = MinSteps;
			else
				USAGE_ERROR("invalid optimization level '" << level << "'");
		}
		else if (s == "-S")
		{
			summary = true;
		}
		else if (s == "-s")
		{
			USAGE_REQUIRE(mode == DefaultMode, "only one mode may be specified");
			USAGE_REQUIRE(++i < argc, "missing argument for option -s");
			mode = StrategyMode;
			strat_name = argv[i];
			if (strat_name == "file")
			{
				USAGE_REQUIRE(++i < argc, "missing input filename for file strategy");
				strat_file = argv[i];
			}
		}
		else if (s == "-t")
		{
			USAGE_REQUIRE(mode == DefaultMode, "only one mode may be specified");
			mode = TestMode;
		}
		else if (s == "-e")
		{
			USAGE_REQUIRE(filter_name.empty(), "only one equivalence filter may be specified");
			USAGE_REQUIRE(++i < argc, "missing argument for option -f");
			filter_name = argv[i];
		}
		else if (s == "-h")
		{
			usage();
			return 0;
		}
		else if (s == "-no")
		{
			constraints.use_obvious = false;
		}
		else if (s == "-nc")
		{
			no_correction = true;
		}
		else if (s == "-q")
		{
			verbose = 0;
		}
		else if (s == "-r")
		{
			USAGE_REQUIRE(++i < argc, "missing argument for option -r");
			USAGE_REQUIRE(secret.empty(), "-r rules must be specified before -p secret");
			std::string name = argv[i];
			if (name == "mm")
				rules = Rules(4, 6, true);
			else if (name == "bc")
				rules = Rules(4, 10, false);
			else if (name == "lg")
				rules = Rules(5, 8, true);
			else
				rules = Rules(name.c_str());
			USAGE_REQUIRE(rules, "invalid rules: " << argv[i]);
		}
		else if (s == "-md")
		{
			USAGE_REQUIRE(++i < argc, "missing argument for option -md");
			std::string cnt(argv[i]);
			int max_depth;
			USAGE_REQUIRE((std::istringstream(cnt) >> max_depth) && (max_depth > 0),
				"positive integer argument expected for option -md");
			constraints.max_depth = (unsigned char)std::min(100, max_depth);
		}
		else if (s == "-po")
		{
			constraints.pos_only = true;
		}
		else if (s == "-prof")
		{
#if ENABLE_CALL_COUNTER
			prof = true;
#else
			std::cerr << "Warning: option -prof is not supported by this build"
				" and is ignored." << std::endl;
#endif
		}
		else if (s == "-mt")
		{
			int n = -1;
			if (i+1 < argc && argv[i+1][0] != '-')
			{
				std::string cnt(argv[++i]);
				USAGE_REQUIRE((std::istringstream(cnt) >> n) && (n > 0),
					"positive integer argument expected for option -mt");
			}
#ifdef _OPENMP
			if (n < 0)
			{
				mt = omp_get_max_threads();
			}
			else
			{
				if (n > omp_get_max_threads())
				{
					std::cerr << "Warning: number of threads set to maximum value "
						<< omp_get_max_threads() << std::endl;
					n = omp_get_max_threads();
				}
				mt = n;
			}
#else
			std::cerr << "Warning: option -mt is not supported by this build"
				" and is ignored." << std::endl;
#endif
		}
		else if (s == "-v")
		{
			version();
			return 0;
		}
#if 0
		else if (s == "-vv")
		{
			verbose = 2;
		}
#endif
		else
		{
			USAGE_REQUIRE(false, "unknown option: " << s);
		}
	}

	// Set number of threads.
#ifdef _OPENMP
	omp_set_num_threads(mt);
	omp_set_nested(0);
#endif

	// Enables or disables profiling according to -prof switch.
	util::call_counter::enable(prof);

	// Create an algorithm engine.
	Engine engine(rules);
	const Engine *e = &engine;

	// Create the specified equivalence filter.
	EquivalenceFilter *filter = NULL;
	if (filter_name == "default" || filter_name == "")
	{
		filter = new CompositeEquivalenceFilter(
			RoutineRegistry<CreateEquivalenceFilterRoutine>::get("Color")(e),
			RoutineRegistry<CreateEquivalenceFilterRoutine>::get("Constraint")(e));
	}
	else if (filter_name == "color")
	{
		filter = RoutineRegistry<CreateEquivalenceFilterRoutine>::get("Color")(e);
	}
	else if (filter_name == "constraint")
	{
		filter = RoutineRegistry<CreateEquivalenceFilterRoutine>::get("Constraint")(e);
	}
	else if (filter_name == "none")
	{
		filter = RoutineRegistry<CreateEquivalenceFilterRoutine>::get("Dummy")(e);
	}
	else
	{
		USAGE_ERROR("unknown equivalence filter: " << filter_name);
	}
	std::unique_ptr<EquivalenceFilter> filter_obj(std::move(filter));

	// Execute the specified action.
	int ret = 0;
	switch (mode)
	{
	case StrategyMode:
		ret = build_strategy(e, filter, verbose, strat_name, strat_file, 
			constraints, no_correction, obj, summary);
		break;
	case PlayerMode:
		ret = interactive_player(e, verbose, secret);
		break;
	case DebugMode:
		ret = interactive_analyst(e, verbose);
		break;
	case TestMode:
		ret = test(rules, verbose ? true : false);
		break;
	default:
		USAGE_ERROR("missing mode");
	}

	// Display available profiling results. It is useful to disgard the 
	// profiling switch here to detect any code that doesn't respect the
	// switch.
#if 0
	if (prof)
#endif
	{
		if (prof)
		{
			std::cout << std::endl << "**** Profiling Details ****" << std::endl;
		}
		auto cr = util::call_counter::registry();
		for (auto it = cr.begin(); it != cr.end(); ++it)
		{
			if (it->second.total_calls() > 0)
				std::cout << it->second << std::endl;
		}
	}

	return ret;
}
