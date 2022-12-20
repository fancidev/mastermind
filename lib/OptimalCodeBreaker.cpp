#include <iostream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <vector>
#include <array>
#include <functional>
#include <numeric>

#include "Engine.hpp"
#include "Strategy.hpp"
#include "Equivalence.hpp"
#include "ObviousStrategy.hpp"
#include "HeuristicStrategy.hpp"
#include "OptimalStrategy.hpp"
#include "StrategyTree.hpp"
#include "util/call_counter.hpp"
#include "util/hr_timer.hpp"
#include "util/io_format.hpp"

using namespace Mastermind;

/**
 * Searches for an obviously optimal strategy for the given remaining secrets.
 *
 * @returns The cost of the obvious strategy if one is found; otherwise zero.
 */
static StrategyCost fill_obviously_optimal_strategy(
	const Engine *e,       // algorithm engine
	CodewordRange secrets, // list of remaining secrets
	StrategyObjective obj,
	StrategyConstraints c,
	StrategyTree &tree,    // Strategy tree that stores the best strategy
	StrategyTree::iterator where // iterator to the current state
	)
{
	// @todo Take into account the lower-bound estimate of non-possible
	// guesses.
	StrategyCost _cost;
	Codeword guess = make_obvious_guess(e, secrets, c.max_depth, obj, _cost, obj);
	if (guess.IsEmpty())
		return StrategyCost();

	//	VERBOSE_COUT << "Found obvious guess: " << obvious << std::endl;

	// Automatically fill the strategy tree using this guess.This requires
	// all cells in the partition to have no more than two possibilities.
	// This is equivalent to Knuth's 'x' notation in writing a strategy.
	Feedback perfect = Feedback::perfectValue(e->rules());
	FeedbackList fbs;
	e->compare(guess, secrets, fbs);
	size_t n = secrets.size();
	unsigned int cost = 0;

	for (size_t j = 0; j < Feedback::size(e->rules()); ++j)
	{
		Codeword first;
		StrategyTree::iterator it;
		for (size_t i = 0; i < n; ++i)
		{
			if (fbs[i] == Feedback(j))
			{
				if (first.IsEmpty())
				{
					++cost;
					first = secrets[i];
					it = tree.insert_child(where, StrategyNode(guess, fbs[i]));
					if (fbs[i] != perfect)
					{
						++cost;
						tree.insert_child(it, StrategyNode(first, perfect));
					}
				}
				else
				{
					cost += 3;
					tree.insert_child(
						tree.insert_child(it, StrategyNode(first, e->compare(secrets[i], first))),
						StrategyNode(secrets[i], perfect));
				}
			}
		}
	}
	assert(cost == _cost.steps);
	return _cost;
}

#define VERBOSE_COUT(text) WRAP_STATEMENTS( \
	if (verbose) { \
		std::cout << std::setw(depth*2) << "" << "[" << (depth+1) << "] " \
			<< text << std::endl; \
	} )

typedef HeuristicStrategy<Heuristics::MinimizeLowerBound> LowerBoundEstimator;

/// @todo: Use heuristic code breaker to estimate an upper bound of total guesses
/// This could be helpful in pruning obvious bad candidates

/**
 * Searches for an optimal strategy for the given set of remaining secrets.
 *
 * @param depth Depth of the current state. This is equal to the number of 
 *      guesses already made.
 * @param threshold Branch pruning threshold. If the cost of revealing all
 *      the secrets would reach or exceed this threshold, the function will
 *      fail. This parameter should be set to the currently best cost achieved
 *      by some known strategy.
 * @returns The cost of the optimal strategy if one is found, or zero if 
 *      one is not found either because some secret would require more than
 *      <code>c.max_depth</code> guesses to reveal, or because the cost of
 *      any strategy would reach or exceed the cut-off threshold.
 * @todo Add progress report (0% - 100%)
 */
// @todo: We might change the equivalence filter interface to operate on
// the input inplace? This could save a few memory copies but may change
// the output.
// all the secrets, or -1 if such optimal will not be less than _best_.
static StrategyCost fill_strategy_tree(
	const Engine *e,
	CodewordRange secrets,            // remaining secrets; will be partitioned
	CodewordRange candidates,         // canonical guesses; may be sorted
	const EquivalenceFilter *filter1, // response-independent equivalence filter
	const EquivalenceFilter *filter2, // response-dependent equivalence filter
	LowerBoundEstimator &estimator,   // lower bound estimator
	const int depth,                  // depth of the current state; root=0
	StrategyObjective obj,            // objective
	StrategyConstraints c,            // constraints
	StrategyCost threshold,           // prunes branch if cost >= threshold
	StrategyTree &tree,               // tree to store the strategy found
	StrategyTree::iterator where      // iterator to the current state
	)
{
	UPDATE_CALL_COUNTER("OptimalRecursion", (int)secrets.size());

	bool verbose = false; // (depth < 1);

	VERBOSE_COUT("Checking " << secrets.size() << " remaining secrets");

	// Fail if the secret set is empty or the max-depth is 0.
	if (secrets.empty() || c.max_depth == 0)
		return StrategyCost();

	// Initialize common variables.
	const Feedback perfect = Feedback::perfectValue(e->rules());
	const unsigned int nsecrets = (int)secrets.size();

	// Short-cut if there is only one secret.
	if (nsecrets == 1)
	{
		tree.insert_child(where, StrategyNode(secrets[0], perfect));
		return StrategyCost(1, 1, 1);
	}

	// From now on, we will need to make at least one guess to reveal any 
	// secret, and at least two guesses to reveal all secrets. This accounts
	// for n total steps and 1 extra step. 
	
	// Update the constraints and threshold so that they will apply to a 
	// situation after the initial guess.
	if (c.max_depth == 1)
		return StrategyCost();
	else
		--c.max_depth;

	if (threshold.steps <= nsecrets)
		return StrategyCost();
	else
		threshold.steps -= nsecrets;

	// @bug The following criteria is problematic.
	if (threshold.depth <= 1)
		return StrategyCost();
	else
		--threshold.depth;

	// Define a strategy cost comparer.
	StrategyCostComparer superior(obj);

#if 0
	// If find_last is true, then we only cut-off a guess when it's
	// strictly worse than the current best; otherwise, we cut-off
	// a guess as long as it's no better than the current best.
	int cut_off_delta = (c.find_last? 1 : 0);
	cut_off += cut_off_delta ???
#endif

	// Filter canonical guesses as candidates.
	//CodewordList candidates = filter->get_canonical_guesses(e.universe());
	//VERBOSE_COUT("Found " << candidates.size() << " canonical guesses.");

	// Compute a lower bound of the cost for each candidate guess, exluding
	// the cost of making the initial guess because this is the same for all
	// candidates. Then sort the candidates by this lower bound so that more
	// "promising" candidates are processed first. This helps to improve the
	// upper bound as early as possible.
	// @todo It might be better to rename scores to extra_cost.
	typedef Heuristics::MinimizeLowerBound::score_t lowerbound_t;
	std::vector<lowerbound_t> scores(candidates.size());
	//estimator.make_guess(secrets, candidates, scores.data());
	estimator.evaluate(secrets, candidates, scores.data());

	// @todo We might opt to remove the need to create an index array.
	// Instead, we could scan for the element in each iteration.
	std::vector<int> order(candidates.size());
	std::iota(order.begin(), order.end(), 0);

	// Define SORT_CANDIDATES to 1 to explicitly sort the candidate guesses.
	// Since many guesses will be pruned right away (especially if we have
	// a good estimate of the lower-bound of the cost), it is usually faster
	// to (linearly) search for the smallest element in each iteration,
	// instead of sorting the whole array at the beginning. However, the
	// results will be the same as long as we perform a stable sort.
#define SORT_CANDIDATES 0

#if SORT_CANDIDATES
	std::sort(order.begin(), order.end(), [&](int i, int j) -> bool {
		if (superior(scores[i], scores[j]))
			return true;
		if (superior(scores[j], scores[i]))
			return false;
		return i < j;
	});
#endif

	// Initialize state variables to store the best guess and its cost so far.
	StrategyCost best;
	StrategyTree best_tree(e->rules());

	// Try each candidate guess.
	size_t candidate_count = candidates.size();
	for (size_t index = 0; index < candidate_count; ++index)
	{
#if !SORT_CANDIDATES
		// Find the guess with the lowest estimated cost in the remaining
		// candidates, and swap it to the front.
		auto min_it = std::min_element(order.begin() + index, order.end(),
			[&](int i, int j) -> bool { 
				return superior(scores[i], scores[j]); 
		});
		std::swap(*min_it, order[index]);
#endif
		size_t i = order[index];
		Codeword guess = candidates[i];

		// Since we keep improving the upper bound dynamically,
		// and we sort the candidates by their lower bound,
		// we need to check here whether the remaining candidates
		// are still worth checking.
		if (!superior(scores[i], threshold))
		{
			VERBOSE_COUT("Pruned " << (candidate_count - index)
				<< " remaining guesses: lower bound (" << scores[i]
				<< ") >= cut-off (" << threshold << ")");
#if 0
			if (i == 0)
			{
				std::cout << "Very first pruned: best lower bound = "
					<< scores[i] << ", cut-off = " << cut_off << std::endl;
			}
#endif
			break;
		}

		VERBOSE_COUT("Checking guess " << (i+1) << " of "
			<< candidate_count << " (" << guess << ") -> ");

		// If there's a limit on the maximum number of guesses allowed,
		// check if we can prune it.
		if (scores[i].depth > c.max_depth)
		{
			if (verbose)
				std::cout << "Skipped: guess will have too many steps"
				<< std::endl;
			continue;
		}

#if 0
		static size_t test_counter = 0;
		++test_counter;
		if (test_counter == 8479)
		{
			__debugbreak();
		}
#endif

		// Partition the remaining secrets using this guess.
		// Note that after successive calls to @c partition,
		// the order of the secrets are shuffled. However,
		// that should not impact the optimality of the result.
		CodewordPartition cells = e->partition(secrets, guess);

		// Sort the partitions by their size, so that smaller partitions
		// (i.e. smaller search trees) are processed first. This helps
		// to improve the lower bound (slack) at an earlier stage.
		std::array<int,Feedback::MaxOutcomes> responses;
		size_t nresponses = cells.size();
		std::iota(responses.begin(), responses.begin() + nresponses, 0);
		std::sort(responses.begin(), responses.begin() + nresponses,
			[&cells](int i, int j) -> bool
		{
			if (cells[i].size() == 0)
				return false;
			if (cells[j].size() == 0)
				return true;
			if (cells[i].size() < cells[j].size())
				return true;
			if (cells[j].size() < cells[i].size())
				return false;
			return i < j;
		});

		// Find the number of availble responses, and skip this guess if it 
		// generates only one response.
		while (nresponses > 0 && cells[responses[nresponses-1]].empty())
			--nresponses;
		if (nresponses <= 1)
		{
			if (verbose)
				std::cout << "Skipped: guess produces unit partition"
				<< std::endl;
			continue;
		}

		// Estimate a lower bound of the cost of revealing the secrets in 
		// each partition, NOT counting the cost of making the initial guess.
		// If the total lower bound reaches or exceeds the cut-off threshold,
		// we can prune this guess.
		// Note: this step is redundant because we have already calculated
		// the same score before.
		StrategyCost lb_part[256];
		StrategyCost lb;
		for (size_t j = 0; j < nresponses; ++j)
		{
			Feedback feedback = Feedback(responses[j]);
			if (feedback != perfect)
			{
				lowerbound_t estimate =
					estimator.heuristic().simple_estimate((int)cells[feedback.value()].size());
				lb_part[j] = estimate;
				lb += lb_part[j];
			}
		}

		// Since this is a double-computation, we shouldn't be pruning
		// this guess here.
		assert(lb == scores[i]);

		if (verbose)
		{
			std::cout << nresponses << " cells:";
			for (size_t j = 0; j < nresponses; ++j)
			{
				if (j > 0)
					std::cout << ',';
				std::cout << cells[responses[j]].size();
			}
			std::cout << "; lower bound = " 	<< lb 
				<< ", best = " << best << std::endl;
		}

		// Find the best guess for each partition. We adopt a two-phase
		// equivalence filtering method. First, we filter all possible
		// codewords by (response-indepedent) constraint equivalence.
		// This can be done once for all response classes. Then, for each
		// individual response class, we apply the response-dependent 
		// color equivalence filter.
		CodewordList pre_filtered;
		std::unique_ptr<EquivalenceFilter> pre_filter(filter1->clone());
		pre_filter->add_constraint(guess, Feedback(), e->universe());
		// @todo we may change the interface of add_constraint to return
		// a new filter.

		bool pruned = false;
		StrategyTree this_tree(e->rules());
		for (size_t j = 0; j < nresponses && !pruned; ++j)
		{
			Feedback feedback = Feedback(responses[j]);
			const CodewordRange &cell = cells[feedback.value()];

			// Add this node to the strategy tree.
			StrategyNode node(guess, feedback);
			StrategyTree::iterator it = this_tree.insert_child(this_tree.root(), node);

			// Do not recurse for a perfect match.
			if (feedback == perfect)
			{
				VERBOSE_COUT("- Checking cell " << feedback
					<< " -> perfect");
				continue;
			}

			VERBOSE_COUT("- Checking cell " << feedback
				<< " -> lower bound = " << lb_part[j]);

			// Short-cut if only one additional guess is allowed
			// but we are left with more than one secret in this
			// partition.
			// @todo such pruning could be improved and consolidated with
			//  the pruning in the beginning of the routine.
			if (c.max_depth == 1 && cell.size() > 1)
			{
				pruned = true;
				break;
			}

			// If there's an obviously optimal guess for this cell, use it.
			// @todo "mastermind -v -s optimal -r mm -md 5" doesn't
			// respect the "-md 5" option.
			StrategyCost cell_cost = fill_obviously_optimal_strategy(
				e, cell, obj, c, this_tree, it);
			if (!!cell_cost)
			{
				//VERBOSE_COUT("- Checking cell " << cell.feedback
				//	<< " -> found obvious guess");
				VERBOSE_COUT("  Found obvious guess");
			}
			else
			{
				//VERBOSE_COUT("- Checking cell " << cell.feedback
				//	<< " -> lower bound = " << lb);
				//VERBOSE_COUT("- Checking cell " << cell.feedback
				//	<< " -> lower bound = " << lb_part[j]);

				// @todo: estimate a lower bound of the cost of any guess.
				// If the lower bound is greater than the cut-off, we don't
				// need to proceed any more.

				// Apply constraint filter on the candidate guesses if not 
				// already done so. This filter does not depend on the response,
				// so a single run can be used for all response classes.
				if (pre_filtered.empty())
				{
					if (c.pos_only)
						pre_filtered = pre_filter->get_canonical_guesses(secrets);
					else
						pre_filtered = pre_filter->get_canonical_guesses(e->universe());
				}

				// Apply color filter on the pre-filtered candidates.
				std::unique_ptr<EquivalenceFilter> new_filter(filter2->clone());
				new_filter->add_constraint(guess, feedback, cell);
				CodewordList canonical = new_filter->get_canonical_guesses(pre_filtered);

				// @todo: Check this. The minus sign doesn't work for complex
				// cost structure.
				cell_cost = fill_strategy_tree(e, cell, canonical,
					pre_filter.get(), new_filter.get(), estimator,
					depth + 1, obj, c, threshold - (lb - lb_part[j]),
					this_tree, it);
			}

			if (!cell_cost) // No strategy was found for this cell
			{
				VERBOSE_COUT("Pruned this guess because the recursion returns -1.");
				pruned = true;
				break;
			}

#if 1
			if (superior(lb_part[j], cell_cost))
				VERBOSE_COUT("  Cell optimal cost is " << cell_cost);
			else if (superior(cell_cost, lb_part[j]))
				VERBOSE_COUT("  ERROR: LOWER BOUND IS HIGHER THAN ACTUAL COST.");
			else
				VERBOSE_COUT("  Lower bound unchanged.");
#endif

			// Refine the lower bound estimate.
			// @bug: the lower bound estimate needs to be amended
			lb += (cell_cost - lb_part[j]);
			lb_part[j] = cell_cost;
			if (!superior(lb, threshold))
			{
				VERBOSE_COUT("Skipping " << (nresponses-j-1) << " remaining "
					<< "partitions because lower bound (" << lb << ") >= best ("
					<< best << ")");
				pruned = true;
				break;
			}
		}

		// Now the guess is either pruned, or is the best guess so far.
		if (!pruned)
		{
			assert(superior(lb, threshold));
			assert(!best || superior(lb, best));
			best = lb;
			threshold = best;
			VERBOSE_COUT("Improved cut-off to " << best);

			// Use this_tree as best_tree.
			// @todo: optimize swap()
			std::swap(this_tree, best_tree);
		}
	}

	// If a best strategy was found, append the strategy to the tree.
	// Also, since 'best' is calculated without accounting for the initial
	// guess, we need to add it back.
	if (!!best)
	{
		tree.insert_child(where, best_tree, false);
		best.steps += nsecrets;
		++best.depth;
		// best.worst;
	}

	return best;
}

StrategyTree build_optimal_strategy_tree(
	const Engine *e, StrategyObjective obj, StrategyConstraints constraints)
{
	CodewordList all = e->generateCodewords();

	// Creates a composite equivalence filter by chaining a
	// response-indepedent filter with a response-dependent filter.
	CompositeEquivalenceFilter filter(
		CreateConstraintEquivalenceFilter(e),
		CreateColorEquivalenceFilter(e));

	// Create a strategy tree.
	StrategyTree tree(e->rules());

	// Set options.
	// StrategyObjective obj = min_depth ? MinDepth : MinSteps;
	//StrategyConstraints c;
	//c.max_depth = (unsigned char)std::min(100, max_depth);
	//c.find_last = false;

	// Create a cost lower-bound estimator.
	LowerBoundEstimator estimator(e, Heuristics::MinimizeLowerBound(e));

	// Filter canonical candidates for the initial guess.
	CodewordList initial = filter.get_canonical_guesses(e->universe());

	// Recursively find an optimal strategy.
	StrategyCost threshold(1000000, 100, 0);
	/* int best = */ fill_strategy_tree(e, all, initial, 
		filter.first(), filter.second(), estimator,
		0, obj, constraints, threshold, tree, tree.root());
	// std::cout << "OPTIMAL: " << best << std::endl;
	return tree;
}

// Call statistics for optimal Mastermind (p4c6r) that finds the FIRST:
// Total # of calls : 5832
// Total # of ops   : 59209
// Avg ops per call : 10.15
// Time: 1.21 s
// Total   Avg    1    2    3    4    5    6
//  5625 4.340    1    8   92  648  542    5
//
// Call statistics for optimal Mastermind (p4c6r) that finds the LAST:
// Total # of calls : 12029
// Total # of ops   : 89738
// Avg ops per call : 7.46
// Time: 2.45
// Total   Avg    1    2    3    4    5    6
//  5625 4.340    1    8   90  648  548    1
//
// This gives a lower bound and upper bound of the time needed.
// And FIND_LAST happens to return a strategy with the least worst-case
// secrets.

#if 0
void test_optimal_strategy(Engine &e)
{
	util::hr_timer t1;
	t1.start();
	StrategyTree tree = 	build_optimal_strategy_tree(e);
	double t = t1.stop();

	StrategyTreeInfo info("optimal", tree, t, tree.root());
	std::cout << std::endl << util::header << info;
}
#endif
