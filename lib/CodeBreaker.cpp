#include <vector>
#include "CodeBreaker.hpp"
#include "ObviousStrategy.hpp"

namespace Mastermind {

// Create a free-standing function to make a guess.
Codeword MakeGuess(
	const Engine *e,
	CodewordConstRange secrets,
	Strategy *strat,
	const EquivalenceFilter *filter,
	const CodeBreakerOptions &options)
{
	size_t count = secrets.size();
	if (count == 0)
		return Codeword();

	// stat->NPossibilities = count;

	// Check for obvious guess.
	if (options.optimize_obvious)
	{
		Codeword guess = ObviousStrategy(e).make_guess(secrets, secrets);
		if (!guess.IsEmpty())
		{
			//stat->NCandidates = (it - possibilities.begin()) + 1;
			//stat->Guess = *it;
			return guess;
		}
	}

	// Initialize the set of candidate guesses.
	CodewordConstRange candidates = options.possibility_only ?
		secrets : e->universe();

	// Filter the candidate set to remove "equivalent" guesses.
	CodewordList canonical = filter->get_canonical_guesses(candidates);

	// Make a guess using the strategy provided.
	Codeword guess = strat->make_guess(secrets, canonical);
	return guess;
}

// Build a partial strategy tree from the given state.
// That state has NOT been output to the tree yet.
// The following fields must be filled in partial_node:
//   - node.depth
//   - node.guess
//   - node.feedback
// Note: This actually can be turned into a tail recursion.
//       We can do that later if necessary.
static void FillStrategy(
	StrategyTree &tree,            // the tree to fill
	StrategyTree::iterator where,  // iterator to the current state
	const Engine *e,               // algorithm engine
	unsigned char depth,           // number of guesses made so far
	const CodewordRange &secrets,  // list of remaining possibilities
	Strategy *strat,
	const EquivalenceFilter *filter,
	const CodeBreakerOptions &options,
	int *progress)
{
	// Make a guess.
	Codeword guess = MakeGuess(e, secrets, strat, filter, options);
	if (guess.IsEmpty())
		return;

	// Partition the possibility set using this guess.
	CodewordPartition cells = e->partition(secrets, guess);

	// Recursively fill the strategy for each possible response.
	Feedback perfect = Feedback::perfectValue(e->rules());
#if _OPENMP
	// index variable in OpenMP 'for' statement must have signed integral type
	#pragma omp parallel for schedule(dynamic)
#endif
	for (int k = 0; k < (int)cells.size(); ++k)
	{
		// Skip if this cell is empty.
		CodewordRange cell = cells[k];
		if (cell.empty())
			continue;

		// Create a subtree rooted from this guess/response pair.
		Feedback response(k);
		StrategyTree subtree(e->rules(), StrategyNode(guess, response)); // , depth + 1);

		if (response == perfect)
		{
			//(*progress)++;
			//double pcnt = (double)(*progress) / m_all.size();
			//printf("\rProgress: %3.0f%%", pcnt*100);
			//fflush(stdout);
		}
		else
		{
			// Create a child filter.
			std::unique_ptr<EquivalenceFilter> new_filter(filter->clone());
			new_filter->add_constraint(guess, response, cell);

			// Recursively build the strategy tree.
			FillStrategy(subtree, subtree.root(), e, depth + 1, cell, strat, 
				new_filter.get(), options, progress);
		}

		// Add the subtree to the big tree.
#if _OPENMP
		#pragma omp critical (CodeBreaker_FillStrategy)
#endif
		{
			tree.insert_child(where, subtree, true);
		}
	}
}

StrategyTree BuildStrategyTree(
	const Engine *e, 
	Strategy *strat, 
	const EquivalenceFilter *filter,
	const CodeBreakerOptions &options)
{
	CodewordList all = e->generateCodewords();

	StrategyTree tree(e->rules());

	int progress = 0;
	FillStrategy(tree, tree.root(), e, 0, all, strat, filter, options, &progress);
	return tree;
}

} // namespace Mastermind
