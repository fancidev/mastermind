#ifndef MASTERMIND_STRATEGY_TREE_HPP
#define MASTERMIND_STRATEGY_TREE_HPP

#include <iostream>
#include <vector>
#include <string>

#include "Rules.hpp"
#include "Codeword.hpp"
#include "Feedback.hpp"
#include "util/simple_tree.hpp"

namespace Mastermind {

/// Encapsulates information about a node on a strategy tree.
class StrategyNode
{
	Codeword::compact_type _guess;    // [unsigned long]
	Feedback::compact_type _response; // [unsigned char]

public:

	/// Constructs a node corresponding to the root state.
	StrategyNode() { }

	/// Constructs a node with the given guess and response.
	StrategyNode(const Codeword &guess, const Feedback &response)
		: _guess(guess.pack()), _response(response.pack())
	{
	}

	/// Returns the guess of the node.
	Codeword guess() const { return Codeword::unpack(_guess); }

	/// Returns the response of the node.
	Feedback response() const { return Feedback::unpack(_response); }
};

/**
 * Represents a guessing strategy of the code breaker.
 *
 * The strategy is internally stored as a tree. Each node in the tree
 * represents a <i>state</i> (or situation) of the game. Each state
 * is identified by a sequence of constraints consisting of the guesses
 * made and responses received so far. The constraints are read from
 * the root node to the given node.
 *
 * For example,
 * [Root] <- this state contains all the codewords
 *   [1234:0A0B] <- this state contains the parent state filtered by
 *                  1234:0A0B
 *     ...
 *     [3456:1A0B] <- this state contains the parent state filtered
 *                    by 3456:1A0B
 *     ...
 *
 * When we build a strategy tree, we use a recursive algorithm which always
 * appends a node or branch to the end of the tree. This is critical for
 * enabling a simple storage pattern. That is, we can store the nodes of the
 * tree linearly in a vector; the structure of the tree is identified by the
 * depth field associated with each node. This simple storage scheme makes
 * the implementation robust and efficient, though perhaps a bit less
 * straightforward than an implementation using linked-lists and pointers.
 *
 * A complete strategy tree starts from the initial state (with empty guess
 * and response); a partial strategy tree starts from a certain state of the
 * game (with at least one guess/response pairs). The root depth of the tree
 * or branch is also stored for integrity-checking purpose.
 */
class StrategyTree : public util::simple_tree<StrategyNode, int>
{
	Rules _rules;

public:

	/// Constructs a strategy tree (or branch) with the given rules.
	/// If <code>root.depth() == 0</code>, a full tree is constructed; otherwise,
	/// a branch is constructed.
	StrategyTree(
		const Rules &rules,
		const StrategyNode &root_data = StrategyNode()
#if 0
		, int root_depth = 0) : simple_tree(root_data, root_depth), _rules(rules)
#else
		)
#endif
		: util::simple_tree<StrategyNode,int>(root_data), _rules(rules)
	{ }

	//const std::string& name() const { return _name; }

	/// Returns the rules that the strategy applies to.
	const Rules& rules() const { return _rules; }
};

/// Encapsulates information of a strategy tree or branch.
class StrategyTreeInfo
{
	const StrategyTree *_tree;

	// Index of the first sub-state within this branch.
	StrategyTree::const_iterator _root;

	// Total number of secrets revealed in this branch. This is calculated as
	// the number of child states with a perfect response.
	unsigned int _total_secrets;

	// Total number of steps used to reveal all secrets in this branch.
	// This count starts from the root state of the branch, not from the
	// root of the strategy tree.
	size_t _total_depth;

	// depth_freq[i] = number of secrets revealed using (i) guesses.
	// depth_freq[0] will always be zero and is ignored.
	std::vector<unsigned int> _depth_freq;

	// _children[j] = index of the child state corresponding to response j,
	// or 0 if that response is not available.
	std::vector<StrategyTree::const_iterator> _children;

	std::string _name;

public:

	/// Gather information from a branch of a strategy tree.
	StrategyTreeInfo(
		const std::string &name,
		const StrategyTree &tree,
		StrategyTree::const_iterator root);

#if 1
	/// Returns the name.
	std::string name() const { return _name; }
#endif

	/// Returns the suggested guess for the root state.
	Codeword suggestion() const
	{
		auto children = _tree->children(_root);
		if (children.empty())
			return Codeword();
		else
			return children.begin()->guess();
	}

#if 1
	/// Returns a node iterator to the child state corresponding to making 
	/// the suggested guess receiving the given response. If the response
	/// is impossible, returns a null iterator.
	StrategyTree::const_iterator child(Feedback response) const
	{
		return _children[response.value()];
	}
#endif

	/// Returns the maximum depth of all leaf nodes.
	int max_depth() const
	{
		return (int)_depth_freq.size() - 1;
	}

	/// Returns the number of leaf nodes with the given depth.
	unsigned int count_depth(int depth) const
	{
		return (depth >= 0 && depth <= max_depth())?
			_depth_freq[depth] : 0;
	}

	/// Returns the sum of the depths of all leaf nodes.
	size_t total_depth() const { return _total_depth; }

	/// Returns the total number of secrets revealed in this branch.
	unsigned int total_secrets() const { return _total_secrets; }

	/// Returns the average number of guesses used to reveal each secret.
	double average_depth() const
	{
		return (double)total_depth() / total_secrets();
	}
};

/// Outputs strategy tree statistics to a stream.
std::ostream& operator << (std::ostream &os, const StrategyTreeInfo &info);

#if 0
/// Defines the file format to serialize a stratey tree.
enum FileFormat
{
	DefaultFormat = 0,
	TextFormat = 0,
	XmlFormat = 1,
	BinaryFormat = 2,
};
#endif

/// Outputs a strategy tree in text format (Irving convention).
void WriteStrategy_TextFormat(std::ostream &os, const StrategyTree &tree);

/// Outputs a strategy tree in XML format.
void WriteStrategy_XmlFormat(std::ostream &os, const StrategyTree &tree);

/// Deserializes a strategy tree from a stream.
/// If the input format is invalid, the stream's fail bit is set,
/// and an error message is written to cerr.
std::istream& operator >> (std::istream &is, StrategyTree &tree);

} // namespace Mastermind

#endif // MASTERMIND_STRATEGY_TREE_HPP
