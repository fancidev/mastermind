#include <iostream>
#include <iomanip>
#include <cassert>
#include <sstream>

#include <algorithm>
#include <utility>
#include <map>
#include "util/io_format.hpp"

#include "StrategyTree.hpp"
#include "Algorithm.hpp"
#include "Engine.hpp"

namespace Mastermind {

///////////////////////////////////////////////////////////////////////////
// StrategyTreeInfo implementation

StrategyTreeInfo::StrategyTreeInfo(
	const std::string &name,
	const StrategyTree &tree,
	StrategyTree::const_iterator root)
	: _tree(&tree), _root(root), _total_secrets(0), _total_depth(0),
	_depth_freq(1), 	_children(Feedback::size(tree.rules())), _name(name)
{
	Feedback perfect = Feedback::perfectValue(tree.rules());
	auto root_depth = root.depth();
	auto nodes = tree.traverse(root);
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		if (it.depth() == root_depth + 1)
		{
			_children[it->response().value()] = it;
		}
		if (it->response() == perfect)
		{
			unsigned int d = it.depth() - root_depth;
			if (d >= _depth_freq.size())
			{
				_depth_freq.resize(d+1);
			}
			++_depth_freq[d];
			++_total_secrets;
			_total_depth += d;
		}
	}
}

std::ostream& operator << (std::ostream &os, const StrategyTreeInfo &info)
{
	const int max_display = 9;

	// Display header.
	if (util::hasheader(os))
	{
		//os << "Stategy Statistics" << std::endl;
		//os << "--------------------" << std::endl;
		os << "Strategy   Total   Avg   1     2     3     4     5     6     7     8    >8" << std::endl;
		os << util::noheader;
	}

	// Display label.
	os  << std::left << std::setw(8) << info.name() << std::right
		<< std::setw(8) << info.total_depth() << " "
		<< std::setw(5) << std::setprecision(3)
		<< std::fixed << info.average_depth();

	// Display each frequency.
	unsigned int running_total = 0;
	for (int d = 1; d <= max_display; d++)
	{
		unsigned int count;
		if (d < max_display)
		{
			count = info.count_depth(d);
			running_total += d * count;
		}
		else
		{
			count = (unsigned int)info.total_depth() - running_total;
		}
		os << std::setw(d == 1? 4 : 6);
		if (count > 0)
			os << count;
		else
			os << '-';
	}

	os << std::endl;

	return os;
}

///////////////////////////////////////////////////////////////////////////
// StrategyTree output (Irving notation)

/*
Example:
1296 (1123: 2 (2311), 44A, 222B, 276C, 81D; 4 (1312), 84E, 230F, 182G; 5 (1321: 1, 0, 0, 0,
0; 2 (1132), 0, 0, 0; 1, 0, 0; 0; 1), 40H, 105I; 20J; 1)

A = (2345: 0, 6 (3412), 10 (3612: 1, 0, 0, 0, 0; 1, 1, 2 (4211), 0; 1, 0, 2 (3411); 1; 1),
    2 (3611), 0; 0, 6 (3215), 8 (2461: 0, 1, 1,0, 0; 0, 0, 1, 1; 0, 1, 2 (2231); 1; 0), 4 (2611); 0,
    2 (2314), 4 (2316); 2 (2315); 0)
*/
void WriteState_TextFormat(
	std::ostream &os, const StrategyTree &tree, StrategyTree::const_iterator root, 
	std::map<char,std::string> &symbols, int symbol_level)
{
	const Rules rules = tree.rules();

	// Get info about this strategy branch.
	StrategyTreeInfo state("", tree, root);

	// Write the number of remaining secrets in this state.
	os << state.total_secrets();

	// If there is zero or one secret, we are done.
	if (state.max_depth() <= 1)
		return;

	// If all secrets can be revealed with one extra guess, we output 
	// the obvious guess without elaborating the responses to it.
	// In addition, if this guess is not among the remaining possibilities,
	// append an asterisk.
	// @todo Use Knuth's notation and compress a guess that partitions
	// remaining possibilities into cells with size <= 2.
	if (state.max_depth() <= 2)
	{
		os << " (" << state.suggestion();
		if (!state.child(Feedback::perfectValue(rules)))
			os << '*';
		os << ")";
		return;
	}

	// Now we need to elaborate the sub-strategy for each possible response.
	// If labelling is allowed, we will create a symbol for this.
	int p = rules.pegs();
	bool use_symbol = (symbol_level == 1) && (root > 0) && 
		(state.total_secrets() >= (unsigned int)(p*(p+3)/2));
	std::ostringstream strs;
	std::ostream &ss = use_symbol ? strs : os;

	ss << " (" << state.suggestion() << ":";
	for (int a = 0; a <= p; ++a)
	{
		for (int b = p - a; b >= 0; --b)
		{
			if (a == p - 1 && b == 1) // skip 3A1B
				continue;
			ss << ' ';

			Feedback response(a, b);
			StrategyTree::const_iterator child = state.child(response);
			if (!child)
			{
				ss << 0;
			}
			else if (a == p)
			{
				ss << 1;
			}
			else
			{
				WriteState_TextFormat(ss, tree, child, symbols, symbol_level-1);
			}
			ss << ((b > 0)? ',' : (a == p)? ')' : ';');
		}
	}

	if (use_symbol)
	{
		char label = (char)('A' + symbols.size());
		symbols.insert(std::pair<char,std::string>(label, strs.str()));
		os << label;
	}
}

/*
Example:
1296 (1123: 2 (2311), 44A, 222B, 276C, 81D; 4 (1312), 84E, 230F, 182G; 5 (1321: 1, 0, 0, 0,
0; 2 (1132), 0, 0, 0; 1, 0, 0; 0; 1), 40H, 105I; 20J; 1)

A = (2345: 0, 6 (3412), 10 (3612: 1, 0, 0, 0, 0; 1, 1, 2 (4211), 0; 1, 0, 2 (3411); 1; 1),
    2 (3611), 0; 0, 6 (3215), 8 (2461: 0, 1, 1,0, 0; 0, 0, 1, 1; 0, 1, 2 (2231); 1; 0), 4 (2611); 0,
    2 (2314), 4 (2316); 2 (2315); 0)
*/
void WriteStrategy_TextFormat(std::ostream &os, const StrategyTree &tree)
{
	// Create a symbol table.
	std::map<char,std::string> symbols;

	// Outputs the root state with symbol enabled.
	WriteState_TextFormat(os, tree, tree.root(), symbols, 2);
	os << std::endl;
	os << std::endl;

	// Outputs each symbol.
	for (auto it = symbols.begin(); it != symbols.end(); ++it)
	{
		os << it->first << " =" << it->second << std::endl;
		os << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////
// StrategyTree output (XML)

void WriteStrategy_XmlFormat(std::ostream &os, const StrategyTree &tree)
{
	Rules rules = tree.rules();

	// Output XML header.
	os << "<mastermind-strategy"
		<< " pegs=" << '"' << rules.pegs() << '"'
		<< " colors=\"" << rules.colors() << '"'
		<< " repeatable=\"" << std::boolalpha << rules.repeatable() << '"'
		<< ">" << std::endl;

#if 0
	// Output summary.
	const int max_depth = 100;
	unsigned int freq[max_depth];
	unsigned int total = tree.getDepthInfo(freq, max_depth);

	os << "<summary totalsteps=\"" << total << "\">" << std::endl;
	for (int i = 1; i <= max_depth; i++)
	{
		if (freq[i-1] > 0)
		{
			os << "  <where"
				<< " steps=\"" << i << '"'
				<< " count=\"" << freq[i-1] << '"'
				<< "/>" << std::endl;
		}
	}
	os << "</summary>" << std::endl;
#endif

	// Output the strategy.
	os << "<details>" << std::endl;

	Feedback perfect = Feedback::perfectValue(tree.rules());
	int indent = 2;
	size_t level = 0;
	auto nodes = tree.traverse(tree.root());
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		size_t depth = it.depth();

		// Close deeper branches.
		for (; level > depth; --level)
		{
			os << std::setw(indent*level) << "" << "</case>" << std::endl;
		}
		level = depth;

		if (it->response() == perfect)
		{
			os << std::setw(indent*level) << "" << "<case"
				<< " guess=" << '"' << it->guess() << '"'
				<< " response=" << '"' << it->response() << '"'
				<< "/>" << std::endl;
		}
		else
		{
			os << std::setw(indent*level) << "" << "<case"
				<< " guess=" << '"' << it->guess() << '"'
				<< " response=" << '"' << it->response() << '"'
				//<< "npos=\"" << node.npossibilities << "\" "
				//<< "next=\"" << node.suggestion 
				<< ">" << std::endl;
		}
	}

	// Write closing tags.
	os << "</details>" << std::endl;
	os << "</mastermind-strategy>" << std::endl;
}

///////////////////////////////////////////////////////////////////////////
// Strategy tree input

/*
Example:
1296 (1123: 2 (2311), 44A, 222B, 276C, 81D; 4 (1312), 84E, 230F, 182G; 5 (1321: 1, 0, 0, 0,
0; 2 (1132), 0, 0, 0; 1, 0, 0; 0; 1), 40H, 105I; 20J; 1)

A = (2345: 0, 6 (3412), 10 (3612: 1, 0, 0, 0, 0; 1, 1, 2 (4211), 0; 1, 0, 2 (3411); 1; 1),
    2 (3611), 0; 0, 6 (3215), 8 (2461: 0, 1, 1,0, 0; 0, 0, 1, 1; 0, 1, 2 (2231); 1; 0), 4 (2611); 0,
    2 (2314), 4 (2316); 2 (2315); 0)
*/

#define PARSING_ERROR(msg) err << is.tellg() << ": " << msg << std::endl
#define PARSING_WARNING(msg) err << is.tellg() << ": warning: " << msg << std::endl

/// Reads a situation from a compact format strategy.
static bool ReadSituation_TextFormat(
	std::istream &is,      // input stream
	StrategyTree &tree,    // the tree to build
	std::ostream &err,     // a stream to display error
	bool read_count,       // whether to read number of remaining possibilities
	Engine &e,             // algorithm engine
	CodewordRange secrets, // list of remaining possibilities
	std::map<char,CodewordRange> &placeholders
	)
{
	// Read number of remaining possibilities.
	if (read_count)
	{
		size_t n;
		if (!(is >> n))
		{
			PARSING_ERROR("expecting number of possibilities");
			return false;
		}
		if (n != secrets.size())
		{
			PARSING_WARNING("mismatch in number of possibilities: got "
				<< n << ", expecting " << secrets.size());
			// return false;
		}
	}

	// Read strategy for this situation.
	// Three formats are allowed following the number of possibilities:
	// - empty: this situation is trivial and not elaborated;
	// - [A-Z]: this situation is expanded later by the given identifier;
	// - (...): this situation is elaborated in-place.
	char c;
	if (!(is >> c))
		return true;

	if (c >= 'A' && c <= 'Z') // [A-Z] identifies a sub-situation
	{
		placeholders.insert(std::make_pair(c, secrets));
		return true;
	}
	else if (c == '(') // in-place strategy display
	{
		// Read the guess.
		Codeword guess;
		if (!(is >> setrules(e.rules()) >> guess))
		{
			is.clear();
			PARSING_ERROR("expecting guess");
			return false;
		}

		// The next character determines the following format:
		// ) means this guess partitions the remaining possibilities
		//   into singleton cells, and this guess is one of the
		//   remaining possibilities;
		// * means this guess partitions the remaining possibilities
		//   into singleton cells, but this guess is not one of the
		//   remaining possibilities;
		// x means this guess partitions the remaining possibilities
		//   into cells with no more than two possibilities;
		// : means none of the above trivial scenarios are true, so
		//   more details follow.
		// Note that * and x are treated as indications only and
		// ignored.
		if (!(is >> c))
		{
			PARSING_ERROR("expecting one of ')', '*', 'x', ':'");
			return false;
		}
		if (c == '*' || c == 'x') // ignore indicative flags
		{
			if (!(is >> c))
			{
				PARSING_ERROR("expecting ')' or ':'");
				return false;
			}
		}

		if (c == ')')
		{
			// Trivial situation, nothing elaborated.
			return true;
		}
		else if (c == ':')
		{
			// List of situations in the following order:
			// a04,a03,a02,a01,a00;a13,a12,a11,a10;a22,a21,a20;a30;a40

			// Let's partition the remaining secrets first.
			CodewordPartition parts = e.partition(secrets, guess);

			int p = e.rules().pegs();
			for (int nA = 0; nA <= p; ++nA)
			{
				for (int nB = (p - nA); nB >= 0; --nB)
				{
					// Read number of remaining possibilities.
					size_t n;
					if (!(is >> n))
					{
						is.clear();
						PARSING_ERROR("expecting number of possibilities");
						return false;
					}

					// Note that (p-1,1) is an impossible situation, in which
					// the number of remaining possibilities is always zero.
					// Thus this situation is omitted by Irving (1979), etc.
					// However, in Knuth(1976), this redundant situation a31=0
					// is explicitly shown; so we tweak a bit to be compatible
					// with his format.
					if (nA == p - 1 && nB == 1)
					{
						// Eat the next token if it is "0,".
						if (n == 0)
						{
							char c1;
							if (is >> c1)
							{
								if (c1 == ',')
									continue;
								else
									is.unget();
							}
						}
						--nB;
					}

					// Check sub-situation size.
					size_t j = Feedback(nA, nB).value();
					CodewordRange cell = parts[j];
					if (n != cell.size())
					{
						PARSING_WARNING("mismatch in number of possibilities: "
							<< "got " << n << ", expecting " << cell.size());
						// return false;
					}

					// Read sub-state.
					if (!ReadSituation_TextFormat(is, tree, err, false, e,
						cell, placeholders))
						return false;

					// Read separator.
					char sp_want = (nB > 0)? ',' : (nA < p)? ';' : ')';
					char sp_read;
					if (!(is >> sp_read) || (sp_read != sp_want))
					{
						PARSING_ERROR("expecting '" << sp_want << "', got '"
							<< sp_read << "'");
						return false;
					}
				}
			}
			return true;
		}
		else // bad syntax
		{
			PARSING_ERROR("expecting ')' or ':', got '" << c << "'");
			return false;
		}
	}
	else // trivial, non-elaborated situation
	{
		is.unget();
		return true;
	}
}

/*
Example:
1296 (1123: 2 (2311), 44A, 222B, 276C, 81D; 4 (1312), 84E, 230F, 182G; 5 (1321: 1, 0, 0, 0,
0; 2 (1132), 0, 0, 0; 1, 0, 0; 0; 1), 40H, 105I; 20J; 1)

A = (2345: 0, 6 (3412), 10 (3612: 1, 0, 0, 0, 0; 1, 1, 2 (4211), 0; 1, 0, 2 (3411); 1; 1),
    2 (3611), 0; 0, 6 (3215), 8 (2461: 0, 1, 1,0, 0; 0, 0, 1, 1; 0, 1, 2 (2231); 1; 0), 4 (2611); 0,
    2 (2314), 4 (2316); 2 (2315); 0)
*/
static bool ReadStrategy_TextFormat(std::istream &is, StrategyTree &tree, std::ostream &err)
{
	// Create a list of all the secrets.
	Engine e(tree.rules());
	CodewordList secrets(e.generateCodewords());
	CodewordRange all(secrets);

	// Create a map of situation placeholders (symbol table).
	std::map<char,CodewordRange> placeholders;
	std::map<char,bool> defined;

	// Read root situation.
	if (!ReadSituation_TextFormat(is, tree, err, true, e, all, placeholders))
		return false;

	// Read each sub-situation identified by [A-Z]. For example,
	// A = (2345: 0, ....
	char id;
	while (is >> id)
	{
		char c;
		if (!(is >> c) || (c != '='))
		{
			PARSING_ERROR("expecting '='");
			return false;
		}
		if (placeholders.find(id) == placeholders.end())
		{
			PARSING_ERROR("undeclared situation identifier '" << id << "'");
			return false;
		}
		if (defined[id])
		{
			PARSING_ERROR("situation " << id << " is already defined");
			return false;
		}
		CodewordRange cell = placeholders.at(id);
		if (!ReadSituation_TextFormat(is, tree, err, false, e, cell, placeholders))
			return false;
		defined[id] = true;
	}

	// Check for declared but undefined sub-situations.
	if (defined.size() < placeholders.size())
	{
		is.clear();
		PARSING_ERROR("some situations are not defined");
		return false;
	}

	return true;
}

std::istream& operator >> (std::istream &is, StrategyTree &tree)
{
	if (!ReadStrategy_TextFormat(is, tree, std::cerr))
	{
		// Set fail bit.
	}
	return is;
}

} // namespace Mastermind

#if 0


void StrategyTreeNode::WriteToFile(std::ostream &os, FileFormat format, int indent) const
{
	assert(indent >= 0 && indent < 200);

	for (int i = 0; i < m_childcount; i++)
	{
		int j = m_childindex[i];
		const StrategyTreeNode *child = m_children[j];
		if (format == XmlFormat)
		{
			if (child == Done())
			{
				os << std::setw(indent) << "" << "<case feedback=\""
					<< Feedback(i) << "\">" << std::endl;
			}
			else
			{
				os << std::setw(indent) << "" << "<case "
					<< "feedback=\"" << Feedback(i) << "\" "
					<< "npos=\"" << child->State.NPossibilities << "\" "
					<< "ncand=\"" << child->State.NCandidates << "\" "
					<< "guess=\"" << child->State.Guess << "\">"
					<< std::endl;
			}
		}
		else
		{
			os << std::setw(indent) << ""
				<< State.Guess << ":" << Feedback(i) << std::endl;
		}
		if (child != Done())
		{
			child->WriteToFile(os, format, indent + 2);
		}
		if (format == XmlFormat)
		{
			os << std::setw(indent) << "" << "</case>" << std::endl;
		}
	}
}

void StrategyTreeNode::WriteToFile(
	std::ostream &os, FileFormat format, const Rules &rules) const
{
	if (format == XmlFormat)
	{
		os << "<mmstrat"
			<< " npegs=" << '"' << rules.pegs() << '"'
			<< " ncolors=\"" << rules.colors() << '"'
			<< " allowrepetition=\"" << std::boolalpha << rules.repeatable() << '"'
			<< ">" << std::endl;

		const int max_depth = 100;
		int freq[max_depth+1];
		int total = GetDepthInfo(freq, max_depth);

		os << "<summary totalsteps=\"" << total << "\">" << std::endl;
		for (int i = 0; i <= max_depth; i++)
		{
			if (freq[i] > 0)
			{
				os << "  <where"
					<< " steps=\"" << i << '"'
					<< " count=\"" << freq[i] << '"'
					<< "/>" << std::endl;
			}
		}
		os << "</summary>\n";
		os << "<details>\n";
	}
	WriteToFile(os, format, 0);
	if (format == XmlFormat)
	{
		os << "</details>" << std::endl;
		os << "</mmstrat>" << std::endl;
	}
}



#endif
