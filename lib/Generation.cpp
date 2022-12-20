#include "Algorithm.hpp"

namespace Mastermind {

static void generate_recursion(
	int npegs, int ncolors, int max_repeat,
	int peg, const Codeword &_partial, Codeword* &output)
{
	Codeword partial(_partial);
	for (int k = 0; k < ncolors; ++k)
	{
		if (partial.count(k) < max_repeat)
		{
			partial.set(peg, k);
			if (peg == npegs - 1)
				*(output++) = partial;
			else
				generate_recursion(npegs, ncolors, max_repeat, peg+1, partial, output);
		}
	}
}

void GenerateCodewords(const Rules &rules, Codeword *results)
{
	int pegs = rules.pegs();
	int colors = rules.colors();
	generate_recursion(pegs, colors, rules.repeatable()? pegs : 1,
		0, Codeword(), results);
}

} // namespace Mastermind
