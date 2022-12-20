#include "Engine.hpp"

namespace Mastermind {

// @todo We could move the implementation of compare() to a header file
// and then implement a custom updater to do the filtering.
CodewordList Engine::filterByFeedback(
	const CodewordList &list,
	const Codeword &guess,
	const Feedback &feedback) const
{
	FeedbackList fblist;
	FeedbackFrequencyTable freq = compare(guess, list, fblist);

	// Count feedbacks equal to feedback.
	size_t count = freq[feedback.value()];

	// Copy elements whose feedback are equal to fb.
	CodewordList result(count);
	size_t j = 0;
	for (size_t i = 0; i < fblist.size(); i++)
	{
		if (fblist[i] == feedback)
			result[j++] = list[i];
	}
	return result;
}

CodewordPartition Engine::partition(
	CodewordRange codewords,
	const Codeword &guess) const
{
	// If there's no element in the list, do nothing.
	if (codewords.empty())
		return CodewordPartition();

	// Compare the guess to each codeword in the list.
	FeedbackList fbl;
	FeedbackFrequencyTable freq = compare(guess, codewords, fbl);

	// Build a table to store the range of each partition.
	struct partition_location
	{
		//size_t begin; // begin of the partition
		size_t end;     // end of the partition
		size_t current; // next location to insert
	} part[256+1];

	size_t i = 0;
	part[0].current = 0;
	for (size_t k = 0; k < freq.size(); k++)
	{
		i += freq[k];
		part[k].end = i;
		part[k+1].current = i;
	}
	part[freq.size()].end = std::numeric_limits<size_t>::max();

	// Find the first non-empty partition.
	int k = 0; // current partition
	while (freq[k] == 0)
		++k;

	// Perform a in-place partitioning.
	CodewordIterator first = codewords.begin();
	size_t count = codewords.size();
	for (size_t i = 0; i < count; )
	{
		int fbv = fbl[i].value();
		//std::cout << "Feedback[" << i << "] = " << fbl[i] << std::endl;
		if (fbv == k)
		{
			// Codeword[i] is in the correct partition.
			// Advance the current partition pointer.
			// If it's reached the end, move to the next partition.
			if (++part[k].current >= part[k].end)
			{
				for (++k; part[k].current >= part[k].end; ++k);
			}
			i = part[k].current;
		}
		else
		{
			// Codeword[i] is NOT in the correct partition.
			// Swap it into the correct partition, and increment
			// the pointer of that partition.
			size_t j = part[fbv].current++;
			std::swap(first[i], first[j]);
			std::swap(fbl[i], fbl[j]);
		}
	}
	//return freq;
	return CodewordPartition(codewords, freq);
}

} // namespace Mastermind
