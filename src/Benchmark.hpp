#ifndef MASTERMIND_BENCHMARK_HPP
#define MASTERMIND_BENCHMARK_HPP

#include <iostream>
#include <iomanip>
#include "Engine.hpp"
#include "util/hr_timer.hpp"

namespace Mastermind {

/// Dummy test driver that does nothing in the test and always returns success.
/// @ingroup prog
template <class Routine>
class TestDriver
{
	const Engine *e;
	Routine f;

public:

	/// Constructs a dummy test driver.
	TestDriver(const Engine *engine, Routine func) : e(engine), f(func) { }

	/// Runs the test.
	bool operator () () const { return true; }
};

/// Compares the running time of two routines.
/// @ingroup prog
template <class Routine>
bool compareRoutines(
	const Engine *e,
	const char *routine1,
	const char *routine2,
	long times)
{
	Routine func1 = RoutineRegistry<Routine>::get(routine1);
	Routine func2 = RoutineRegistry<Routine>::get(routine2);

	TestDriver<Routine> drv1(e, func1);
	TestDriver<Routine> drv2(e, func2);

	// Verify computation results.
	drv1();
	drv2();
	if (!(drv1 == drv2))
		return false;

	// Display results in debug mode.
	if (times == 0)
	{
		std::cout << "Result 1: " << std::endl << drv1 << std::endl;
		std::cout << "Result 2: " << std::endl << drv2 << std::endl;
	}

	// Time it.
	util::hr_timer timer;
	double t1 = 0, t2 = 0;

	for (int pass = 0; pass < 10; pass++)
	{
		timer.start();
		for (int k = 0; k < times / 10; k++)
			drv1();
		t1 += timer.stop();

		timer.start();
		for (int k = 0; k < times / 10; k++)
			drv2();
		t2 += timer.stop();
	}

	std::cout << "Algorithm 1: " << std::setw(6) << std::setprecision(3) 
		<< t1 << std::endl;
	std::cout << "Algorithm 2: " << std::setw(6) << std::setprecision(3) 
		<< t2 << std::endl;
	std::cout << "Throughput Ratio: " << std::setw(5) << std::setprecision(2) 
		<< (t1/t2) << std::endl;
	return true;
}

#if 0
// Codeword generation benchmark.
// Test: Generate all codewords of 4 pegs, 10 colors, and no repeats.
//       Total 5040 items in each run.
// Results: (100,000 runs, Release mode)
// LexOrder: 4.43 s
// CombPerm: 8.54 s [legacy]
// CombPermParallel:  0.96 s [ASM][legacy]
// CombPermParallel2: 0.68 s [ASM][legacy]
template <> class TestDriver<GenerationRoutine>
{
	Engine &e;
	GenerationRoutine f;
	size_t count;
	CodewordList list;

public:
	TestDriver(Engine &engine, GenerationRoutine func)
		: e(engine), f(func), count(f(e.rules(), 0)), list(count) { }

	void operator()() { 	f(e.rules(), list.data()); }

	bool operator == (const TestDriver &r) const
	{
		return list == r.list;
	}

	friend std::ostream& operator << (std::ostream &os, const TestDriver &r)
	{
		// Display first 10 items.
		os << "First 10 of " << r.list.size() << " items:" << std::endl;
		for (size_t i = 0; i < 10 && i < r.list.size(); ++i)
			os << r.list[i] << std::endl;
		return os;
	}
};
#endif

#if 1
/**
 * Test driver that executes a codeword comparison routine.
 *
 * Benchmark results:
 * <pre>
 * Test:     Compare a given codeword to 5040 non-repeatable codewords, and
 *           update a frequency table of the responses.
 * Results:  (100,000 runs, x64, VC++ 2010)
 * generic:  1.68 s (feedback) / 2.46 s (freq)
 * norepeat: 0.62 s (feedback) / 1.18 s (freq)
 * </pre>
 *
 * @ingroup prog
 */
template <> class TestDriver<ComparisonRoutine>
{
	const Engine *e;
	ComparisonRoutine f;
	CodewordList codewords;
	size_t count;
	FeedbackFrequencyTable freq;
	Codeword secret;

public:

	/// Constructs the test driver.
	TestDriver(const Engine *engine, ComparisonRoutine func)
		: e(engine), f(func), codewords(e->generateCodewords()),
		count(codewords.size()), secret(codewords[count/2])
		{ }

	/// Runs the test.
	void operator()()
	{
		freq.resize(Feedback::size(e->rules()));
		f(secret, codewords.data(), codewords.size(), 0, freq.data());
	}

	/// Compares the results of two runs.
	bool operator == (const TestDriver &r) const
	{
		if (freq.size() != r.freq.size())
		{
			std::cout << "**** ERROR: Different sizes." << std::endl;
			return false;
		}
		for (size_t i = 0; i < freq.size(); i++)
		{
			if (freq[i] != r.freq[i])
			{
				std::cout << "**** ERROR: Inconsistent frequency for ["
					<< Feedback(i) << "]: " << freq[i] << " v " << r.freq[i]
					<< std::endl;
				return false;
			}
		}
		return true;
	}

	/// Outputs the result.
	friend std::ostream& operator << (std::ostream &os, const TestDriver &r)
	{
		return os << r.freq;
	}
};
#endif

#if 0
// Color-mask scanning benchmark.
// Test: Scan 5040 codewords for 100,000 times.
//
// **** Old Results ****
// These results are for "short" codewords.
//
// ScanDigitMask_v1 (C):              5.35 s
// ScanDigitMask_v2 (16-bit ASM):     2.08 s
// ScanDigitMask_v3 (v2 improved):    1.43 s
// ScanDigitMask_v4 (v3 improved):    1.12 s
// ScanDigitMask_v5 (32-bit ASM):     2.09 s
// ScanDigitMask_v6 (v5 improved):    1.10 s
// ScanDigitMask_v7 (v6 generalized): 1.10 s
//
// Observations:
//   - ASM with parallel execution and loop unrolling performs the best.
//   - There is little performance difference between 16-bit ASM and 32-bit ASM.
//   - Loop unrolling has limited effect. Seems 1.10s is the lower bound
//     the current algorithm can improve to.
//
// **** New Results ****
// These results are for "long" codewords.
// ScanDigitMask_v1 (SSE2): 0.40 s
template <> class TestDriver<MaskRoutine>
{
	Engine &e;
	MaskRoutine f;
	CodewordList list;
	unsigned short mask;

public:
	TestDriver(Engine &env, MaskRoutine func)
		: e(env), f(func), list(e.generateCodewords()), mask(0) { }

	void operator()()
	{
		mask = list.empty()? 0 : f(&list.front(), &list.back() + 1);
	}

	bool operator == (const TestDriver &r) const
	{
		if (mask != r.mask)
		{
			std::cout << "**** Inconsistent color mask ****" << std::endl;
			return false;
		}
		return true;
	}

	friend std::ostream& operator << (std::ostream &os, const TestDriver &r)
	{
		os << "Present digits:";
		unsigned short t = r.mask;
		for (int i = 0; i < 16; i++)
		{
			if (t & 1)
				os << " " << i;
			t >>= 1;
		}
		return os << std::endl;
	}
};
#endif

#if 0
/// Compares frequency counting algorithms.
///
/// Test: Compute the frequency table of 5040 feedbacks.
/// Run the test for 500,000 times.
///
/// Results:
/// <pre>
/// count_freq_v1 (plain C):              5.55 s
/// count_freq_v2 (ASM, 2-parallel OOE):  3.85 s
/// count_freq_v3 (C, 4-parallel OOE):    5.51 s
/// count_freq_v4 (ASM, loop-unroll):     4.54 s
/// count_freq_v5 (ASM, 4-parallel OOE):  2.99 s
/// count_freq_v6 (improved ASM from v5): 2.88 s
/// </pre>
/// Conclusion:
/// Due to the memory-intensive nature of the algorithm, the performance
/// cannot be improved much. count_freq_v6() is the chosen implementation.
int TestFrequencyCounting(const Rules &rules, long times)
{
	CodewordList list = generateCodewords(rules);
	FeedbackList fblist = compare(rules, list[0], list.cbegin(), list.cend());
	int count = fblist.size();
	const unsigned char *fbl = (const unsigned char *)fblist.data();
	const unsigned char maxfb = 63;
	unsigned int freq[(int)maxfb+1];

	FREQUENCY_COUNTING_ROUTINE *func1 = CountFrequenciesImpl->GetRoutine("c");
	FREQUENCY_COUNTING_ROUTINE *func2 = CountFrequenciesImpl->GetRoutine("c_luf4");

	if (times == 0) {
		int total = 0;
		//count = 11;
		func2(fbl, count, freq, maxfb);
		for (int i = 0; i <= maxfb; i++) {
			if (freq[i] > 0)
			{
				std::cout << Feedback(i) << " = " << freq[i] << std::endl;
				total += freq[i];
			}
		}
		std::cout << "Expected count: " << count << std::endl;
		std::cout << "Actual total:   " << total << std::endl;
		system("PAUSE");
		return 0;
	}

	HRTimer timer;
	double t1, t2;
	t1 = t2 = 0;

	for (int pass = 0; pass < 10; pass++) {
		timer.start();
		for (int j = 0; j < times / 10; j++) {
			func1(fbl, count, freq, maxfb);
		}
		t1 += timer.stop();

		timer.start();
		for (int j = 0; j < times / 10; j++) {
			func2(fbl, count, freq, maxfb);
		}
		t2 += timer.stop();
	}

	printf("Algorithm 1: %6.3f\n", t1);
	printf("Algorithm 2: %6.3f\n", t2);
	printf("Improvement: %5.1f%%\n", (t1/t2-1)*100);

	system("PAUSE");
	return 0;
}
#endif

} // namespace Mastermind

#endif // MASTERMIND_BENCHMARK_HPP
