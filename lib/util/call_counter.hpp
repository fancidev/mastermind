/**
 * @defgroup CallCounter Call Counter
 * Utilities to collect function call statistics.
 *
 * @ingroup util
 */

#ifndef UTILITIES_CALL_COUNTER_HPP
#define UTILITIES_CALL_COUNTER_HPP

#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include "intrinsic.hpp"

namespace util {

/**
 * Represents a counter that collects function call statistics.
 * It supports to up to <code>2^32-2</code> operations per call for 32-bit
 * and <code>2^64-2</code> operations per call for 64-bit. Under either
 * architecture, the total number of calls (or operations) supported is 
 * up to <code>2^64-1</code>.
 * @ingroup CallCounter
 */
class call_counter
{
	// Name of the routine being profiled.
	std::string _name; 

	// Classifies the call statistics into groups, grouped by the 
	// number of operations in a call as follows:
	// [0,1), [1,3), [3,7), ..., [2^31-1,2^32-1)
	struct group_t
	{
		unsigned long long calls; // number of calls
		unsigned long long ops;   // number of operations
		group_t() : calls(0), ops(0) { }
	};

	// N = number of groups, determined by the size of size_t.
	static const int N = sizeof(size_t)*8;
	group_t _stat[N];

	typedef std::map<std::string, call_counter> registry_type;

	/// Returns a map of registered call counters.
	static registry_type& get_registry()
	{
		static registry_type ccs;
		return ccs;
	}

#if 0
	// Use static variable which probably saves a branching instruction.
	// However, the variable must be explicitly defined in some CPP file.
	static bool _enabled;
	static bool &get_enabled() { return _enabled; }
#else
	// Use static local variable. This doesn't require explicit definition
	// in another CPP file.
	static bool &get_enabled()
	{
		static bool _enabled = true;
		return _enabled;
	}
#endif

public:

	/// Constructs a call counter of the given display name.
	call_counter(const std::string &name) : _name(name) { }

	/// Returns the name of the routine being profiled.
	std::string name() const { return _name; }

	/// Returns the total number of calls recorded by this counter.
	unsigned long long total_calls() const 
	{
		unsigned long long total = 0;
		for (int i = 0; i < N; ++i)
			total += _stat[i].calls;
		return total;
	}

	/// Returns the total number of operations recorded by this counter.
	unsigned long long total_ops() const 
	{
		unsigned long long total = 0;
		for (int i = 0; i < N; ++i)
			total += _stat[i].ops;
		return total;
	}

	/// Records one function call with a given number of operations.
	/// @param ops Number of operations in this call.
	void add_call(size_t ops)
	{
		assert(ops != (size_t)-1);

		// The group is the location of the most significant 1.
		unsigned int pos = util::intrinsic::bit_scan_reverse(ops + 1);
		_stat[pos].calls++;
		_stat[pos].ops += ops;
	}

	/// Outputs the call counter statistics to a stream.
	/// @param os The stream.
	/// @param cc The call counter.
	/// @returns The stream.
	friend std::ostream& operator << (std::ostream &os, const call_counter &cc)
	{
		unsigned long long ncalls = cc.total_calls();
		unsigned long long nops = cc.total_ops();

		os << "==== Call Statistics for " << cc.name() << "() ====" << std::endl;
		os << "Total # of calls : " << ncalls << std::endl;
		os << "Total # of ops   : " << nops << std::endl;

		if (ncalls > 0)
		{
			double ops_per_call = (double)nops / ncalls;
			os << "Avg ops per call : " << ops_per_call << std::endl;
			os << "#Ops             %Calls    %Ops     #Calls         #Ops" << std::endl;
			for (int k = call_counter::N - 1; k >= 0; k--)
			{
				if (cc._stat[k].calls > 0)
				{
					os << "[" << std::setw(6) << (1<<k)-1 << " -"
						<< std::setw(6) << (1<<(k+1))-2 << "] "
						<< std::setw(6) << std::setprecision(2) << std::fixed
						<< (double)cc._stat[k].calls / ncalls * 100.0 << "  "
						<< std::setw(6) << std::setprecision(2) << std::fixed
						<< (double)cc._stat[k].ops / nops * 100.0 << " "
						<< std::setw(10) << cc._stat[k].calls << " "
						<< std::setw(12) << cc._stat[k].ops << std::endl;
				}
			}
		}
		return os;
	}

	/// Returns a reference to the global call counter of the given name.
	/// If no call counter of this name exists, a new one with the name
	/// is created and returned.
	static call_counter& get(const std::string &name)
	{
		registry_type &ccs = get_registry();
		auto it = ccs.find(name);
		if (it == ccs.end())
		{
			it = ccs.insert(std::make_pair(name, call_counter(name))).first;
		}
		return it->second;
	}

	/// Enables or disables call counter statistics globally.
	static void enable(bool flag)
	{
		get_enabled() = flag;
	}

	/// Checks whether call counter is enabled at a global level.
	static bool enabled()
	{
		return get_enabled();
	}

	/// Returns a map of all registered call counters.
	static const registry_type& registry()
	{
		return get_registry();
	}
};

} // namespace util

#ifndef ENABLE_CALL_COUNTER
/** 
 * Define <code>ENABLE_CALL_COUNTER = 0</code> explicitly to disable the 
 * <code>UPDATE_CALL_COUNTER</code> macro at compile time. This completely
 * removes the run-time overhead of updating call counters. 
 *
 * Alternatively, <code>UPDATE_CALL_COUNTER</code> can be enabled or disabled
 * at run-time via <code>call_counter::enable()</code>. However, some small 
 * overhead still remains in this case.
 *
 * If <code>ENABLE_CALL_COUNTER</code> is not explicitly defined and 
 * <code>call_counter::enable()</code> is not called, call counters are 
 * enabled by default.
 *
 * @remarks This macro MUST be defined the same in all translation units.
 *      Otherwise the compiler may cache part of the generated code and the
 *      final result is unpredictable.
 *
 * @ingroup CallCounter
 */
#define ENABLE_CALL_COUNTER 1
#endif

/// Wraps a block of statements in as a single statement.
#ifdef _MSC_VER
#define WRAP_STATEMENTS(stmt) \
	__pragma(warning(push)) \
	__pragma(warning(disable:4127)) /* conditional expression is constant */ \
	do { stmt ; } while (0) \
	__pragma(warning(pop))
#else
#define WRAP_STATEMENTS(stmt) do { stmt ; } while (0) 
#endif

/// Registers a global call counter of the given name (if not already
/// registered) and updates its call statistics.
/// @ingroup CallCounter
#if ENABLE_CALL_COUNTER
#define UPDATE_CALL_COUNTER(id,nops) WRAP_STATEMENTS( \
	if (util::call_counter::enabled()) { \
		static util::call_counter& _call_counter = util::call_counter::get(id); \
		_call_counter.add_call(nops); \
	} )
#else
#define UPDATE_CALL_COUNTER(id,nops) do {} while (0)
#endif

#endif // UTILITIES_CALL_COUNTER_HPP
