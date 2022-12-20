/// @defgroup HRTimer High-Resolution Timer
/// @ingroup util

#ifndef UTILITIES_HR_TIMER_HPP
#define UTILITIES_HR_TIMER_HPP

#if _OPENMP
#include <omp.h>
#else
#include <ctime>
#endif

namespace util
{

/// High-resolution timer used to time program execution.
/// @ingroup HRTimer
class hr_timer
{
	double _start;

public:

	/// Creates a timer.
	hr_timer() : _start(0.0) { }

	/// Starts timing.
	void start()
	{
#if _OPENMP
		_start = omp_get_wtime();
#else
		_start = (double)time(NULL);
#endif
	}

	/// Stops timing and returns the number of seconds elapsed
	/// since the last call to <code>start()</code>.
	double stop()
	{
#if _OPENMP
		return omp_get_wtime() - _start;
#else
		return (double)time(NULL) - _start;
#endif
	}
};

} // namespace util

#endif // UTILITIES_HR_TIMER_HPP
