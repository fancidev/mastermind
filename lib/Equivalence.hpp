#ifndef MASTERMIND_EQUIVALENCE_HPP
#define MASTERMIND_EQUIVALENCE_HPP

#include <memory>
#include "Engine.hpp"

namespace Mastermind {

/// Defines an interface for an equivalence filter that filters canonical 
/// guesses from a set of candidate codewords. 
/// @ingroup equiv
struct EquivalenceFilter
{
	/// Destroys the equivalence filter.
	virtual ~EquivalenceFilter() { }

	/// Allocates and initializes an identical filter to this one.
	/// The allocated memory must be freed with @c delete.
	virtual EquivalenceFilter* clone() const = 0;

	/// Returns a list of canonical guesses from a set of candidates.
	virtual CodewordList get_canonical_guesses(
		CodewordConstRange candidates
		) const = 0;

	/// Adds a constraint to the current state.
	virtual void add_constraint(
		const Codeword &guess,
		Feedback response, 
		CodewordConstRange remaining
		) = 0;
};

/// Typedef of pointer to function that creates an equivalence filter.
/// @ingroup equiv
// typedef EquivalenceFilter* (*CreateEquivalenceFilterRoutine)(const Engine *e);

extern EquivalenceFilter* CreateDummyEquivalenceFilter(const Engine *e);
extern EquivalenceFilter* CreateColorEquivalenceFilter(const Engine *e);
extern EquivalenceFilter* CreateConstraintEquivalenceFilter(const Engine *e);

/// Composite equivalence filter which chains two underlying filters.
/// @ingroup equiv
class CompositeEquivalenceFilter : public EquivalenceFilter
{
	std::unique_ptr<EquivalenceFilter> _filter1, _filter2;

public:

	/// Constructs a composite filter that chains the supplied individual 
	/// filters. The supplied filters are cloned so that the supplied 
	/// pointers can be disposed by the caller.
	CompositeEquivalenceFilter(
		const EquivalenceFilter *filter1, 
		const EquivalenceFilter *filter2)
		: _filter1(filter1->clone()), _filter2(filter2->clone())
	{
	}

#if 0
	// Constructs a composite filter that chains two individual filters
	/// of the given names in the routine registry of
	/// @c CreateEquivalenceFilterRoutine.
	CompositeEquivalenceFilter(const Engine *e, const char *name1, const char *name2)
		: _filter1(RoutineRegistry<CreateEquivalenceFilterRoutine>::get(name1)(e)),
		 _filter2(RoutineRegistry<CreateEquivalenceFilterRoutine>::get(name2)(e))
	{
	}
#endif

	virtual EquivalenceFilter* clone() const 
	{
		return new CompositeEquivalenceFilter(_filter1.get(), _filter2.get());
	}

	virtual CodewordList get_canonical_guesses(
		CodewordConstRange candidates
		) const 
	{
		CodewordList temp = _filter1->get_canonical_guesses(candidates);
		return _filter2->get_canonical_guesses(temp);
	}

	virtual void add_constraint(
		const Codeword & guess,
		Feedback response, 
		CodewordConstRange remaining)
	{
		_filter1->add_constraint(guess, response, remaining);
		_filter2->add_constraint(guess, response, remaining);
	}

	/// Returns the first filter.
	const EquivalenceFilter* first() const { return _filter1.get(); }

	/// Returns the second filter.
	const EquivalenceFilter* second() const { return _filter2.get(); }
};

} // namespace Mastermind

#endif // MASTERMIND_EQUIVALENCE_HPP
