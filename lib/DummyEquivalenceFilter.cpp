#include "Equivalence.hpp"

namespace Mastermind {

/// Represents a dummy equivalence filter that doesn't filter anything.
class DummyEquivalenceFilter : public EquivalenceFilter
{
public:

	virtual EquivalenceFilter* clone() const 
	{
		return new DummyEquivalenceFilter();
	}

	virtual CodewordList get_canonical_guesses(
		CodewordConstRange candidates
		) const 
	{
		return CodewordList(candidates.begin(), candidates.end());
	}

	virtual void add_constraint(
		const Codeword & /* guess */,
		Feedback /* response */, 
		CodewordConstRange /* remaining */)
	{
	}
};

EquivalenceFilter* CreateDummyEquivalenceFilter(const Engine *)
{
	return new DummyEquivalenceFilter();
}

#if 0
REGISTER_ROUTINE(CreateEquivalenceFilterRoutine,
				 "Dummy",
				 CreateDummyEquivalenceFilter)
#endif

} // namespace Mastermind
