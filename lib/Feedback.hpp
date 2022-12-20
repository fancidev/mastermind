#ifndef MASTERMIND_FEEDBACK_HPP
#define MASTERMIND_FEEDBACK_HPP

#include <cassert>
#include <iostream>
#include <utility>
#include <algorithm>
#include <string>

#include "Rules.hpp"

namespace Mastermind {

/**
 * Represents the feedback from comparing two codewords.
 *
 * For a feedback of the form <code>xAyB</code>, this class stores
 * the pair <code>(x,y)</code> and name it as <code>(nA,nB)</code>.
 * To improve memory access locality, the feedback is internally
 * represented by its ordinal position in an order described below.
 *
 * The set of all distinct feedback values can be arranged in the
 * form of a triangle as follows:
 * <pre>
 * 0A0B  1A0B  2A0B  3A0B  4A0B
 * 0A1B  1A1B  2A1B  3A1B*
 * 0A2B  1A2B  2A2B
 * 0A3B  1A3B
 * 0A4B
 * </pre>
 * The arrangement looks like a Pythagoras triangle, in that the
 * feedbacks in any diagonal have the same sum <code>nA+nB</code>.
 * The feedbacks are then ordered like this:
 *   0A0B -> 0A1B -> 1A0B -> 0A2B -> 1A1B -> 2A0B -> ...
 * That is, they are first sorted by <code>nA+nB</code>, and then
 * by <code>nA</code>. Such ordering can be conveniently interpreted
 * as increasing level of matching. For example, 0A0B is no match,
 * and 4A0B is a perfect match.
 *
 * A game with <code>p</code> pegs can have a range of feedbacks
 * up to the <code>p</code>-th diagonal (starting from zero).
 * Thus the total number of valid feedbacks in a <code>p</code>-peg
 * game is <code>1+2+...+(p+1)</code>, which is equal to
 * <code>(p+1)*(p+2)/2</code>. 
 *
 * Note, however, that the last feedback in the row marked with [*]
 * is a practically-impossible feedback for a <code>p</code>-peg game. 
 * For example, in a game with 4 pegs, an outcome can never be 3A1B. 
 * Therefore, the total number of practically possible feedbacks in a
 * <code>p</code>-peg game is <code>(p+1)*(p+2)-1</code>, which is 
 * equal to <code>p*(p+3)/2</code>.
 *
 * To convert a feedback <code>(nA,nB)</code> to its ordinal position
 * in the above arrangement, first compute <code>nAB=nA+nB</code> and
 * then use the following formula
 *   <code>pos = nAB*(nAB+1)/2+nA</code>.
 * The largest feedback value for a <code>p</code>-peg game is
 * <code>p*(p+3)/2</code>. This means 7-bits storage can represent
 * up to 14 pegs, which corresponds to a maximum ordinal index of 119.
 * The highest bit is reserved to indicate an error condition.
 *
 * To convert a feedback from its ordinal position, there is no
 * simple formula. Therefore we build a lookup table for this purpose.
 *
 * @ingroup Feedback
 */
class Feedback
{
public:

	/// Type of the internal representation of the feedback.
	typedef signed char value_type;

	/// Maximum number of distinct feedback outcomes.
	static const int MaxOutcomes = (MM_MAX_PEGS+1)*(MM_MAX_PEGS+2)/2;

private:

	/// Ordinal position of the feedback. The special value @c -1
	/// is reserved to represent an empty feedback.
	value_type _value;

	/// Lookup table that maps an ordinal position to a pair of
	/// <code>(nA,nB)</code>.
	struct outcome_table
	{
		std::pair<int,int> table[MaxOutcomes];
		
		outcome_table()
		{
			for (int nAB = 0; nAB <= MM_MAX_PEGS; ++nAB)
			{
				for (int nA = 0; nA <= nAB; ++nA)
				{
					int nB = nAB - nA;
					int k = (nAB+1)*nAB/2+nA;
					assert(k < MaxOutcomes);
					table[k] = std::pair<int,int>(nA, nB);
				}
			}
		}

		static std::pair<int,int> lookup(value_type value)
		{
			static const outcome_table ot;
			return (value >= 0 && value < MaxOutcomes)?
				ot.table[(size_t)value] : std::pair<int,int>(-1,-1);
		}
	};

public:

	/// Creates an empty feedback.
	Feedback() : _value(-1) { }

	/// Creates a feedback from its ordinal position.
	explicit Feedback(size_t index) : _value((value_type)index) { }

	/// Creates a feedback with the given <code>nA</code> and <code>nB</code>.
	/// If the arguments are not valid, an empty feedback is created.
	explicit Feedback(int nA, int nB)
	{
		int nAB = nA + nB;
		if (nA >= 0 && nB >= 0 && nAB <= MM_MAX_PEGS)
			_value = (value_type)((nAB+1)*nAB/2+nA);
		else
			_value = -1;
	}

	/// Creates a feedback from a string of the form "1A2B". If the
	/// input string is invalid, the feedback is set to an empty value.
	explicit Feedback(const char *s) : _value(-1)
	{
		if (s &&
			(s[0] >= '0' && s[0] <= '9') &&
			(s[1] == 'A' || s[1] == 'a') &&
			(s[2] >= '0' && s[2] <= '9') &&
			(s[3] == 'B' || s[3] == 'b') &&
			(s[4] == '\0'))
		{
			*this = Feedback(s[0]-'0', s[2]-'0');
		}
	}

	/// Returns the internal representation of the feedback.
	value_type value() const { return _value; }

	/// Tests whether the feedback is empty.
	bool empty() const { return _value < 0; }

	/// Tests whether the feedback is empty.
	bool operator ! () const { return empty(); }

	/// Tests whether the feedback is non-empty.
	operator void* () const { return empty() ? 0 : (void*)this; }

	/// Returns <code>nA</code>, the number of correct colors in the 
	/// correct pegs. If this feedback is empty, returns -1.
	int nA() const { return outcome_table::lookup(_value).first; }

	/// Returns <code>nB</code>, the number of correct colors in the
	/// wrong pegs. If this feedback is empty, returns -1.
	int nB() const { return outcome_table::lookup(_value).second; }

	/// Tests whether the feedback conforms to the given set of rules.
	bool conforming(const Rules &rules) const
	{
		if (!rules)
			return false;
		
		if (empty())
			return false;

		int a = nA(), b = nB(), p = rules.pegs();
		return (a >= 0 && b >= 0 && a + b <= p);
	}

	/**
	 * Compact format of a feedback.
	 * In this format, bits 0-3 stores @c nB, and bits 4-7 stores
	 * @c nA. To convert a feedback into the compact format, use
	 * the formula <code>x = (nA << 4) | nB</code>.
	 * To restore a feedback from compact format, use the formula
	 * <code>nA = x >> 4</code> and <code>nB = x & 0x0F</code>.
	 */
	typedef unsigned char compact_type;

	/// Converts the feedback into compact form.
	compact_type pack() const
	{
		std::pair<int,int> ab = outcome_table::lookup(_value);
		return (compact_type)((ab.first << 4) | ab.second);
	}

	/// Restores a feedback from compact form.
	/// @todo check malformed input
	static Feedback unpack(compact_type ab)
	{
		int nA = (ab >> 4), nB = (ab & 0x0F);
		return Feedback(nA, nB);
	}

	/// Returns the feedback for a perfect match under a given set
	/// of rules. This is the largest valid feedback value for this
	/// set of rules.
	static Feedback perfectValue(const Rules &rules)
	{
		return Feedback(rules.pegs(), 0);
	}

	/// Returns the size of the set of distinct feedback values under
	/// a given set of rules. The practically impossible feedback
	/// <code>(p-1,1)</code> is included.
	static size_t size(const Rules &rules)
	{
		return (size_t)perfectValue(rules)._value + 1;
	}
};

/// Tests whether two feedbacks are equal.
/// @ingroup Feedback
inline bool operator == (const Feedback &a, const Feedback &b)
{
	return a.value() == b.value();
}

/// Tests whether two feedbacks are not equal.
/// @ingroup Feedback
inline bool operator != (const Feedback &a, const Feedback &b)
{
	return a.value() != b.value();
}

/// Outputs the feedback to a stream in the form "1A2B".
/// @ingroup Feedback
inline std::ostream& operator << (std::ostream &os, const Feedback &feedback)
{
	char s[5] = "-A-B";
	if (feedback)
	{
		s[0] = (char)('0' + feedback.nA());
		s[2] = (char)('0' + feedback.nB());
	}
	return os << s;
}

/// Inputs the feedback from a stream in the form "1A2B".
/// @ingroup Feedback
inline std::istream& operator >> (std::istream &is, Feedback &feedback)
{
	char s[5];

	// Read the first character, skipping whitespaces.
	if (!(is >> s[0]))
		return is;

	// Read the next three characters without skipping whitespaces.
	for (int i = 1; i <= 3; ++i)
	{
		if ((s[i] = (char)is.get()) == EOF)
			return is;
	}
	s[4] = '\0';

	// Get the rules associated with the stream.
	Rules rules = getrules(is);

	// Parse the input.
	Feedback ret(s);
	if (ret && (!rules || ret.conforming(rules)))
	{
		feedback = ret;
	}
	else
	{
		is.setstate(std::ios_base::failbit);
	}
	return is;
}

} // namespace Mastermind

#endif // MASTERMIND_FEEDBACK_HPP
