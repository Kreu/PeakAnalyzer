#include <algorithm>

#include "range.h"

namespace bioscripts
{
	//Range::Range() : start(0), end(1) {}

	Range::Range(std::size_t start, std::size_t end) : start(start), end(end) {}

	Length length(const Range& r)
	{
		return (r.end - r.start);
	}

	bool overlap(const Range& r1, const Range& r2)
	{
		return std::max(r1.start, r2.start) < std::min(r1.end, r2.end);
	}

	bool overlap(const Position& p, const Range& r)
	{
		return r.start <= p && r.end > p;
	}

	Position centre(const Range& r)
	{
		//1 is removed from the end because the range is [start, end).
		//E.g. a range covering [0,5) has a centre at position 2, not 2.5
		auto median = (r.start + r.end - 1) / 2.0;
		return std::llround(median);
	}

	Distance distance(const Range& r1, const Range& r2)
	{
		if (overlap(r1, r2)) {
			return 0;
		}
		return std::max(r1.start, r2.start) - std::min(r1.end, r2.end);
	}

	Distance distance(const Position& p, const Range& r)
	{
		if (overlap(p, r)) {
			return 0;
		}

		std::size_t distance_to_record = 0;
		if (r.end < p) {
			//We add one because range is half-closed [start, end)
			distance_to_record = p - r.end + 1;
		}
		else if (r.start > p) {
			distance_to_record = r.start - p;
		}
		return distance_to_record;
	}


}