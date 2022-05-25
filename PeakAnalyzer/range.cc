#include <algorithm>

#include "range.h"

namespace bioscripts
{
	Range::Range(std::size_t start, std::size_t end) : start(start), end(end) {}

	bool overlap(const Range& r1, const Range& r2)
	{
		return std::max(r1.start, r2.start) < std::min(r1.end, r2.end);
	}

	bool overlap(const Point& p, const Range& r)
	{
		return r.start >= p && r.end < p;
	}

}