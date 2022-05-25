#ifndef BIOSCRIPTS_RANGE_H
#define BIOSCRIPTS_RANGE_H

#include <cstddef>

namespace bioscripts
{
	using Point = std::size_t;

	/**
	 * @brief  Represents an a half-open range [start, end).
	 * @pre  @a start >= 0 && @a end >= 0
	 * @pre  @a start < @a end.
	 */
	struct Range
	{
		Range(std::size_t start, std::size_t end);
		Point start;
		Point end;
	};

	/**
	 * @brief  Determines whether two ranges overlap.
	 * @return  True if the ranges overlap, false otherwise.
	 */
	bool overlap(const Range& r1, const Range& r2);

	/**
	 * @brief  Determines whether a point is contained within a range.
	 * @return  True if @a p is contained within @a r, else return false.
	 */
	bool overlap(const Point& p, const Range& r);
}
#endif // !BIOSCRIPTS_RANGE_H
