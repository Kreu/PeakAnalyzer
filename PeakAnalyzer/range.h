#ifndef BIOSCRIPTS_RANGE_H
#define BIOSCRIPTS_RANGE_H

#include <cmath>
#include <cstddef>

namespace bioscripts
{
	using Position = std::size_t;

	/**
	 * @brief  Represents a 0-based half-open range [start, end).
	 * @pre  @a start >= 0 && @a end >= 0
	 * @pre  @a start < @a end.
	 */
	struct Range
	{
		Range(std::size_t start, std::size_t end);
		Position start;
		Position end;
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
	bool overlap(const Position& p, const Range& r);

	/**
	 * @brief  Returns the centre position of the range
	 * 
	 * Return value will be rounded to the nearest integer.
	 */
	Position centre(const Range& r);
}
#endif // !BIOSCRIPTS_RANGE_H
