#ifndef BIOSCRIPTS_RANGE_H
#define BIOSCRIPTS_RANGE_H

#include <cmath>
#include <cstddef>

namespace bioscripts
{
	using Position = std::size_t;
	using Length = std::size_t;
	using Distance = std::size_t;

	/**
	 * @brief  Represents a 0-based half-open range [start, end).
	 * @pre  @a start >= 0 && @a end >= 0
	 * @pre  @a start < @a end.
	 */
	struct Range
	{
		////Creates an empty range.
		//Range();

		Range(std::size_t start, std::size_t end);
		Position start;
		Position end;
		auto operator<=>(const Range& rhs) const = default;
	};

	/**
	 * @return  The length of the range.
	 */
	Length length(const Range& r);

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
	 * @brief  Calculates the centre position of the range
	 * @return  Position of the central point. Position will be rounded to the nearest integer.
	 */
	Position centre(const Range& r);

	/**
	 * @brief  Calculates the distance between two ranges.
	 * @return  Distance between the two ranges. If the ranges overlap, returns 0.
	 */
	Distance distance(const Range& r1, const Range& r2);

	Distance distance(const Position& p, const Range& r);


}
#endif // !BIOSCRIPTS_RANGE_H
