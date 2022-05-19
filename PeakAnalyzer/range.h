#ifndef BIOSCRIPTS_RANGE_H
#define BIOSCRIPTS_RANGE_H

#include <cstddef>

namespace bioscripts
{
	/**
	 * @brief  Represents an a half-open range [start, end).
	 * @pre  @a start < @a end.
	 */
	struct Range
	{
		Range(std::size_t start, std::size_t end);
		std::size_t start;
		std::size_t end;
	};

	/**
	 * @brief  Determines whether two ranges overlap.
	 * @return  True if the ranges overlap, false otherwise.
	 */
	bool doRangesOverlap(const Range& r1, const Range& r2);
}
#endif // !BIOSCRIPTS_RANGE_H
