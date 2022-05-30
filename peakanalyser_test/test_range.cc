#include "pch.h"

#include "../PeakAnalyzer/range.h"


TEST(TestRange, overlap_TwoRangesOverlap_ReturnsTrue)
{
	auto r1 = bioscripts::Range{ 0, 100 };
	auto r2 = bioscripts::Range{ 50, 100 };
	auto r3 = bioscripts::Range{ 99, 500 };
	EXPECT_EQ(bioscripts::overlap(r1, r2), true);
	EXPECT_EQ(bioscripts::overlap(r1, r3), true);
}

TEST(TestRange, overlap_TwoRangesDoNotOverlap_ReturnsFalse)
{
	auto r1 = bioscripts::Range{ 0, 100 };
	auto r2 = bioscripts::Range{ 100, 150 };
	EXPECT_EQ(bioscripts::overlap(r1, r2), false);
}

TEST(TestRange, overlap_PositionIsWithinRange_ReturnsTrue)
{
	auto r1 = bioscripts::Range{ 0, 100 };
	bioscripts::Position position = 50;
	bioscripts::Position position2 = 99;
	EXPECT_EQ(bioscripts::overlap(position, r1), true);
	EXPECT_EQ(bioscripts::overlap(position2, r1), true);
}

TEST(TestRange, overlap_PositionIsNotWithinRange_ReturnsFalse)
{
	auto r1 = bioscripts::Range{ 0, 100 };
	bioscripts::Position position = 100;
	bioscripts::Position position2 = 300;
	EXPECT_EQ(bioscripts::overlap(position, r1), false);
	EXPECT_EQ(bioscripts::overlap(position2, r1), false);
}

TEST(TestRange, centre_CentreIsWholeInteger_ReturnsExpectedCentre)
{
	auto r1 = bioscripts::Range{ 0, 5 };
	EXPECT_EQ(bioscripts::centre(r1), 2);
}

TEST(TestRange, centre_CentreIsBetweenTwoNumbers, ReturnsCentreRoundedUp)
{
	auto r1 = bioscripts::Range{ 0, 6 };
	auto r2 = bioscripts::Range{ 2, 7 };
	EXPECT_EQ(bioscripts::centre(r1), 3);
	EXPECT_EQ(bioscripts::centre(r2), 4);
}

TEST(TestRange, length_NonZeroLengthRange_ReturnsCorrectLength)
{
	auto r1 = bioscripts::Range{ 0, 3 };
	auto r2 = bioscripts::Range{ 2, 8 };
	EXPECT_EQ(bioscripts::length(r1), 3);
	EXPECT_EQ(bioscripts::length(r2), 6);
}

TEST(TestRange, distance_TwoRangesDoNotOverlap_ReturnsCorrectDistance)
{
	auto r1 = bioscripts::Range{ 0, 4 };
	auto r2 = bioscripts::Range{ 7, 12 };
	EXPECT_EQ(bioscripts::distance(r1, r2), 3);
}

TEST(TestRange, distance_ParameterOrderDoesNotInfluenceResult_ReturnsSameDistance)
{
	auto r1 = bioscripts::Range{ 0, 4 };
	auto r2 = bioscripts::Range{ 7, 12 };
	EXPECT_EQ(bioscripts::distance(r1, r2), 3);
	EXPECT_EQ(bioscripts::distance(r2, r1), 3);
}

TEST(TestRange, distance_TwoRangesOverlap_ReturnsZero)
{
	auto r1 = bioscripts::Range{ 0, 4 };
	auto r2 = bioscripts::Range{ 3, 10 };
	EXPECT_EQ(bioscripts::distance(r1, r2), 0);
}

TEST(TestRange, distance_TwoRangesAreAdjacent_ReturnsZero)
{
	auto r1 = bioscripts::Range{ 0, 4 };
	auto r2 = bioscripts::Range{ 4, 10 };
	EXPECT_EQ(bioscripts::distance(r1, r2), 0);
}

TEST(TestRange, distance_CalculateDistanceToPrecedingRecord_ReturnsDistanceToEndOfRecord)
{
	auto r2 = bioscripts::Range{ 0, 5 };
	auto p = bioscripts::Position{ 10 };

	EXPECT_EQ(distance(p, r2), 6);
}

TEST(TestRange, distance_CalculateDistanceToFollowingRecord_ReturnsDistanceToStartOfRecord)
{
	auto r2 = bioscripts::Range{ 20, 30 };
	auto p = bioscripts::Position{ 10 };

	EXPECT_EQ(distance(p, r2), 10);
}


//TEST(TestRange, length_ZeroLengthRange, ReturnsZero)
//{
//
//}