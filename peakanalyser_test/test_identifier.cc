#include "pch.h"

#include "../PeakAnalyzer/identifier.h"


class IdentifierTest : public ::testing::Test
{
protected:
	bioscripts::Identifier<bioscripts::Gene> gene_identifier{"AT1G12345"};
	bioscripts::Identifier<bioscripts::Transcript> transcript_identifier{"AT1G12345.1"};
	bioscripts::Identifier<bioscripts::Full> full_identifier{"AT1G12345.1.exon1"};
};

TEST_F(IdentifierTest, to_string_ReturnsCorrectString)
{
	EXPECT_EQ(gene_identifier.to_string(), "AT1G12345");
	EXPECT_EQ(transcript_identifier.to_string(), "AT1G12345.1");
	EXPECT_EQ(full_identifier.to_string(), "AT1G12345.1.exon1");
}

TEST_F(IdentifierTest, gene_ReturnsCorrectString)
{
	EXPECT_EQ(gene_identifier.gene(), "AT1G12345");
	EXPECT_EQ(transcript_identifier.gene(), "AT1G12345");
}

TEST_F(IdentifierTest, version_ReturnsCorrectString)
{
	EXPECT_EQ(transcript_identifier.version(), "1");
}
TEST(CreateTranscriptIdentifier, ThrowExceptionWhenIdentifierIsInvalid) {
	std::string gene_id = "AT1G12345";
	EXPECT_THROW(bioscripts::Identifier<bioscripts::Transcript> gene_identifier{ gene_id }, std::logic_error);
}

TEST(TestIdentifierEquality, SameTranscriptIdentifiers_ComparesEqual)
{
	std::string transcript_id = "AT1G12345.1";
	auto identifier_1 = bioscripts::Identifier<bioscripts::Transcript>{ transcript_id };
	auto identifier_2 = bioscripts::Identifier<bioscripts::Transcript>{ transcript_id };
	EXPECT_EQ(identifier_1, identifier_2);
}

TEST(TestIdentifierEquality, DifferentTranscriptIdentifiers_ComparesNotEqual)
{
	std::string isoform_1 = "AT1G12345.1";
	std::string isoform_2 = "AT1G12345.2";
	auto identifier_1 = bioscripts::Identifier<bioscripts::Transcript>{ isoform_1 };
	auto identifier_2 = bioscripts::Identifier<bioscripts::Transcript>{ isoform_2 };
	EXPECT_NE(identifier_1, identifier_2);
}


