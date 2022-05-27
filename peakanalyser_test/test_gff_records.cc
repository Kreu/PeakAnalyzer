#include "pch.h"

#include "../PeakAnalyzer/gff.h"


class RecordTest : public ::testing::Test
{
protected:
	void SetUp()
	{
		auto record_1 = bioscripts::gff::Record{
			.type = bioscripts::gff::Record::Type::CDS,
			.strand = bioscripts::Strand::Sense,
			.span = bioscripts::Range{0, 100},
			.sequence_id = std::string{"Chromosome_1"},
			.attributes = "ID=CDS:ATMG00180.1;Parent=transcript:ATMG00180.1;protein_id=ATMG00180.1"
		};

		auto record_2 = bioscripts::gff::Record{
			.type = bioscripts::gff::Record::Type::mRNA,
			.strand = bioscripts::Strand::Sense,
			.span = bioscripts::Range{0, 300},
			.sequence_id = std::string{"Chromosome_1"},
			.attributes = "ID=CDS:ATMG00180.1;Parent=transcript:ATMG00180.1;protein_id=ATMG00180.1"
		};

		records.add(record_1);
		records.add(record_2);
	}

	//bioscripts::Identifier<bioscripts::Gene> gene_identifier{ "AT1G12345" };
	//bioscripts::Identifier<bioscripts::Transcript> transcript_identifier{ "AT1G12345.1" };
	//bioscripts::Identifier<bioscripts::Full> full_identifier{ "AT1G12345.1.exon1" };

	bioscripts::gff::Records records;
};

TEST_F(RecordTest, getRecordsAt_RecordExistsAtGivenGenomicPosition_ReturnAllUnderlyingRecords)
{
	auto found_records = records.getRecordsAt(50, bioscripts::Identifier<bioscripts::Full>{ "Chromosome_1" });
	EXPECT_EQ(found_records.size(), 2);
};

TEST_F(RecordTest, getRecordsAt_NoRecordAtGivenGenomicPosition_ReturnEmptyCollection)
{
	auto found_records = records.getRecordsAt(500, bioscripts::Identifier<bioscripts::Full>{ "Chromosome_1" });
	EXPECT_EQ(found_records.size(), 0);
};

TEST_F(RecordTest, getRecordsAt_RecordAtPositionIsWrongType_ReturnEmptyCollection)
{
	auto found_records = records.getRecordsAt(50, bioscripts::Identifier<bioscripts::Full>{ "Chromosome_1" }, bioscripts::gff::Record::Type::gene);
	EXPECT_EQ(found_records.size(), 0);
};

TEST_F(RecordTest, getRecordsAt_RecordAtPositionButNotOnSameSequence_ReturnEmptyCollection)
{
	auto found_records = records.getRecordsAt(50, bioscripts::Identifier<bioscripts::Full>{ "Chromosome_Unknown" });
	EXPECT_EQ(found_records.size(), 0);
};


//TEST(findClosestRecordTest, FindClosestRecordBeforeCurrentRecord_ReturnRecordWithClosestEndPoint)
//{
//
//}
//
//
//TEST(findClosestRecordTest, FindClosestRecordBeforeCurrentRecord_ReturnRecordWithClosestEndPoint)
//{
//
//}
//
//

TEST(RecordsTest, extractAttribute_ExtractIdentifiersFromAttributes_ReturnCorrectIdentifiers)
{
	auto record = bioscripts::gff::Record{
	.type = bioscripts::gff::Record::Type::CDS,
	.strand = bioscripts::Strand::Sense,
	.span = bioscripts::Range{0, 100},
	.sequence_id = std::string{"Chromosome_1"},
	.attributes = "ID=CDS:ATMG00180.1;Parent=transcript:ATMG00180.1;protein_id=ATMG00180.1"
	};

	auto transcript_id = bioscripts::gff::extractAttribute(record, "ID=CDS");
	EXPECT_EQ(transcript_id, "ATMG00180.1");

	auto parent_id = bioscripts::gff::extractAttribute(record, "ID=CDS");
	EXPECT_EQ(parent_id, "ATMG00180.1");
	
	auto protein_id = bioscripts::gff::extractAttribute(record, "protein_id");
	EXPECT_EQ(protein_id, "ATMG00180.1");
}

TEST(RecordsTest, add_AddRecord_SizeIncreaseByOne)
{
	bioscripts::gff::Records records;

	auto record = bioscripts::gff::Record{
	.type = bioscripts::gff::Record::Type::CDS,
	.strand = bioscripts::Strand::Sense,
	.span = bioscripts::Range{0, 100},
	.sequence_id = std::string{""},
	.attributes = ""
	};

	records.add(record);
	EXPECT_EQ(records.size(), 1);
}

TEST(RecordsTest, size_EmptyRecordsObject_SizeIsZero)
{
	bioscripts::gff::Records records;
	EXPECT_EQ(records.size(), 0);
}
