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



class FindRecordTest : public ::testing::Test
{
protected:
	bioscripts::gff::Records records;

	void SetUp()
	{
		auto r1 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::CDS,
		.strand = bioscripts::Strand::Antisense,
		.span = bioscripts::Range{9436697, 9436890},
		.sequence_id = std::string{"4"},
		.attributes = "ID=CDS:AT4G16770.1;Parent=transcript:AT4G16770.1;protein_id=AT4G16770.1"
		};

		auto r1_1 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::CDS,
		.strand = bioscripts::Strand::Antisense,
		.span = bioscripts::Range{9437061, 9437184},
		.sequence_id = std::string{"4"},
		.attributes = "ID=CDS:AT4G16770.1;Parent=transcript:AT4G16770.1;protein_id=AT4G16770.1"
		};

		auto r1_2 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::exon,
		.strand = bioscripts::Strand::Antisense,
		.span = bioscripts::Range{9437061, 9437207},
		.sequence_id = std::string{"4"},
		.attributes = "Parent=transcript:AT4G16770.1;Name=AT4G16770.1.exon1;constitutive=1;ensembl_end_phase=0;ensembl_phase=-1;exon_id=AT4G16770.1.exon1;rank=1"
		};

		auto r1_3 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::five_prime_UTR,
		.strand = bioscripts::Strand::Antisense,
		.span = bioscripts::Range{9437184, 9437207},
		.sequence_id = std::string{"4"},
		.attributes = "Parent=transcript:AT4G16770.1"
		};

		auto r2 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::gene,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9449114, 9450906},
		.sequence_id = std::string{"4"},
		.attributes = "ID=gene:AT4G16780;Name=HAT4;biotype=protein_coding;description=Homeobox-leucine zipper protein HAT4 [Source:UniProtKB/Swiss-Prot%3BAcc:Q05466];gene_id=AT4G16780;logic_name=araport11"
		};

		auto r2_1 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::mRNA,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9449114, 9450906},
		.sequence_id = std::string{"4"},
		.attributes = "ID=transcript:AT4G16780.1;Parent=gene:AT4G16780;biotype=protein_coding;transcript_id=AT4G16780.1"
		};

		auto r2_2 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::five_prime_UTR,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9449114, 9449291},
		.sequence_id = std::string{"4"},
		.attributes = "Parent=transcript:AT4G16780.1"
		};

		auto r2_3 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::exon,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9449114, 9449451},
		.sequence_id = std::string{"4"},
		.attributes = "Parent=transcript:AT4G16780.1;Name=AT4G16780.1.exon1;constitutive=1;ensembl_end_phase=1;ensembl_phase=-1;exon_id=AT4G16780.1.exon1;rank=1"
		};

		auto r2_4 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::CDS,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9449291, 9449451},
		.sequence_id = std::string{"4"},
		.attributes = "ID=CDS:AT4G16780.1;Parent=transcript:AT4G16780.1;protein_id=AT4G16780.1"
		};

		auto r2_5 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::exon,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9449554, 9449853},
		.sequence_id = std::string{"4"},
		.attributes = "Parent=transcript:AT4G16780.1;Name=AT4G16780.1.exon2;constitutive=1;ensembl_end_phase=0;ensembl_phase=1;exon_id=AT4G16780.1.exon2;rank=2"
		};


		auto r2_6 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::CDS,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9449554, 9449853},
		.sequence_id = std::string{"4"},
		.attributes = "ID=CDS:AT4G16780.1;Parent=transcript:AT4G16780.1;protein_id=AT4G16780.1"
		};

		auto r2_7 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::exon,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9450121, 9450201},
		.sequence_id = std::string{"4"},
		.attributes = "Parent=transcript:AT4G16780.1;Name=AT4G16780.1.exon3;constitutive=1;ensembl_end_phase=2;ensembl_phase=0;exon_id=AT4G16780.1.exon3;rank=3"
		};

		auto r2_8 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::CDS,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9450121, 9450201},
		.sequence_id = std::string{"4"},
		.attributes = "ID=CDS:AT4G16780.1;Parent=transcript:AT4G16780.1;protein_id=AT4G16780.1"
		};

		auto r2_9 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::CDS,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9450289, 9450906},
		.sequence_id = std::string{"4"},
		.attributes = "ID=CDS:AT4G16780.1;Parent=transcript:AT4G16780.1;protein_id=AT4G16780.1"
		};

		auto r2_10 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::exon,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9450289, 9450906},
		.sequence_id = std::string{"4"},
		.attributes = "Parent=transcript:AT4G16780.1;Name=AT4G16780.1.exon4;constitutive=1;ensembl_end_phase=-1;ensembl_phase=2;exon_id=AT4G16780.1.exon4;rank=4"
		};

		auto r2_11 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::three_prime_UTR,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9450605, 9450906},
		.sequence_id = std::string{"4"},
		.attributes = "Parent=transcript:AT4G16780.1"
		};



		auto r3 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::gene,
		.strand = bioscripts::Strand::Antisense,
		.span = bioscripts::Range{9451476, 9453301},
		.sequence_id = std::string{"4"},
		.attributes = "ID=gene:AT4G16790;biotype=protein_coding;description=AT4g16790/dl4420c [Source:UniProtKB/TrEMBL%3BAcc:Q9SUK6];gene_id=AT4G16790;logic_name=araport11"
		};

		auto r3_1 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::mRNA,
		.strand = bioscripts::Strand::Antisense,
		.span = bioscripts::Range{9451476, 9453301},
		.sequence_id = std::string{"4"},
		.attributes = "ID=transcript:AT4G16790.1;Parent=gene:AT4G16790;biotype=protein_coding;transcript_id=AT4G16790.1"
		};

		auto r3_2 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::three_prime_UTR,
		.strand = bioscripts::Strand::Antisense,
		.span = bioscripts::Range{9451476, 9451747},
		.sequence_id = std::string{"4"},
		.attributes = "Parent=transcript:AT4G16790.1"
		};

		auto r3_3 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::exon,
		.strand = bioscripts::Strand::Antisense,
		.span = bioscripts::Range{9451476, 9453301},
		.sequence_id = std::string{"4"},
		.attributes = "Parent=transcript:AT4G16790.1;Name=AT4G16790.1.exon1;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=AT4G16790.1.exon1;rank=1"
		};

		auto r3_4 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::CDS,
		.strand = bioscripts::Strand::Antisense,
		.span = bioscripts::Range{9451747, 9453169},
		.sequence_id = std::string{"4"},
		.attributes = "ID=CDS:AT4G16790.1;Parent=transcript:AT4G16790.1;protein_id=AT4G16790.1"
		};

		auto r3_5 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::five_prime_UTR,
		.strand = bioscripts::Strand::Antisense,
		.span = bioscripts::Range{9453169, 9453301},
		.sequence_id = std::string{"4"},
		.attributes = "Parent=transcript:AT4G16790.1"
		};

		auto r4 = bioscripts::gff::Record{
		.type = bioscripts::gff::Record::Type::CDS,
		.strand = bioscripts::Strand::Sense,
		.span = bioscripts::Range{9467563, 9469117},
		.sequence_id = std::string{"4"},
		.attributes = "ID=CDS:AT4G16820.1;Parent=transcript:AT4G16820.1;protein_id=AT4G16820.1"
		};


		records.add(r1);
		records.add(r1_1);
		records.add(r1_2);
		records.add(r1_3);
		records.add(r2);
		records.add(r2_1);
		records.add(r2_2);
		records.add(r2_3);
		records.add(r2_4);
		records.add(r2_5);
		records.add(r2_6);
		records.add(r2_7);
		records.add(r2_8);
		records.add(r2_9);
		records.add(r2_10);
		records.add(r2_11);
		records.add(r3);
		records.add(r3_1);
		records.add(r3_2);
		records.add(r3_3);
		records.add(r3_4);
		records.add(r3_5);
		records.add(r4);
	}
};

TEST_F(FindRecordTest, FindClosestRecordNearCurrentRecord_ReturnClosestRecord)
{
	bioscripts::Position pos_after_last_CDS = 9450920;
	bioscripts::Position pos_between_two_records = 9450200;
	bioscripts::Position pos_closest_to_last_CDS_record = 9450260;

	auto sequence_id = bioscripts::Identifier<bioscripts::Full>{"4"};
	auto identifier_to_search = bioscripts::Identifier<bioscripts::Gene>{ "AT4G16780.1"};

	auto sequence_records = records.data(sequence_id);

	auto found_record = records.findClosestRecord(pos_after_last_CDS, sequence_id, identifier_to_search, bioscripts::gff::Record::Type::CDS);
	auto found_record2 = records.findClosestRecord(pos_between_two_records, sequence_id, identifier_to_search, bioscripts::gff::Record::Type::CDS);
	auto found_record3 = records.findClosestRecord(pos_closest_to_last_CDS_record, sequence_id, identifier_to_search, bioscripts::gff::Record::Type::CDS);

	auto expected_record = sequence_records[13];
	auto expected_record2 = sequence_records[12];

	EXPECT_EQ(*found_record, expected_record);
	EXPECT_EQ(*found_record2, expected_record2);
	EXPECT_EQ(*found_record3, expected_record);
};

TEST_F(FindRecordTest, CollectCodingSequences_FindAllCodingSequencesOnSenseStrand_ReturnCorrectNumberOfRecords)
{
	auto sequence_id = bioscripts::Identifier<bioscripts::Full>{ "4" };
	auto sequence_records = records.data(sequence_id);
	auto starting_record = sequence_records[12];
	auto all_coding_sequence_records = bioscripts::gff::collectCodingSequenceRecords(starting_record, records);
	EXPECT_EQ(all_coding_sequence_records.size(), 2);

	starting_record = sequence_records[8];
	all_coding_sequence_records = bioscripts::gff::collectCodingSequenceRecords(starting_record, records);
	EXPECT_EQ(all_coding_sequence_records.size(), 4);
}

TEST_F(FindRecordTest, CollectCodingSequences_FindAllCodingSequencesOnSenseStrand_SequencesAreInCorrectOrder)
{
	auto sequence_id = bioscripts::Identifier<bioscripts::Full>{ "4" };
	auto sequence_records = records.data(sequence_id);
	auto starting_record = sequence_records[12];
	auto all_coding_sequence_records = bioscripts::gff::collectCodingSequenceRecords(starting_record, records);
	EXPECT_EQ(all_coding_sequence_records[0], starting_record);
	EXPECT_EQ(all_coding_sequence_records[1], sequence_records[13]);

	starting_record = sequence_records[8];
	all_coding_sequence_records = bioscripts::gff::collectCodingSequenceRecords(starting_record, records);
	EXPECT_EQ(all_coding_sequence_records[0], starting_record);
	EXPECT_EQ(all_coding_sequence_records[1], sequence_records[10]);
	EXPECT_EQ(all_coding_sequence_records[2], sequence_records[12]);
	EXPECT_EQ(all_coding_sequence_records[3], sequence_records[13]);
}

TEST_F(FindRecordTest, CollectCodingSequences_FindAllCodingSequencesOnAntiSenseStrand_ReturnCorrectNumberOfRecords)
{
	auto sequence_id = bioscripts::Identifier<bioscripts::Full>{ "4" };
	auto sequence_records = records.data(sequence_id);
	auto starting_record = sequence_records[1];
	auto all_coding_sequence_records = bioscripts::gff::collectCodingSequenceRecords(starting_record, records);
	EXPECT_EQ(all_coding_sequence_records.size(), 2);
}

TEST_F(FindRecordTest, CollectCodingSequences_FindAllCodingSequencesOnAntiSenseStrand_SequencesAreInCorrectOrder)
{
	auto sequence_id = bioscripts::Identifier<bioscripts::Full>{ "4" };
	auto sequence_records = records.data(sequence_id);
	auto starting_record = sequence_records[1];
	auto all_coding_sequence_records = bioscripts::gff::collectCodingSequenceRecords(starting_record, records);
	//EXPECT_EQ(all_coding_sequence_records[0], sequence_records[1]);
	//EXPECT_EQ(all_coding_sequence_records[1], sequence_records[0]);
}


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
