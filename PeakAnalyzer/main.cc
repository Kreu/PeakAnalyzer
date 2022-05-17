#include "gff.h"
#include "peak.h"

#include <filesystem>
#include <fstream>
#include <iostream>

struct TranscriptData
{
	std::size_t id;
	std::size_t start_pos;
	std::size_t end_pos;
	std::string transcript_id;
};

namespace
{
	void writeOutputFile(const std::filesystem::path& output_file, const std::vector<TranscriptData>& data) {
		std::ofstream of;
		of.open(output_file);
		for (const auto& elem : data) {
			of << elem.id << "\t" << elem.transcript_id << "\t" << elem.start_pos << "\t" << elem.end_pos << "\n";
		}
	}
}


int main(int argc, char* argv[])
{
	if (argc != 3) {
		std::cerr << "Unknown arguments deteced.\n";
		std::cerr << "Usage: " << argv[0] << " [peaks_file] [gff_file]\n";
		return 1;
	}

	auto gff_file = argv[2];
	auto gff_records = bioscripts::gff::Records{ gff_file };

	//We are only interested in the CDS records.
	auto isNotCdsRecord = [](const auto& elem)
	{
		return elem.type != bioscripts::gff::Record::Type::CDS;
	};

	for (auto& [id, records] : gff_records) {
		std::erase_if(records, isNotCdsRecord);
	}
	std::cout << "Finished deleting all records that is not CDS\n";

	auto peaks_file = argv[1];
	auto peaks = bioscripts::peak::Peaks{ peaks_file };

	std::vector<TranscriptData> data_to_write;
	std::cout << "Analysing peaks\n";
	std::size_t peak_id = 0;
	for (const auto& peak : peaks)
	{
		auto midpoint = bioscripts::peak::midpoint(peak);

		auto records_under_the_peak = gff_records.findUnderlyingRecords(midpoint, peak.sequence_id, bioscripts::gff::Record::Type::CDS);

		//Remove those whose identifier is different from the called peak identifier
		auto recordIdentifierDoesNotMatchPeakIdentifier = [&peak](const auto& elem)
		{
			return peak.feature.identifier != bioscripts::Identifier<bioscripts::Transcript>{ bioscripts::gff::extractAttribute(elem, "ID=CDS") };
		};
		
		std::erase_if(records_under_the_peak, recordIdentifierDoesNotMatchPeakIdentifier);

		//If there are no records underneath the peak midpoint, find the closest record instead

		if (records_under_the_peak.empty()) {
			auto closest_record = gff_records.findClosestRecord(midpoint, peak.sequence_id, peak.feature.identifier, bioscripts::gff::Record::Type::CDS);
			if (closest_record == nullptr) {
				continue;
			}

			auto all_cds_records = bioscripts::gff::collectCodingSequenceRecords(*closest_record, gff_records);
			for (const auto& rec : all_cds_records) {
				data_to_write.push_back(TranscriptData{
					.id = peak_id,
					.start_pos = rec.start_pos,
					.end_pos = rec.end_pos,
					.transcript_id = bioscripts::gff::extractAttribute(rec, "ID=CDS")
					});
			}
		}
		else {
			for (auto& record : records_under_the_peak) {
				auto all_cds_records = bioscripts::gff::collectCodingSequenceRecords(record, gff_records);

				for (const auto& rec : all_cds_records) {
					data_to_write.push_back(TranscriptData{
						.id = peak_id,
						.start_pos = rec.start_pos,
						.end_pos = rec.end_pos,
						.transcript_id = bioscripts::gff::extractAttribute(rec, "ID=CDS")
						});
				}
			}
		}
		++peak_id;
	}
	std::cout << "Data to write: " << data_to_write.size() << "\n";
	writeOutputFile("transcript_data.txt", data_to_write);
}