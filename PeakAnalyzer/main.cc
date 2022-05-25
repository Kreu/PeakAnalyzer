#include "gff.h"
#include "peak.h"

#include <filesystem>
#include <fstream>
#include <iostream>

#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

struct TranscriptData
{
	std::size_t id;
	std::size_t start_pos;
	std::size_t end_pos;
	std::string transcript_id;
};

namespace
{
	void configureLogger(bool debug_enabled)
	{
		el::Configurations default_logging;
		default_logging.setToDefault();
		if (debug_enabled) {
			default_logging.set(el::Level::Debug, el::ConfigurationType::Enabled, "true");
			default_logging.set(el::Level::Debug, el::ConfigurationType::Filename, "debug.log");
			default_logging.set(el::Level::Debug, el::ConfigurationType::ToStandardOutput, "false");
		}
		default_logging.set(el::Level::Info, el::ConfigurationType::ToStandardOutput, "false");
		default_logging.set(el::Level::Info, el::ConfigurationType::Filename, "peaks.log");
		default_logging.set(el::Level::Info, el::ConfigurationType::Format, "%datetime %level %msg");
		el::Loggers::reconfigureLogger("default", default_logging);
	}

	void writeOutputFile(const std::filesystem::path& output_file, const std::vector<TranscriptData>& data) {
		std::ofstream of;
		of.open(output_file);
		for (const auto& elem : data) {
			of << elem.id << "\t" << elem.transcript_id << "\t" << elem.start_pos << "\t" << elem.end_pos << "\n";
		}
	}


	void analysePeaks()
	{

	};


}


int main(int argc, char* argv[])
{
	if (argc != 3) {
		std::cerr << "Unknown arguments deteced.\n";
		std::cerr << "Usage: " << argv[0] << " [peaks_file] [gff_file]\n";
		return 1;
	}

	configureLogger(true);

	auto gff_file = argv[2];
	auto gff_records = bioscripts::gff::Records{ gff_file };
	//We are only interested in CDS records because we want to reconsitute the protein-coding parts and nothing else
	auto cds_gff_records = bioscripts::gff::fetchRecords(gff_records, bioscripts::gff::Record::Type::CDS);

	LOG(INFO) << "Parsing peak file";
	auto peaks_file = argv[1];
	auto peaks = bioscripts::peak::Peaks{ peaks_file };


	std::vector<TranscriptData> data_to_write;
	LOG(INFO) << "Analysing peaks";
	std::size_t peak_id = 0;
	const auto total_nr_of_peaks = peaks.size();
	for (const auto& peak : peaks)
	{
		LOG(DEBUG) << "Analysing peak " << peak_id << "\\" << total_nr_of_peaks;

		auto midpoint = bioscripts::peak::midpoint(peak);
		//LOG(DEBUG) << "Analysing peak with gene ID: " << peak.feature.identifier.to_string() << ", midpoint at " << midpoint << "\n";

		auto records_under_the_peak = gff_records.getRecordsAt(midpoint, peak.sequence_id, bioscripts::gff::Record::Type::CDS);

		//Remove those whose identifier is different from the called peak identifier
		auto recordIdentifierDoesNotMatchPeakIdentifier = [&peak](const auto& elem)
		{
			return peak.associated_identifier != bioscripts::Identifier<bioscripts::Transcript>{ bioscripts::gff::extractAttribute(elem, "ID=CDS") };
		};
		std::erase_if(records_under_the_peak, recordIdentifierDoesNotMatchPeakIdentifier);
		//LOG(DEBUG) << records_under_the_peak.size() << " GFF records found under the peak";

		//If there are no records underneath the peak midpoint, find the closest record instead

		if (records_under_the_peak.empty()) {
			auto closest_record = gff_records.findClosestRecord(midpoint, peak.sequence_id, peak.associated_identifier, bioscripts::gff::Record::Type::CDS);
			if (closest_record == nullptr) {
				continue;
			}
			//LOG(DEBUG) << "Closest record to peak has attributes " << closest_record->attributes;
			auto all_cds_records = bioscripts::gff::collectCodingSequenceRecords(*closest_record, gff_records);
			//LOG(DEBUG) << "Writing all CDS records belonging to the same gene as peak";
			for (const auto& rec : all_cds_records) {
				//LOG(DEBUG) << "Record attribute: " << rec.attributes;
				data_to_write.push_back(TranscriptData{
					.id = peak_id,
					.start_pos = rec.start(),
					.end_pos = rec.end(),
					.transcript_id = bioscripts::gff::extractAttribute(rec, "ID=CDS")
					});
			}
		}
		else {
			for (auto& record : records_under_the_peak) {
				auto all_cds_records = bioscripts::gff::collectCodingSequenceRecords(record, gff_records);
				//LOG(DEBUG) << "Writing all CDS records belonging to the same gene as peak";
				for (const auto& rec : all_cds_records) {
					//LOG(DEBUG) << "Record attribute: " << rec.attributes;

					data_to_write.push_back(TranscriptData{
						.id = peak_id,
						.start_pos = rec.start(),
						.end_pos = rec.end(),
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