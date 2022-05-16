#include "gff.h"
#include "peak.h"

#include <iostream>t


int main(int argc, char* argv[])
{
	if (argc != 3) {
		std::cerr << "Unknown arguments deteced.\n";
		std::cerr << "Usage: " << argv[0] << " [peaks_file] [gff_file]\n";
		return 1;
	}

	auto peaks_file = argv[1];
	auto gff_file = argv[2];

	//For the analysis we only need CDS records.
	auto isNotCdsRecord = [](const auto& elem)
	{
		return elem.type != bioscripts::gff::Record::Type::CDS;
	};
	auto gff_records = bioscripts::gff::Records{ gff_file };

	for (auto& [id, records] : gff_records) {
		std::cout << "Size before deleting: " << records.size() << "\n";
		std::erase_if(records, isNotCdsRecord);
		std::cout << "Size after deleting: " << records.size() << "\n";
	}
	std::cout << "Finished deleting all records that is not CDS\n";

	auto peaks = bioscripts::peak::Peaks{ peaks_file };
	
	std::cout << "Analysing peaks\n";
	for (const auto& peak : peaks)
	{
		auto midpoint = bioscripts::peak::midpoint(peak);

		auto overlapping_records = gff_records.findUnderlyingRecords(midpoint, peak.sequence_id, bioscripts::gff::Record::Type::CDS);

		//Remove those whose identifier is different from the called peak identifier
		auto recordIdentifierDoesNotMatchPeakIdentifier = [&peak](const auto& elem)
		{
			return peak.feature.identifier != bioscripts::Identifier<bioscripts::Transcript>{ bioscripts::gff::extractAttribute(elem, "ID=CDS") };
		};
		
		std::erase_if(overlapping_records, recordIdentifierDoesNotMatchPeakIdentifier);
		//std::cout << "After filtering there are " << overlapping_records.size() << " records left\n";

		//std::vector<bioscripts::gff::Record> subsequent_gff_records;
		if (overlapping_records.empty()) {
			auto closest_record = gff_records.findClosestRecord(midpoint, peak.sequence_id, peak.feature.identifier, bioscripts::gff::Record::Type::CDS);
			if (closest_record == nullptr) {
				continue;
			}
			auto all_cds_records = bioscripts::gff::collectCodingSequenceRecords(*closest_record, gff_records);
			//std::cout << "Found another " << all_cds_records.size() << " exons corresponding to the peak CDS\n";
			//if (all_cds_records.size() >= 10) {
			//	std::cout << closest_record->attributes << "\n";
			//}

		}
		else {
			for (auto& record : overlapping_records) {
				auto all_cds_records = bioscripts::gff::collectCodingSequenceRecords(record, gff_records);
				//std::cout << "Found another " << all_cds_records.size() << " exons corresponding to the peak CDS\n";
				//if (all_cds_records.size() >= 10) {
				//	std::cout << record.attributes << "\n";
				//}
			}
		}



		//Recreate the full CDS from the peak.
		//auto last_peak_cds = gff_records.findLastRecord(peak.feature.ensembl_id, peak.sequence_id, bioscripts::gff::Record::Type::CDS);




		//Write the output
	}

}