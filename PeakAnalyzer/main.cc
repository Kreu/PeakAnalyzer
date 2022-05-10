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

	std::cout << "Reading in GFF records\n";
	auto gff_records = bioscripts::gff::Records{ gff_file };
	std::cout << "Found " << gff_records.size() << " records\n";
	std::cout << "Reading in peak file\n";
	auto peaks = bioscripts::peak::Peaks{ peaks_file };
	std::cout << "Found " << peaks.size() << " records\n";
	
	std::cout << "Analysing peaks\n";
	for (const auto& peak : peaks) {
		auto midpoint = bioscripts::peak::midpoint(peak);

		//std::cout << "Peak belongs to sequence ID of " << peak.feature.identifier.to_string() << "\n";

		//Find all underlying records on the same chromosome
		auto overlapping_records = gff_records.findUnderlyingRecords(midpoint, peak.sequence_id, bioscripts::gff::Record::Type::CDS);
		std::cout << "Found " << overlapping_records.size() << " records that are at peak midpoint at " << midpoint << "\n";

		if (overlapping_records.empty()) {
			//TODO: Find closest record instead
			continue;
		}

		//Remove those whose identifier is different from the called peak identifier
		auto recordIdentifierDoesNotMatchPeakIdentifier = [&peak](const auto& elem)
		{
			return peak.feature.identifier.to_string() != bioscripts::gff::extractAttribute(elem, "Name");
		};

		std::erase_if(overlapping_records, recordIdentifierDoesNotMatchPeakIdentifier);
		std::cout << "After filtering there are " << overlapping_records.size() << " records left\n";


		//if (overlapping_records.size() != 1) {
		//	throw std::runtime_error("After filtering more than one record remains");
		//}



		//std::cout << "After filtering there are " << overlapping_records.size() << " records left\n";

		//if (overlapping_records.size() != 1) {
		//	throw std::runtime_error("After filtering more than one record remains");
		//}


		//bioscripts::gff::Records::pointer closest_record = nullptr;

		//else if (overlapping_records.empty()) {
		//	closest_record = gff_records.findClosestRecord(midpoint, peak.sequence_id, bioscripts::gff::Record::Type::CDS);
		//}

		////Recreate the full CDS from the peak.
		//auto last_peak_cds = gff_records.findLastRecord(peak.feature.ensembl_id, peak.sequence_id, bioscripts::gff::Record::Type::CDS);




		//Write the output
	}

}