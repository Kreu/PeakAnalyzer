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

	auto gff_records = bioscripts::gff::Records{ gff_file };
	auto peaks = bioscripts::peak::Peaks{ peaks_file };
	
	for (const auto& peak : peaks) {
		auto midpoint = bioscripts::peak::midpoint(peak);

		//Find all underlying records on the same chromosome
		auto overlapping_records = gff_records.findUnderlyingRecords(midpoint, peak.sequence_id, bioscripts::gff::Record::Type::CDS);
		bioscripts::gff::Records::pointer closest_record = nullptr;

		auto recordIdentifierDoesNotMatchPeakIdentifier = [&peak](const auto& elem)
		{
			return peak.feature.ensembl_id.to_string() != bioscripts::gff::extractAttribute(elem, "Name");
		};

		if (!overlapping_records.empty()) {
			std::erase_if(overlapping_records, recordIdentifierDoesNotMatchPeakIdentifier);
		}



		else if (overlapping_records.empty()) {
			closest_record = gff_records.findClosestRecord(midpoint, peak.sequence_id, bioscripts::gff::Record::Type::CDS);
		}

		//Recreate the full CDS from the peak.
		auto last_peak_cds = gff_records.findLastRecord(peak.feature.ensembl_id, peak.sequence_id, bioscripts::gff::Record::Type::CDS);




		//Write the output
	}

}