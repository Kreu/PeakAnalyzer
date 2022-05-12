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
	for (const auto& peak : peaks)
	{
		auto midpoint = bioscripts::peak::midpoint(peak);

		auto overlapping_records = gff_records.findUnderlyingRecords(midpoint, peak.sequence_id, bioscripts::gff::Record::Type::CDS);

		//Remove those whose identifier is different from the called peak identifier
		auto recordIdentifierDoesNotMatchPeakIdentifier = [&peak](const auto& elem)
		{
			return peak.feature.identifier.to_string() != bioscripts::Identifier<bioscripts::Transcript>{ bioscripts::gff::extractAttribute(elem, "ID=CDS") };
		};
		
		std::erase_if(overlapping_records, recordIdentifierDoesNotMatchPeakIdentifier);
		//std::cout << "After filtering there are " << overlapping_records.size() << " records left\n";

		if (overlapping_records.empty()) {
			auto closest_record = gff_records.findClosestRecord(midpoint, peak.sequence_id, peak.feature.identifier, bioscripts::gff::Record::Type::CDS);
			if (closest_record == nullptr) {
				continue;
			}
		}

		//Recreate the full CDS from the peak.
		//auto last_peak_cds = gff_records.findLastRecord(peak.feature.ensembl_id, peak.sequence_id, bioscripts::gff::Record::Type::CDS);




		//Write the output
	}

}