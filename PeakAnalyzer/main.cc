#include "gff.h"
#include "peak.h"

int main()
{
	auto peaks_file = "peaks.txt";
	auto gff_file = "gff_records.txt";

	auto gff_records = bioscripts::gff::Records{ gff_file };
	auto peaks = bioscripts::peak::Peaks{ peaks_file };
	
	for (const auto& peak : peaks) {
		auto midpoint = bioscripts::peak::midpoint(peak);

		//Find all underlying records on the same chromosome
		auto overlapping_records = gff_records.findUnderlyingRecords(midpoint, peak.sequence_id, bioscripts::gff::Record::Type::CDS);
		if (overlapping_records.empty()) {
			overlapping_records = gff_records.findClosestRecord(midpoint, peak.sequence_id, bioscripts::gff::Record::Type::CDS);
		}

		//Recreate the full CDS from the peak.
		auto last_peak_cds = gff_records.findLastRecord(peak.feature, peak.sequence_id, bioscripts::gff::Record::Type::CDS);





		//Write the output
	}

}