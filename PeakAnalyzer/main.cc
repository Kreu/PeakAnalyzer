#include "gff.h"
#include "peak.h"

int main()
{
	auto peaks_file = "peaks.txt";
	auto gff_file = "gff_records.txt";

	auto gff_records = bioscripts::gff::Records{ gff_file };
	auto peaks = bioscripts::peak::Peaks{ peaks_file };
	
	for (const auto& peak : peaks) {

	}

}