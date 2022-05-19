#ifndef BIOSCRIPTS_PEAK_H
#define BIOSCRIPTS_PEAK_H

#include <cstdint>
#include <filesystem>
#include <string>

#include "identifier.h"
#include "strand.h"

namespace bioscripts
{
	namespace peak
	{
		struct Feature
		{
			std::size_t start_pos;
			std::size_t end_pos;
			Strand strand;
			Identifier<Gene> identifier; //Denotes the identifier under the peak. Usually it is the gene name, e.g. AT1G12345
			Identifier<Transcript> ensembl_id; //This is the identifier associate with the peak after R analysis. More detailed "identifier", e.g. AT1G12345.2.exon4
		};

		struct Peak
		{
			std::size_t start_pos;
			std::size_t end_pos;
			uint16_t width;
			uint16_t score;
			double pvalue;
			double qvalue;
			std::string sequence_id; //Denotes the ID of the sequence the peak is found on. Usually it is the chromosome identifier.
			Feature feature;
		};

		class Peaks
		{
		public:
			using iterator = std::vector<Peak>::iterator;
			using const_iterator = std::vector<Peak>::const_iterator;

			Peaks(const std::filesystem::path& peaks_file);

			iterator begin();
			const_iterator begin() const;
			const_iterator cbegin() const;

			iterator end();
			const_iterator end() const;
			const_iterator cend() const;

			std::size_t size() const;

		private:
			std::vector<Peak> peaks;
		};

		double midpoint(const Peak& peak);

	}
}
#endif