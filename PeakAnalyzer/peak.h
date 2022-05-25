#ifndef BIOSCRIPTS_PEAK_H
#define BIOSCRIPTS_PEAK_H

#include <cstdint>
#include <filesystem>
#include <string>

#include "identifier.h"
#include "strand.h"
#include "range.h"

namespace bioscripts
{
	namespace peak
	{
		struct Peak
		{
			Range span;
			Strand strand;
			Identifier<Gene> associated_identifier; //Denotes the identifier under the peak. Usually it is the gene name, e.g. AT1G12345
			//Identifier<Transcript> ensembl_id; //This is the identifier associate with the peak after R analysis. More detailed "identifier", e.g. AT1G12345.2.exon4
			std::string sequence_id; //Denotes the ID of the sequence the peak is found on. Usually it is the chromosome identifier.

			Position start() const
			{
				return span.start;
			}

			Position end() const
			{
				return span.end;
			}
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