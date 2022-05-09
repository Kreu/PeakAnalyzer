#ifndef BIOSCRIPTS_PEAK_H
#define BIOSCRIPTS_PEAK_H

#include <cstdint>
#include <filesystem>
#include <string>

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
			std::string identifier;
			std::string ensembl_id;
		};

		struct Peak
		{
			std::size_t start_pos;
			std::size_t end_pos;
			uint16_t width;
			uint16_t score;
			double pvalue;
			double qvalue;
			std::string sequence_id;
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
		private:
			std::vector<Peak> peaks;
		};

		double midpoint(const Peak& peak);

	}
}
#endif