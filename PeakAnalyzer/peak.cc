#include <fstream>

#include "helpers.h"
#include "peak.h"

#include "easylogging++.h"

namespace bioscripts
{
	namespace peak
	{
		Peaks::Peaks(const std::filesystem::path& peak_file)
		{
            std::ifstream f{ peak_file };
            if (!f.is_open()) {
                LOG(ERROR) << "Could not open peak file";
                return;
            }

            std::string line;
            std::getline(f, line); /* Skips the header in database file */
            while (std::getline(f, line)) {
                const auto tokens = helper::tokenise(line, '\t');

                if (tokens.size() != 21) {
                    continue;
                }

                std::string sequence_identifier = tokens[1];
                std::size_t start_pos = std::stoull(tokens[2]);
                std::size_t end_pos = std::stoull(tokens[3]);

                //Add one to end-pos because Range is 0-based [start, end)
                //but GFF coordinates are [start, end]
                auto peak_range = Range{ start_pos, end_pos + 1 };

                //uint16_t width = std::stoul(tokens[4]);
                //uint16_t score = std::stoul(tokens[6]);
                //double pvalue = std::stod(tokens[8]);
                //double qvalue = std::stod(tokens[9]);

                std::string feature_identifier = tokens[11];
                //std::size_t feature_start_pos = std::stoull(tokens[12]);
                //std::size_t feature_end_pos = std::stoull(tokens[13]);
                //std::string feature_ensembl_id = tokens[19];
                auto strand = bioscripts::deduceStrand(tokens[14]);

                peaks.push_back(Peak{
                    .span = peak_range,
                    //.start_pos = start_pos,
                    //.end_pos = end_pos,
                    //.width = width,
                    //.score = score,
                    //.pvalue = pvalue,
                    //.qvalue = qvalue,
                    .sequence_id = sequence_identifier,
                    .feature = Feature {
                        //.start_pos = feature_start_pos,
                        //.end_pos = feature_end_pos,
                        .strand = strand,
                        .identifier = feature_identifier,
                        //.ensembl_id = "",
                        }
                    });

            }
		}

        std::size_t Peaks::size() const
        {
            return peaks.size();
        }

        double midpoint(const Peak& peak)
        {
            return (peak.start() + peak.end()) / 2.0;
        }

        Peaks::iterator Peaks::begin()
        {
            return std::begin(peaks);
        }

        Peaks::const_iterator Peaks::begin() const
        {
            return std::begin(peaks);
        }

        Peaks::const_iterator Peaks::cbegin() const
        {
            return std::cbegin(peaks);
        }

        Peaks::iterator Peaks::end()
        {
            return std::end(peaks);
        }

        Peaks::const_iterator Peaks::end() const
        {
            return std::end(peaks);
        }

        Peaks::const_iterator Peaks::cend() const
        {
            return std::end(peaks);
        }

	}
}