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
                std::string feature_identifier = tokens[11];

                auto strand = bioscripts::deduceStrand(tokens[14]);

                peaks.push_back(Peak{
                    .span = peak_range,
                    .strand = strand,
                    .associated_identifier = feature_identifier,
                    .sequence_id = sequence_identifier,
                });
            }
		}

        std::size_t Peaks::size() const
        {
            return peaks.size();
        }

        double midpoint(const Peak& peak) 
        {
            return (centre(peak.span));
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