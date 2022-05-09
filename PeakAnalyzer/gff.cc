#include <fstream>

#include "helpers.h"
#include "gff.h"

namespace {
    bioscripts::gff::Record::Type deduceType(const std::string& str)
    {
        if (str == "chromosome") {
            return bioscripts::gff::Record::Type::chromosome;
        }
        else if (str == "mRNA") {
            return bioscripts::gff::Record::Type::mRNA;
        }
        else if (str == "five_prime_UTR") {
            return bioscripts::gff::Record::Type::five_prime_UTR;
        }
        else if (str == "exon") {
            return bioscripts::gff::Record::Type::exon;
        }
        else if (str == "CDS") {
            return bioscripts::gff::Record::Type::CDS;
        }
        else if (str == "ncRNA_gene") {
            return bioscripts::gff::Record::Type::ncRNA_gene;
        }
        else if (str == "lnc_RNA") {
            return bioscripts::gff::Record::Type::lnc_RNA;
        }
        else if (str == "miRNA") {
            return bioscripts::gff::Record::Type::miRNA;
        }
        else if (str == "tRNA") {
            return bioscripts::gff::Record::Type::tRNA;
        }
        else if (str == "ncRNA") {
            return bioscripts::gff::Record::Type::ncRNA;
        }
        else if (str == "snoRNA") {
            return bioscripts::gff::Record::Type::snoRNA;
        }
        else if (str == "snRNA") {
            return bioscripts::gff::Record::Type::snRNA;
        }
        else if (str == "rRNA") {
            return bioscripts::gff::Record::Type::rRNA;
        }
        else {
            return bioscripts::gff::Record::Type::Unknown;
        }
    }
}

namespace bioscripts
{
    namespace gff
    {
        Records::Records(const std::filesystem::path& gff_records)
        {
            std::ifstream f{ gff_records };
            if (!f.is_open()) {
                return;
            }

            std::string line;
            std::getline(f, line); /* Skips the header in database file */
            while (std::getline(f, line)) {
                const auto tokens = helper::tokenise(line, '\t');

                if (tokens.size() != 9) {
                    continue;
                }

                const auto& sequence_id = tokens[0];
                const auto& source = tokens[1];

                const auto type = deduceType(tokens[2]);
                if (type == bioscripts::gff::Record::Type::Unknown) {
                    return;
                }

                std::size_t start_pos = std::stoull(tokens[3]);
                std::size_t end_pos = std::stoull(tokens[4]);

                std::optional<double> score;
                if (auto score_string = tokens[5]; score_string == ".") {
                    score = std::nullopt;
                }
                else {
                    score = std::stod(score_string);
                }

                const auto strand = deduceStrand(tokens[6]);
                if (strand == bioscripts::Strand::Unknown) {
                    return;
                }

                std::optional<uint8_t> phase;
                if (auto phase_string = tokens[7]; phase_string == ".") {
                    phase = std::nullopt;
                }
                else {
                    phase = std::stod(phase_string);
                }

                const auto& attributes = tokens[8];

                this->records[sequence_id].push_back(Record{
                    .type = type,
                    .strand = strand,
                    .start_pos = start_pos,
                    .end_pos = end_pos,
                    .phase = phase,
                    .score = score,
                    .sequence_id = sequence_id,
                    .source = source,
                    .attributes = attributes
                    });
            }

            auto Comparator = [](const auto& first, const auto& second) {
                return (first.sequence_id < second.sequence_id)
                    && (first.start_pos < second.start_pos);
            };

            for (auto& [chromosome_id, entries] : records) {
                std::sort(std::begin(entries), std::end(entries), Comparator);
            }
        }

        Records::reference Records::findClosestRecord(const bioscripts::peak::Peak& peak)
        {
            auto peak_sequence_identifier = peak.sequence_id;

            if (!records.contains(peak_sequence_identifier)) {
                return;
            }

            auto& peak_chromosome_records = records.at(peak_sequence_identifier);






            //auto isEndPosLessThanPeakPos = [](const auto& e1, const std::size_t& peak_pos)
            //{
            //    return (e1.end_pos < peak_pos);
            //};

            //auto isStartPosLessThanPeakPos = [](const auto& e1, const std::size_t& peak_pos)
            //{
            //    return (e1.start_pos < peak_pos);
            //};

            //const auto peak_midpoint = bioscripts::peak::midpoint(peak);
            //auto first_overlapping_gff_record = std::lower_bound(std::begin(peak_chromosome_records), std::end(peak_chromosome_records), peak_midpoint, isEndPosLessThanPeakPos);
            //auto first_non_overlapping_gff_record = std::upper_bound(std::begin(peak_chromosome_records), std::end(peak_chromosome_records), peak_midpoint, isStartPosLessThanPeakPos);

            //if (first_overlapping_gff_record == std::end(peak_chromosome_records)) {
            //    return;
            //}

            //if (first_non_overlapping_gff_record == std::end(peak_chromosome_records)) {
            //    return;
            //}


            
            return peak_chromosome_records[0];



        }
    }
}