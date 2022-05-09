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

        //Record& Records::findClosestRecord(std::size_t genomic_position, std::string sequence_id)
        //{
        //    for (const auto& record : records[sequence_id]) {

        //    }
        //}

        Records::pointer Records::findClosestRecord(std::size_t genomic_position, std::string sequence_id, Record::Type type)
        {
            //WARNING: Parentheses around max? Macro expansion
            std::size_t record_distance_to_genomic_position = std::numeric_limits<std::size_t>::max();
            Records::pointer closest_record = nullptr;
            for (const auto& record : records[sequence_id]) {
                if (record.type != type) {
                    continue;
                }

                if (std::abs(record_distance_to_genomic_position - record.start_pos) < record_distance_to_genomic_position) {
                    record_distance_to_genomic_position = std::abs(record_distance_to_genomic_position - record.start_pos);
                    closest_record = &record;
                }
            }
            return closest_record;
        }

        std::vector<Record> Records::findUnderlyingRecords(const std::size_t genomic_position, const std::string& sequence_id)
        {
            if (!records.contains(sequence_id)) {
                return;
            }

            std::vector<Record> results;
            for (const auto& record : records[sequence_id]) {
                if (record.start_pos <= genomic_position && record.end_pos >= genomic_position) {
                    results.push_back(record);
                }
            }

            return records;
        }

        std::vector<Record> Records::findUnderlyingRecords(const std::size_t genomic_position, const std::string& sequence_id, const Record::Type type)
        {
            auto results = findUnderlyingRecords(genomic_position, sequence_id);
            auto IsWrongFeatureType = [&type](const auto& elem) {
                return (elem.type != type);
            };

            for (const auto& record : results) {
                std::erase_if(results, IsWrongFeatureType);

            }

            return results;
        }

        Records::reference Records::findLastRecord(std::string sequence_id, std::string feature_id, Record::Type type)
        {

        }

    }
}