#include <fstream>
#include <cmath>
#include <cassert>


#include <iostream>

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
        else if (str == "gene") {
            return bioscripts::gff::Record::Type::gene;
        }
        else if (str == "five_prime_UTR") {
            return bioscripts::gff::Record::Type::five_prime_UTR;
        }
        else if (str == "three_prime_UTR") {
            return bioscripts::gff::Record::Type::three_prime_UTR;
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

    /**
     * @brief  Calculate the absolute distance between record and @a genomic_position 
     */
    std::size_t distanceToRecord(const std::size_t genomic_position, const bioscripts::gff::Record& record)
    {
        auto distance_to_record = 0;
        if (record.end_pos < genomic_position) {
            distance_to_record = genomic_position - record.end_pos;
        }
        else if (record.start_pos > genomic_position) {
            distance_to_record = record.start_pos - genomic_position;
        }
        assert(distance_to_record > 0, "Distance to record can not negative");
        return distance_to_record;
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
                static constexpr auto comment_token = '#';
                if (line.starts_with(comment_token)) {
                    continue;
                }

                const auto tokens = helper::tokenise(line, '\t');
                if (tokens.size() != 9) {
                    continue;
                }

                const auto& sequence_id = tokens[0];
                const auto& source = tokens[1];

                const auto type = deduceType(tokens[2]);
                if (type == bioscripts::gff::Record::Type::Unknown) {
                    continue;
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
                    continue;
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
                return (first.sequence_id.to_string() < second.sequence_id.to_string())
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

        std::size_t Records::size() const
        {
            return records.size();
        }

        std::string extractAttribute(const Record& record, const std::string& attribute_name)
        {
            auto attribute_name_start_pos = record.attributes.find(attribute_name);
            if (attribute_name_start_pos == std::string::npos) {
                return "";
            }

            //The + 1 is for the '=' character that is always preceded by the attribute value
            //e.g. it is in the format of "attribute_name=attribute_value;"
            auto attribute_value_start_pos = attribute_name_start_pos + attribute_name.length() + 1;

            constexpr auto delimiter = ';';
            auto next_delimiter_pos = record.attributes.find(delimiter, attribute_value_start_pos);
            auto substring_length = next_delimiter_pos - attribute_value_start_pos;
            return record.attributes.substr(attribute_value_start_pos, substring_length);
        }

        Records::pointer Records::findClosestRecord(std::size_t genomic_position, const Identifier<Full>& sequence_id, const Identifier<Gene>& gene_id, Record::Type type)
        {
            /* Closest record is defined as follows :
                a) Record whose end position is closest to the genomic_position if that end position < genomic_position OR
                b) Record whose start position is closest to the genomic_position if that start position > genomic_position
             */
            //std::cout << "Finding the closest record to " << gene_id.to_string() << " at position " << genomic_position << "\n";
            std::int64_t current_smallest_distance = std::numeric_limits<std::int64_t>::max();
            Records::pointer closest_record = nullptr;
            for (auto& record : records[sequence_id.to_string()]) {
                if (record.type != type) {
                    //std::cout << "Record has the wrong type\n";
                    continue;
                }

                auto record_id = Identifier<Transcript>{ extractAttribute(record, "ID=CDS") };
                if (record_id != gene_id) {
                    //std::cout << "Record ID of " << record_id.to_string() << " does not match " << gene_id.to_string() << "\n";
                    continue;
                }                

                auto distance_to_record = distanceToRecord(genomic_position, record);
                //std::cout << "Distance to current record: " << current_smallest_distance << "\n";
                if (distance_to_record < current_smallest_distance) {
                    //std::cout << "Found a closer record at a distance of " << distance_to_record << "\n";
                    current_smallest_distance = distance_to_record;
                    closest_record = &record;
                }
            }
            return closest_record;
        }

        std::vector<Record> Records::findUnderlyingRecords(const std::size_t genomic_position, const Identifier<Full>& sequence_id)
        {
            if (!records.contains(sequence_id.to_string())) {
                return {};
            }

            std::vector<Record> results;
            for (const auto& record : records[sequence_id.to_string()]) {
                if (record.start_pos <= genomic_position && record.end_pos >= genomic_position) {
                    results.push_back(record);
                }
            }
            return results;
        }

        std::vector<Record> Records::findUnderlyingRecords(const std::size_t genomic_position, const Identifier<Full>& sequence_id, const Record::Type type)
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

        void Records::findLastRecord(const Identifier<Full>& sequence_id, const Identifier<Gene>& feature_id, const Identifier<Gene>& gene_id, Record::Type type)
        {

        }

    }
}