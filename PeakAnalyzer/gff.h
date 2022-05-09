#ifndef BIOSCRIPTS_GFF_H
#define BIOSCRIPTS_GFF_H

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "peak.h"

namespace bioscripts
{
    namespace gff
    {
        struct Record
        {
            enum class Type : uint8_t
            {
                chromosome,
                gene,
                mRNA,
                five_prime_UTR,
                exon,
                CDS,
                three_prime_UTR,
                ncRNA_gene,
                lnc_RNA,
                miRNA,
                tRNA,
                ncRNA,
                snoRNA,
                snRNA,
                rRNA,
                Unknown
            };

            Type type;
            Strand strand;
            std::size_t start_pos;
            std::size_t end_pos;
            std::optional<uint8_t> phase;
            std::optional<double> score;

            std::string sequence_id; 
            std::string source;
            std::string attributes;

            auto operator<=>(const Record& rhs) const = default;
        };

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

        class Records
        {
        public:
            Records(const std::filesystem::path& gff_records);

            using iterator = std::vector<Record>::iterator;
            using const_iterator = std::vector<Record>::const_iterator;
            using reference = Record&;
            using const_reference = const Record;
            using pointer = Record*;
            using const_pointer = const Record*;

            //reference findClosestRecord(std::size_t genomic_position, std::string sequence_id);
            pointer findClosestRecord(std::size_t genomic_position, std::string sequence_id, Record::Type type);

            reference findLastRecord(Identifier sequence_id, std::string feature_id, Record::Type type);


            std::vector<Record> findUnderlyingRecords(const std::size_t genomic_position, const std::string& sequence_id);
            std::vector<Record> findUnderlyingRecords(const std::size_t genomic_position, const std::string& sequence_id, const Record::Type type);


            //iterator find(const std::string& sequence_id, );
        private:
            std::unordered_map<std::string, std::vector<Record>> records;
        };

 
    }

    
}

#endif // !BIOSCRIPTS_GFF_H