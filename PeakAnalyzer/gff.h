#ifndef BIOSCRIPTS_GFF_H
#define BIOSCRIPTS_GFF_H

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace bioscripts
{
    enum class Strand : uint8_t
    {
        Sense,
        Antisense,
        Unknown
    };

    Strand deduceStrand(const std::string& str);

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

        class Records
        {
        public:
            Records(const std::filesystem::path& gff_records);

            using iterator = std::vector<Record>::iterator;
            using const_iterator = std::vector<Record>::const_iterator;

            //iterator find(const std::string& sequence_id, );
        private:
            std::unordered_map<std::string, std::vector<Record>> records;
        };
    }
}

#endif // !BIOSCRIPTS_GFF_H