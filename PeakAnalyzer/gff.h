#ifndef BIOSCRIPTS_GFF_H
#define BIOSCRIPTS_GFF_H

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "identifier.h"
#include "peak.h"
#include "range.h"

namespace bioscripts
{
	namespace gff
	{
		struct Record
		{
			using Position = std::size_t;
			using Database = std::string;
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
			Range span;
			//Position start_pos;
			//Position end_pos;
			//std::optional<uint8_t> phase;
			//std::optional<double> score;

			Identifier<Full> sequence_id;
			//Database source;
			std::string attributes;

			auto operator<=>(const Record& rhs) const = default;

			Position start() const
			{
				return span.start;
			}

			Position end() const
			{
				return span.end;
			}

		};

		std::string extractAttribute(const Record& record, const std::string& attribute_name);

		class Records
		{
		public:
			Records(const std::filesystem::path& gff_records);

			using iterator = std::unordered_map<std::string, std::vector<Record>>::iterator;
			using const_iterator = std::unordered_map<std::string, std::vector<Record>>::const_iterator;
			using reference = Record&;
			using const_reference = const Record;
			using pointer = Record*;
			using const_pointer = const Record*;

			iterator begin();
			const_iterator begin() const;
			const_iterator cbegin() const;
			iterator end();
			const_iterator end() const;
			const_iterator cend() const;

			std::vector<Record>& data(const Identifier<Full>& sequence_id);
			const std::vector<Record>& data(const Identifier<Full>& sequence_id) const;
			
			/**
			 * @brief  Find all GFF records that overlap the @a genomic position with the same @sequence id
			 */
			std::vector<Record> getRecordsAt(const std::size_t genomic_position, const Identifier<Full>& sequence_id);

			/**
			 * @brief  Find all GFF records that overlap the @a genomic position with the same @sequence id and @a type
			 */
			std::vector<Record> getRecordsAt(const std::size_t genomic_position, const Identifier<Full>& sequence_id, const Record::Type type);

			/**
			 * @brief  Find the closest record to @a genomic_position with the same gene identifier as @a peak_gene_id; as well as same @a sequence_id and @a type
			 * @return  The closest record, or nullptr if no record could not be 
			 */
			pointer findClosestRecord(std::size_t genomic_position, const Identifier<Full>& sequence_id, const Identifier<Gene>& peak_gene_id, Record::Type type);

			/**
			 * @brief  Return the number of records currently held.
			 */
			std::size_t size() const;


			//iterator find(const std::string& sequence_id, );
		private:
			std::unordered_map<std::string, std::vector<Record>> records;
		};

		Records fetchRecords(Records records, Record::Type type);

		std::vector<bioscripts::gff::Record> collectCodingSequenceRecords(const bioscripts::gff::Record& starting_record, const bioscripts::gff::Records& records);
	}
}

#endif // !BIOSCRIPTS_GFF_H