#include <fstream>
#include <cmath>
#include <cassert>


#include <iostream>

#include "identifier.h"
#include "helpers.h"
#include "gff.h"

#include "easylogging++.h"

namespace {
	/**
	 * @brief  Turn a string representation of a sequence type into an actual type.
	 */
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
	 * 
	 * The distance is calculated from the end of the record to the genomic position if the
	 * record appears before the @a genomic_position; otherwise the distance is calculated
	 * to the genomic position until the start position of the record.
	 */
	std::size_t distanceToRecord(const std::size_t genomic_position, const bioscripts::gff::Record& record)
	{
		std::size_t distance_to_record = 0;
		if (record.end_pos < genomic_position) {
			distance_to_record = genomic_position - record.end_pos;
		}
		else if (record.start_pos > genomic_position) {
			distance_to_record = record.start_pos - genomic_position;
		}
		return distance_to_record;
	}
}

namespace bioscripts
{
	namespace gff
	{
		Records::Records(const std::filesystem::path& gff_records)
		{
			LOG(DEBUG) << "Parsing GFF records from " << gff_records.string();
			std::ifstream f{ gff_records };
			if (!f.is_open()) {
				LOG(ERROR) << "Failed to open " << gff_records.string();
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

		std::size_t Records::size() const
		{
			std::size_t records_number = 0;
			for (const auto& [sequence_id, entries] : records) {
				records_number += entries.size();
			}
			return records_number;
		}

		/**
		 * @brief  Extract the attribute value of the given @a attribute_name
		 * @return  Value of the given @a attribute_name associated with the @record,
		 *			or an empty string if no @a attribute_name is present.
		 */
		std::string extractAttribute(const Record& record, const std::string& attribute_name)
		{
			//LOG(DEBUG) << "Extracting " << attribute_name << " from record with sequence ID of " << record.sequence_id.to_string() << " and attributes of " << record.attributes << "\n";
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
			//LOG(DEBUG) << "Extracted attribute value is \"" << record.attributes.substr(attribute_value_start_pos, substring_length) << "\"\n";
			return record.attributes.substr(attribute_value_start_pos, substring_length);
		}

		Records::pointer Records::findClosestRecord(std::size_t genomic_position, const Identifier<Full>& sequence_id, const Identifier<Gene>& peak_gene_id, Record::Type type)
		{
			/* Closest record is defined as follows :
				a) Record whose end position is closest to the genomic_position if that end position < genomic_position OR
				b) Record whose start position is closest to the genomic_position if that start position > genomic_position
			 */
			//LOG(DEBUG) << "Finding closest record to " << genomic_position << " on sequence \"" << sequence_id.to_string() << "\", peak gene identifier is " << peak_gene_id.to_string() << "\n";
			std::int64_t current_smallest_distance = (std::numeric_limits<std::int64_t>::max)();
			Records::pointer closest_record = nullptr;
			for (auto& record : records[sequence_id.to_string()]) {
				if (record.type != type) {
					continue;
				}

				auto record_id = Identifier<Transcript>{ extractAttribute(record, "ID=CDS") };
				if (record_id != peak_gene_id) {
					continue;
				}

				auto distance_to_record = distanceToRecord(genomic_position, record);
				//LOG(DEBUG) << "Calculated distance to record: " << distance_to_record << "\n";
				if (distance_to_record < current_smallest_distance) {
					current_smallest_distance = distance_to_record;
					closest_record = &record;
					//LOG(DEBUG) << "Found a closer record on sequence \"" << closest_record->sequence_id.to_string() << "\", record attributes are " << closest_record->attributes << "\n";
				}
			}
			return closest_record;
		}

		std::vector<Record> Records::getRecordsAt(const std::size_t genomic_position, const Identifier<Full>& sequence_id)
		{
			//LOG(DEBUG) << "Finding underlying record at position " << genomic_position << " on sequence \"" << sequence_id.to_string() << "\n";
			if (!records.contains(sequence_id.to_string())) {
				return {};
			}

			std::vector<Record> results;
			for (const auto& record : records[sequence_id.to_string()]) {
				if (record.start_pos <= genomic_position && record.end_pos >= genomic_position) {
					//LOG(DEBUG) << "Found a record, attributes are " << record.attributes << "\n";
					results.push_back(record);
				}
			}
			//LOG(DEBUG) << "Found " << results.size() << " underlying records\n";
			return results;
		}

		std::vector<Record> Records::getRecordsAt(const std::size_t genomic_position, const Identifier<Full>& sequence_id, const Record::Type type)
		{
			auto results = getRecordsAt(genomic_position, sequence_id);

			auto IsWrongFeatureType = [&type](const auto& elem) {
				return (elem.type != type);
			};
			//LOG(DEBUG) << "Found " << results.size() << " underlying records\n";
			for (const auto& record : results) {
				std::erase_if(results, IsWrongFeatureType);
			}
			//LOG(DEBUG) << "After removing records of incorrect type, there are " << results.size() << " underlying records\n";
			return results;
		}

		Records::iterator Records::begin()
		{
			return records.begin();
		}

		Records::const_iterator Records::begin() const
		{
			return records.begin();
		}

		Records::const_iterator Records::cbegin() const
		{
			return records.cbegin();
		}

		Records::iterator Records::end()
		{
			return records.end();
		}

		Records::const_iterator Records::end() const
		{
			return records.end();
		}

		Records::const_iterator Records::cend() const
		{
			return records.cend();
		}

		std::vector<Record>& Records::data(const Identifier<Full>& sequence_id)
		{
			auto& corresponding_records = records.at(sequence_id.to_string());
			return corresponding_records;
		}

		const std::vector<Record>& Records::data(const Identifier<Full>& sequence_id) const
		{
			auto& corresponding_records = records.at(sequence_id.to_string());
			return corresponding_records;
		}


		Records fetchRecords(Records records, Record::Type type)
		{
			auto isWrongRecordType = [&type](const auto& elem)
			{
				return elem.type != type;
			};

			for (auto& [id, records] : records) {
				std::erase_if(records, isWrongRecordType);
			}

			return records;
		}

		/**
		 * @brief  Find all CDS type GFF records that belong to the same transcript as @a starting_record. Only records
		 *		   after the @a starting_record are considered.
		 *
		 * @return  All CDS GFF records corresponding to the same transcript ID of the @starting_record, including the @a starting_record. 
		 */
		std::vector<bioscripts::gff::Record> collectCodingSequenceRecords(const bioscripts::gff::Record& starting_record, const bioscripts::gff::Records& records)
		{
			auto record_sequence = starting_record.sequence_id.to_string();
			//LOG(DEBUG) << "Collecting all CDS records of " << record_sequence << "\n";
			auto hasWrongSequenceType = [&starting_record](const auto& record)
			{
				return starting_record.type != record.type;
			};

			auto same_chromosome_records = records.data(record_sequence);
			std::erase_if(same_chromosome_records, hasWrongSequenceType);

			//If the record is on the sense strand, the subsequent CDS records are upstream of the starting record.
			//If the record is on the antisense strand, the subsequent CDS records are before.
			std::vector<bioscripts::gff::Record>::iterator end_iterator;
			if (starting_record.strand == bioscripts::Strand::Sense) {
				end_iterator = std::end(same_chromosome_records);
			}
			else if (starting_record.strand == bioscripts::Strand::Antisense) {
				end_iterator = std::begin(same_chromosome_records);
			}
			else {
				return {};
			}

			auto Comparator = [](const auto& gff_record, const auto& start_pos) {
				return (gff_record.start_pos < start_pos);
			};
			auto start_looking_from = std::lower_bound(std::begin(same_chromosome_records), std::end(same_chromosome_records), starting_record.start_pos, Comparator);
			const auto starting_record_id = bioscripts::Identifier<bioscripts::Transcript>{ bioscripts::gff::extractAttribute(starting_record, "ID=CDS") };



			std::vector<bioscripts::gff::Record> final_records;
			final_records.push_back(starting_record);

			//Find the first record with a starting position larger than the current starting_record.
			for (auto it = start_looking_from; it != end_iterator;) {
				const auto& record = *it;
				if (record.start_pos <= starting_record.start_pos) {
					continue;
				}

				const auto record_transcript_id = bioscripts::Identifier<bioscripts::Transcript>{ bioscripts::gff::extractAttribute(record, "ID=CDS") };
				if (record_transcript_id.gene() != starting_record_id.gene()) {
					//std::cout << "Record transcript does not match starting record\n";
					break;
				}
				//LOG(DEBUG) << "Found a record corresponding to the CDS with sequence id \"" << record.sequence_id.to_string() << "\" with attributes " << record.attributes << "\n";
				final_records.push_back(record);

				if (starting_record.strand == bioscripts::Strand::Sense) {
					++it;
				}
				else if (starting_record.strand == bioscripts::Strand::Antisense) {
					--it;
				}
			}
			return final_records;
		}

	}
}