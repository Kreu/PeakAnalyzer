#ifndef BIOSCRIPTS_IDENTIFIER_H
#define BIOSCRIPTS_IDENTIFIER_H

#include "helpers.h"

#include <string>

namespace bioscripts
{
	struct Full;
	struct Transcript;
	struct Gene;

	template <typename Type>
	class Identifier {};

	template <>
	class Identifier<Full>
	{
	public:
		Identifier(const std::string& id) : identifier{ id } {};

		std::string to_string() const
		{
			return identifier;
		}
	private:
		std::string identifier;
	};

	template <>
	class Identifier<Gene>
	{
	public:
		Identifier(const std::string& id)
		{
			auto tokens = helper::tokenise(id, '.');

			if (tokens.empty()) {
				return;
			}
			gene_name = tokens[0];
		}

		std::string gene() const
		{
			return gene_name;
		}

		std::string to_string() const
		{
			return gene_name;
		}
	private:
		std::string gene_name;
	};

	template <>
	class Identifier<Transcript>
	{
	public:
		Identifier(const std::string& id)
		{
			auto tokens = helper::tokenise(id, '.');

			if (tokens.size() < 2) {
				throw std::logic_error("Cannot create a transcript identifier from " + id);
			}
			gene_name = tokens[0];
			transcript_version = tokens[1];
		}

		std::string gene() const
		{
			return gene_name;
		}

		std::string version() const
		{
			return transcript_version;
		}

		std::string to_string() const
		{
			static constexpr auto delimiter = '.';
			return gene_name + delimiter + transcript_version;
		}
	private:
		std::string gene_name;
		std::string transcript_version;
	};

	bool operator==(const Identifier<Full>& lhs, const Identifier<Full>& rhs);

	bool operator==(const Identifier<Gene>& lhs, const Identifier<Transcript>& rhs);
	bool operator==(const Identifier<Transcript>& lhs, const Identifier<Gene>& rhs);

	bool operator==(const Identifier<Gene>& lhs, const Identifier<Gene>& rhs);
	bool operator!=(const Identifier<Gene>& lhs, const Identifier<Gene>& rhs);

	bool operator==(const Identifier<Transcript>& lhs, const Identifier<Transcript>& rhs);
	bool operator!=(const Identifier<Transcript>& lhs, const Identifier<Transcript>& rhs);

}

#endif