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
			identifier = tokens[0];
		}

		//Identifier(const Identifier<Transcript>& id) : identifier(id.gene)

		std::string to_string() const
		{
			return identifier;
		}
	private:
		std::string identifier;
	};

	template <>
	class Identifier<Transcript>
	{
	public:
		Identifier(const std::string& id)
		{
			auto tokens = helper::tokenise(id, '.');

			if (tokens.size() < 2) {
				return;
			}
			identifier = tokens[0];
			version = tokens[1];
		}

		std::string to_string() const
		{
			static constexpr auto delimiter = '.';
			return identifier + delimiter + version;
		}
	private:
		std::string identifier;
		std::string version;
	};

	//template <typename T>
	//bool operator==(const T& lhs, const T& rhs)
	//{
	//	return rhs.to_string() == lhs.to_string();
	//}


	//bool operator==(const Identifier<Transcript>& rhs, const Identifier<Transcript>& lhs)
	//{
	//	return rhs.to_string() == lhs.to_string();
	//}



}

#endif