#ifndef BIOSCRIPTS_IDENTIFIER_H
#define BIOSCRIPTS_IDENTIFIER_H

#include <string>

namespace bioscripts
{
	class Identifier
	{
	public:
		Identifier(const std::string& identifier);

		std::string to_string() const;
	private:
		std::string identifier;
		std::string version;
		bool is_version_specified;
	};
}

#endif