#include "identifier.h"
#include "helpers.h"

namespace bioscripts
{
	Identifier::Identifier(const std::string& input_identifier)
	{
		auto tokens = helper::tokenise(input_identifier, '.');

		if (tokens.empty()) {
			return;
		}
		else if (tokens.size() == 1) {
			identifier = tokens[0];
			is_version_specified = false;
		}
		else {
			identifier = tokens[0];
			version = tokens[1];
			is_version_specified = true;
		}
	}

	std::string Identifier::to_string() const
	{
		if (is_version_specified) {
			return identifier + "." + version;
		}
		else {
			return identifier;
		}
	}
}