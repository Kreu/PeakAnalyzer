#include "helpers.h"

namespace helper
{
    std::vector<std::string> tokenise(std::string str, char delimiter) {
        std::vector<std::string> tokens;
        std::stringstream s{ str };
        std::string string_token{};
        while (std::getline(s, string_token, delimiter)) {
            if (string_token.length() != 0) {
                tokens.emplace_back(string_token);
            }
        }
        return tokens;
    }
}