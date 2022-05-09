#include "strand.h"

namespace bioscripts
{
    Strand deduceStrand(const std::string& str) {
        if (str == "+") {
            return bioscripts::Strand::Sense;
        }
        else if (str == "-") {
            return bioscripts::Strand::Antisense;
        }
        else {
            return bioscripts::Strand::Unknown;
        }
    }
}