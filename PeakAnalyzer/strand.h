#ifndef BIOSCRIPTS_STRAND_H
#define BIOSCRIPTS_STRAND_H

#include <string>

namespace bioscripts
{

    enum class Strand : uint8_t
    {
        Sense,
        Antisense,
        Unknown
    };

    Strand deduceStrand(const std::string& str);
}

#endif