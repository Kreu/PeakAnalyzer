#include "identifier.h"

namespace bioscripts
{
	bool operator==(const Identifier<Gene>& lhs, const Identifier<Transcript>& rhs)
	{
		return lhs.gene() == rhs.gene();
	}

	bool operator==(const Identifier<Transcript>& lhs, const Identifier<Gene>& rhs)
	{
		return rhs == lhs;
	}

	bool operator==(const Identifier<Gene>& lhs, const Identifier<Gene>& rhs)
	{
		return lhs.to_string() == rhs.to_string();
	}

	bool operator==(const Identifier<Transcript>& lhs, const Identifier<Transcript>& rhs)
	{
		return lhs.to_string() == rhs.to_string();
	}

}