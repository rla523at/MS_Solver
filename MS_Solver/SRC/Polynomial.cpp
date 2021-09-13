#include "../INC/Polynomial.h"

namespace ms {


	bool is_positive_odd_number(const double val) {
		if (val < 0)
			return false;

		if (val - std::floor(val) == 0)
			return static_cast<size_t>(val) % 2 == 0;
		else
			return false;
	}

	bool is_natural_number(const double val) {
		if (val < 0)
			return false;

		if (val - std::floor(val) == 0)
			return true;
		else
			return false;
	}

}

