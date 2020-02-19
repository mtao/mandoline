#include "mandoline/vertex.hpp"

#if defined(_MSC_VER) || defined(__APPLE__)
const static double threshold_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
#endif
