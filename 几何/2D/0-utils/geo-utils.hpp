#pragma once
#include "aizalib.h"

namespace Geo2D {

const ld PI = acos(-1.0);

inline int sign(ld x) { return x < -EPS ? -1 : x > EPS; }
inline int dcmp(ld x, ld y) { return sign(x - y); }

} // namespace Geo2D