#ifndef URBAN_RAD_COMMON_HPP
#define URBAN_RAD_COMMON_HPP

#include <DataTypes.hpp>

namespace URBANXX {
struct RadIndices {
  int c;
  int ib;
  int rtype;
};

struct RadOutput {
  Real Abs; // absorbed
  Real Ref; // reflected
  Real Emi; // emitted
};

struct RadRefComponents {
  Real ToSky;
  Real ToRoad;
  Real ToSunwall;
  Real ToShadewall;
  Real ToOtherwall;
};

} // namespace URBANXX
#endif