#ifndef URBAN_RAD_COMMON_HPP
#define URBAN_RAD_COMMON_HPP

namespace URBANXX {
struct RadIndices {
  int c;
  int ib;
  int rtype;
};

struct RadOutput {
  Real Abs;
  Real Ref;
};

struct RadRefComponents {
  Real ToSky;
  Real ToRoad;
  Real ToSunwall;
  Real ToShadewall;
  Real ToOtherwall;
};

}
#endif