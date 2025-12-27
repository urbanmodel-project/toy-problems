#ifndef URBAN_RAD_COMMON_HPP
#define URBAN_RAD_COMMON_HPP

#include <private/DataTypesImpl.h>

namespace URBANXX {
struct RadIndices {
  int c;
  int ib;
  int rtype;
};

struct RadOutput {
  Real Abs, AbsByWt; // absorbed
  Real Ref, RefByWt; // reflected
  Real Emi, EmiByWt; // emitted
};

struct RadRefComponents {
  Real ToSky, ToSkyByWt;
  Real ToRoad, ToRoadByWt;
  Real ToSunwall, ToSunwallByWt;
  Real ToShadewall, ToShadewallByWt;
  Real ToOtherwall, ToOtherwallByWt;
};

} // namespace URBANXX
#endif