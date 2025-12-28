#ifndef URBAN_RAD_COMMON_IMPL_H
#define URBAN_RAD_COMMON_IMPL_H

#include <private/UrbanDataTypesImpl.h>

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
#endif // URBAN_RAD_COMMON_IMPL_H