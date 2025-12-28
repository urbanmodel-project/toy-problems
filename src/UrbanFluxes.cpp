#include <Kokkos_Core.hpp>
#include <cmath>
#include <iostream>
#include <private/UrbanDataImpl.h>
#include <private/UrbanFluxesImpl.h>

#define ALLOCATE_VIEW(viewname, type, ...)                                     \
  viewname = type(#viewname, __VA_ARGS__);

#define PERVIOUS_ROAD_FRACTION 0.16666667163372040
#define ROOF_FRACTION 0.69999998807907104

namespace URBANXX {

UrbanSurface ::UrbanSurface(int N_LUN, Real temperature) {

  ALLOCATE_VIEW(Temp, Array1DR8, N_LUN);
  ALLOCATE_VIEW(es, Array1DR8, N_LUN);
  ALLOCATE_VIEW(esdT, Array1DR8, N_LUN);
  ALLOCATE_VIEW(qs, Array1DR8, N_LUN);
  ALLOCATE_VIEW(qsdT, Array1DR8, N_LUN);

  HostArray1DR8 TempH;
  ALLOCATE_VIEW(TempH, HostArray1DR8, N_LUN);
  for (int i = 0; i < N_LUN; ++i) {
    TempH(i) = temperature;
  }
  Kokkos::deep_copy(Temp, TempH);
}

UrbanSurfaceFluxes::UrbanSurfaceFluxes(UrbanSharedDataBundle &bundle)
    : data_bundle(bundle), Roof(bundle.N_LUN, 292.0),
      SunlitWall(bundle.N_LUN, 292.0), ShadeWall(bundle.N_LUN, 292.0),
      ImperviousRoad(bundle.N_LUN, 274.0), PerviousRoad(bundle.N_LUN, 274.0) {
  ForcTemp = bundle.input.ForcTemp;
  ForcTh = bundle.input.ForcPotTemp;
  ForcRho = bundle.input.ForcRho;
  ForcQ = bundle.input.ForcSpcHumd;
  ForcPbot = bundle.input.ForcPress;
  ForcU = bundle.input.ForcWindU;
  ForcV = bundle.input.ForcWindV;

  Hwr = bundle.geometry.CanyonHwr;

  int N_LUN = data_bundle.N_LUN;
  ALLOCATE_VIEW(Taf, Array1DR8, N_LUN);
  ALLOCATE_VIEW(Qaf, Array1DR8, N_LUN);
  ALLOCATE_VIEW(fPervRoad, Array1DR8, N_LUN);
  ALLOCATE_VIEW(fRoof, Array1DR8, N_LUN);

  Array1DR8 qafH, tafH, fPervRoadH, fRoofH;
  ALLOCATE_VIEW(qafH, HostArray1DR8, N_LUN);
  ALLOCATE_VIEW(tafH, HostArray1DR8, N_LUN);
  ALLOCATE_VIEW(fPervRoadH, HostArray1DR8, N_LUN);
  ALLOCATE_VIEW(fRoofH, HostArray1DR8, N_LUN);

  Real TAF = 283.0;
  Real QAF = 1e-4;
  for (int i = 0; i < N_LUN; ++i) {
    tafH(i) = TAF;
    qafH(i) = QAF;
    fPervRoadH(i) = PERVIOUS_ROAD_FRACTION;
    fRoofH(i) = ROOF_FRACTION;
  }

  Kokkos::deep_copy(Taf, tafH);
  Kokkos::deep_copy(Qaf, qafH);
  Kokkos::deep_copy(fPervRoad, fPervRoadH);
  Kokkos::deep_copy(fRoof, fRoofH);
}

#define GRAVITY 9.80616
#define VKC 0.4
#define ZETAM 1.574 // transition point of flux-gradient relation (wind profile)
#define ZETAT                                                                  \
  0.465 // transition point of flux-gradient relation (temperature profile)
#define CPAIR 1004.64 // specific heat of dry air [J/kg/K]
#define SHR_CONST_TKFRZ 273.15

KOKKOS_INLINE_FUNCTION
void QSat(Real T, Real p, Real &es, Real &esdT, Real &qs, Real &qsdT) {
  // For water vapor (temperature range 0C-100C)
  Real a0 = 6.11213476;
  Real a1 = 0.444007856;
  Real a2 = 0.143064234e-01;
  Real a3 = 0.264461437e-03;
  Real a4 = 0.305903558e-05;
  Real a5 = 0.196237241e-07;
  Real a6 = 0.892344772e-10;
  Real a7 = -0.373208410e-12;
  Real a8 = 0.209339997e-15;
  // For derivative:water vapor
  Real b0 = 0.444017302;
  Real b1 = 0.286064092e-01;
  Real b2 = 0.794683137e-03;
  Real b3 = 0.121211669e-04;
  Real b4 = 0.103354611e-06;
  Real b5 = 0.404125005e-09;
  Real b6 = -0.788037859e-12;
  Real b7 = -0.114596802e-13;
  Real b8 = 0.381294516e-16;
  // For ice (temperature range -75C-0C)
  Real c0 = 6.11123516;
  Real c1 = 0.503109514;
  Real c2 = 0.188369801e-01;
  Real c3 = 0.420547422e-03;
  Real c4 = 0.614396778e-05;
  Real c5 = 0.602780717e-07;
  Real c6 = 0.387940929e-09;
  Real c7 = 0.149436277e-11;
  Real c8 = 0.262655803e-14;
  // For derivative:ice
  Real d0 = 0.503277922;
  Real d1 = 0.377289173e-01;
  Real d2 = 0.126801703e-02;
  Real d3 = 0.249468427e-04;
  Real d4 = 0.313703411e-06;
  Real d5 = 0.257180651e-08;
  Real d6 = 0.133268878e-10;
  Real d7 = 0.394116744e-13;
  Real d8 = 0.498070196e-16;

  Real T_limit = T - SHR_CONST_TKFRZ;
  if (T_limit > 100.0)
    T_limit = 100.0;
  if (T_limit < -75.0)
    T_limit = -75.0;

  Real td = T_limit;

  if (td >= 0.0) {
    es =
        a0 +
        td * (a1 +
              td * (a2 +
                    td * (a3 +
                          td * (a4 +
                                td * (a5 + td * (a6 + td * (a7 + td * a8)))))));

    esdT =
        b0 +
        td * (b1 +
              td * (b2 +
                    td * (b3 +
                          td * (b4 +
                                td * (b5 + td * (b6 + td * (b7 + td * b8)))))));

  } else {
    es =
        c0 +
        td * (c1 +
              td * (c2 +
                    td * (c3 +
                          td * (c4 +
                                td * (c5 + td * (c6 + td * (c7 + td * c8)))))));

    esdT =
        d0 +
        td * (d1 +
              td * (d2 +
                    td * (d3 +
                          td * (d4 +
                                td * (d5 + td * (d6 + td * (d7 + td * d8)))))));
  }

  es = es * 100.0;     // pa
  esdT = esdT * 100.0; // pa/K

  Real vp = 1.0 / (p - 0.378 * es);
  Real vp1 = 0.622 * vp;
  Real vp2 = vp1 * vp;

  qs = es * vp1;         // kg/kg
  qsdT = esdT * vp2 * p; // 1 / K
}

KOKKOS_INLINE_FUNCTION
void UrbanSurface::ComputeQsat(int c, Real p) {
  Real T = Temp(c);
  QSat(T, p, es(c), esdT(c), qs(c), qsdT(c));
}

KOKKOS_INLINE_FUNCTION
void MoninObukIni(Real ur, Real thv, Real dthv, Real zldis, Real z0m, Real &um,
                  Real &obu) {

  const Real ustar = 0.06;
  const Real wc = 0.5;

  if (dthv > 0) {
    um = std::max(ur, 1.0);
  } else {
    um = std::pow(std::pow(ur, 2.0) + std::pow(wc, 2.0), 0.5);
  }

  const Real rib = GRAVITY * zldis * dthv / (thv * std::pow(um, 2.0));

  Real zeta;
  if (rib >= 0) {
    zeta = rib * std::log(zldis / z0m) / (1.0 - 0.5 * std::min(rib, 0.19));
    zeta = std::min(2.0, std::max(zeta, 0.01));
  } else {
    zeta = rib * std::log(zldis / z0m);
    zeta = std::max(-100.0, std::min(zeta, -0.01));
  }

  obu = zldis / zeta;
}

KOKKOS_INLINE_FUNCTION
void StabilityFunc1(Real zeta, Real &value) {

  Real chik2 = std::pow(1.0 - 16.0 * zeta, 0.5);
  Real chik = std::pow(chik2, 0.5);

  Real term1 = 2.0 * std::log((1.0 + chik) * 0.5);
  Real term2 = std::log((1.0 + chik2) * 0.5);
  Real term3 = 2.0 * std::atan(chik);
  value = term1 + term2 + term3 + M_PI * 0.5;
}

KOKKOS_INLINE_FUNCTION
void StabilityFunc2(Real zeta, Real &value) {

  Real chik2 = std::pow(1.0 - 16.0 * zeta, 0.5);
  value = 2.0 * std::log((1.0 + chik2) * 0.5);
}

KOKKOS_INLINE_FUNCTION
void ComputeUstar(Real zldis, Real obu, Real z0m, Real um, Real &ustar) {

  const Real zeta = zldis / obu;

  if (zeta < -ZETAM) {
    const Real term1 = std::log(ZETAM * obu / z0m);
    const Real term4 =
        1.14 * (std::pow(-zeta, 0.333) - std::pow(-ZETAM, 0.333));

    Real term2, term3;
    StabilityFunc1(-ZETAM, term2);
    StabilityFunc1(z0m / obu, term3);

    ustar = VKC * um / (term1 - term2 + term3 + term4);

  } else if (zeta < 0.0) {

    const Real term1 = std::log(zldis / z0m);

    Real term2, term3;
    StabilityFunc1(zeta, term2);
    StabilityFunc1(z0m / obu, term3);

    ustar = VKC * um / (term1 - term2 + term3);

  } else if (zeta <= 1.0) {

    const Real denom = std::log(zldis / z0m) + 5.0 * zeta - 5.0 * z0m / obu;
    ustar = VKC * um / denom;

  } else {

    const Real denom = std::log(obu / z0m) + 5.0 - 5.0 * z0m / obu +
                       (5.0 * std::log(zeta) + zeta - 1.0);
    ustar = VKC * um / denom;
  }
}

KOKKOS_INLINE_FUNCTION
void ComputeU10m(Real zldis, Real obu, Real z0m, Real um, Real ustar,
                 Real &u10) {

  const Real zeta = zldis / obu;

  if (zldis - z0m <= 10.0) {

    u10 = um;

  } else {

    if (zeta < -ZETAM) {

      const Real term1 = std::log(-ZETAM * obu / (10.0 + z0m));
      const Real term4 =
          1.14 * (std::pow(-zeta, 0.333) - std::pow(-ZETAM, 0.333));

      Real term2, term3;
      StabilityFunc1(-ZETAM, term2);
      StabilityFunc1((10.0 + z0m) / obu, term2);

      u10 = um - ustar / VKC * (term1 - term2 + term3 + term4);

    } else if (zeta < 0.0) {

      const Real term1 = std::log(zldis / (10.0 + z0m));

      Real term2, term3;
      StabilityFunc1(zeta, term2);
      StabilityFunc1((10.0 + z0m) / obu, term3);

      u10 = um - ustar / VKC * (term1 - term2 + term3);

    } else if (zeta <= 1.0) {

      const Real term1 = std::log(zldis / (10.0 + z0m));
      const Real term2 = -5.0 * zeta;
      const Real term3 = -5.0 * (10.0 + z0m) / obu;

      u10 = um - ustar / VKC * (term1 - term2 + term3);

    } else {
      const Real term1 = std::log(obu / (10.0 + z0m)) + 5.0;
      const Real term2 = 5.0 * (10.0 + z0m) / obu;
      const Real term3 = (5.0 * std::log(zeta) + zeta - 1.0);

      u10 = um - ustar / VKC * (term1 - term2 + term3);
    }
  }
}

KOKKOS_INLINE_FUNCTION
void ComputeRelationForOtherScalarProfiles(Real zldis, Real obu, Real z0,
                                           Real &value) {

  const Real zeta = zldis / obu;

  if (zeta < -ZETAT) {

    const Real term1 = std::log(-zeta * obu / z0);
    const Real term4 = 0.8 * (std::pow(-ZETAT, 0.333) - std::pow(zeta, 0.0333));

    Real term2, term3;
    StabilityFunc2(-ZETAT, term2);
    StabilityFunc2(z0 / obu, term3);

    value = VKC / (term1 - term2 + term3 + term4);

  } else if (zeta < 0.0) {

    const Real term1 = std::log(zldis / z0);

    Real term2, term3;
    StabilityFunc2(zeta, term2);
    StabilityFunc2(z0 / obu, term3);

    value = VKC / (term1 - term2 + term3);

  } else if (zeta < 1.0) {

    const Real term1 = std::log(zldis / z0) + 5.0 * zeta;
    const Real term2 = 5.0 * z0 / obu;

    value = VKC / (term1 - term2);

  } else {

    const Real term1 = std::log(obu / z0) + 5.0;
    const Real term2 = (5.0 * z0 / obu);
    const Real term3 = 5.0 * std::log(zeta) + zeta - 1.0;

    value = VKC / (term1 - term2 + term3);
  }
}

KOKKOS_INLINE_FUNCTION
void FrictionVelocity(int iter, Real forc_hgt_u, Real displa, Real z0, Real obu,
                      Real ur, Real um, Real &ustar, Real &temp1, Real &temp12m,
                      Real &temp2, Real &temp22m, Real &fm) {

  const Real zldis = forc_hgt_u - displa;
  const Real zeta = zldis / obu;

  const Real z0m = z0;
  const Real z0h = z0;
  const Real z0q = z0;

  ComputeUstar(zldis, obu, z0m, um, ustar);

  Real vds;
  if (zeta < 0) {
    vds = 2.0e-3 * ustar * (1.0 + std::pow(300.0 / (-obu), 0.666));
  } else {
    vds = 2.0e-3 * ustar;
  }

  Real u10;
  ComputeU10m(zldis, obu, z0m, um, ustar, u10);

  // Calcualte temp1 for the temperature profile
  ComputeRelationForOtherScalarProfiles(zldis, obu, z0h, temp1);

  // Since z0q == z0h, temp2 for the humidity profile is same as
  // that for temperature profile
  temp2 = temp1;

  // Calcualte temp1 at 2m for the temperature profile
  ComputeRelationForOtherScalarProfiles(2.0 + z0h, obu, z0h, temp12m);

  // similarly, set temp2 at 2m for humidity profile to be same that
  // for the temperature profile
  temp22m = temp12m;

  Real fmnew;
  if (std::min(zeta, 1.0) < 0.0) {
    const Real x = std::pow((1.0 - 16.0 * std::min(zeta, 1.0)), 0.25);
    const Real tmp2 = std::log((1.0 + x * x) / 2.0);
    const Real tmp3 = std::log((!.0 + x) / 2.0);
    fmnew = 2.0 * tmp3 + tmp2 - std::atan(x) + M_PI / 2.0;
  } else {
    fmnew = -5.0 * std::min(zeta, 1.0);
  }

  if (iter == 0)
    fm = fmnew;
  else
    fm = 0.5 * (fm + fmnew);
}

KOKKOS_INLINE_FUNCTION
void ComputeCanyonUWind(Real ht_roof, Real z_d_town, Real z_0_town,
                        Real forc_hgt_u, Real wind_hgt_canyon, Real hwr,
                        Real ur, Real &canyon_u_wind) {

  // wind at cannyon top
  const Real term1 = std::log((ht_roof - z_d_town) / z_0_town);
  const Real term2 = std::log((forc_hgt_u - z_d_town) / z_0_town);
  const Real canyontop_wind = ur * term1 / term2;

  const Real factor = std::exp(-0.5 * hwr * (1.0 - wind_hgt_canyon / ht_roof));

  if (hwr < 0.5) {

    canyon_u_wind = canyontop_wind * factor; // eqn 3.61

  } else if (hwr < 1.0) {

    const Real factor2 = 1.0 + 2.0 * (2.0 / M_PI - 1.0) * (hwr - 0.5);
    canyon_u_wind = canyontop_wind * factor2 * factor; // eqn 3.62

  } else {

    canyon_u_wind = canyontop_wind * 2.0 / M_PI * factor; // eqn 3.60
  }
}

KOKKOS_INLINE_FUNCTION
void UrbanSurfaceFluxes::computeQsatForSurfaces(int c) {

  const Real forc_p = ForcPbot(c);

  Roof.ComputeQsat(c, forc_p);
  SunlitWall.ComputeQsat(c, forc_p);
  ShadeWall.ComputeQsat(c, forc_p);
  ImperviousRoad.ComputeQsat(c, forc_p);
  PerviousRoad.ComputeQsat(c, forc_p);
}

KOKKOS_INLINE_FUNCTION
void UrbanSurfaceFluxes::computeNewTafAndQaf(int c, Real canyon_wind, Real thm,
                                             Real rahu, Real rawu, Real &taf,
                                             Real &qaf) {

  const Real hwr = data_bundle.geometry.CanyonHwr(c);

  const Real q_roof = Roof.qs(c);
  const Real q_road_imperv = ImperviousRoad.qs(c);
  const Real q_road_perv = PerviousRoad.qs(c);
  const Real q_sunwall = 0.0;
  const Real q_shadewall = 0.0;

  const Real T_roof = Roof.Temp(c);
  const Real T_road_imperv = ImperviousRoad.Temp(c);
  const Real T_road_perv = PerviousRoad.Temp(c);
  const Real T_sunwall = SunlitWall.Temp(c);
  const Real T_shadewall = ShadeWall.Temp(c);

  const Real forc_q = ForcQ(c);
  const Real forc_rho = ForcRho(c);
  Real canyon_resistance = CPAIR * forc_rho / (11.8 + 4.2 * canyon_wind);

  const Real wt_road_perv = fPervRoad(c);
  const Real wt_roof = fRoof(c);

  const Real fwet_roof = 0.0;
  const Real wtus_roof = wt_roof / canyon_resistance;
  const Real wtuq_roof = fwet_roof * wt_roof / canyon_resistance;

  const Real wtus_road_perv =
      wt_road_perv * (1.0 - wt_roof) / canyon_resistance;
  const Real wtuq_road_perv =
      wt_road_perv * (1.0 - wt_roof) / canyon_resistance;

  const Real fwet_road_imperv = (qaf > q_road_imperv) ? 1.0 : 0.0;
  const Real wtus_road_imperv =
      (1.0 - wt_road_perv) * (1.0 - wt_roof) / canyon_resistance;
  const Real wtuq_road_imperv = fwet_road_imperv * (1.0 - wt_road_perv) *
                                (1.0 - wt_roof) / canyon_resistance;

  const Real wtus_sunwall = hwr * (1.0 - wt_roof) / canyon_resistance;
  const Real wtuq_sunwall = 0.0;

  const Real wtus_shadewall = hwr * (1.0 - wt_roof) / canyon_resistance;
  const Real wtuq_shadewall = 0.0;

  const Real taf_numer =
      thm / rahu + T_roof * wtus_roof + T_road_perv * wtus_road_perv +
      T_road_imperv * wtus_road_imperv + T_sunwall * wtus_sunwall +
      T_shadewall * wtus_shadewall;

  const Real taf_denom = 1.0 / rahu + wtus_roof + wtus_road_perv +
                         wtus_road_imperv + wtus_sunwall + wtus_shadewall;

  const Real qaf_numer =
      forc_q / rawu + q_roof * wtuq_roof + q_road_perv * wtuq_road_perv +
      q_road_imperv * wtuq_road_imperv + q_sunwall * wtuq_sunwall +
      q_shadewall * wtuq_shadewall;

  const Real qaf_denom = 1.0 / rawu + wtuq_roof + wtuq_road_perv +
                         wtuq_road_imperv + wtuq_sunwall + wtuq_shadewall;

  taf = taf_numer / taf_denom;
  qaf = qaf_numer / qaf_denom;
}

void UrbanSurfaceFluxes::computeSurfaceFluxes() {
  std::cout << "In computeSurfaceFluxes \n";

  int N_LUN = data_bundle.N_LUN;
  Kokkos::parallel_for(
      "computeNetSurfaceFluxes", N_LUN, KOKKOS_LAMBDA(const int c) {
        const Real forc_hgt_t = 144.44377627618979; // observational height (m)
        const Real forc_hgt_u = forc_hgt_t;         // observational height (m)
        const Real z_d_town = 113.96331622200367;   // displacement height (m)
        const Real z_0_town =
            0.48046005418613641;           // momentum roughness length (m)
        const Real ht_roof = 120.0;        // height of roof (m)
        const Real wind_hgt_canyon = 60.0; // height above road at which in
                                           // canyon needs to be computed (m)
        const Real lapse_rate = 0.0098;    // dry adiabatic lapse rate (K/m)

        Real taf = Taf(c);
        Real qaf = Qaf(c);

        const Real forc_u = ForcU(c);
        const Real forc_v = ForcV(c);
        const Real u2_plus_v2 = std::pow(forc_u, 2.0) + std::pow(forc_v, 2.0);
        const Real velocity = std::pow(u2_plus_v2, 0.5);
        const Real ur = std::max(1.0, velocity);

        const Real hwr = data_bundle.geometry.CanyonHwr(c);

        Real um, obu;
        Real thm, thv, zldis;
        const Real forc_q = ForcQ(c);
        const Real forc_th = ForcTh(c);

        {
          const Real forc_t = ForcTemp(c);

          thm = forc_t + lapse_rate * forc_hgt_t;
          thv = forc_th * (1.0 + 0.61 * forc_q);
          const Real dth = thm - taf;
          const Real dqh = forc_q - qaf;
          const Real dthv = dth * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * dqh;
          zldis = forc_hgt_u - z_d_town;

          MoninObukIni(ur, thv, dthv, zldis, z_0_town, um, obu);
        }

        computeQsatForSurfaces(c);
        Real canyon_u_wind;
        ComputeCanyonUWind(ht_roof, z_d_town, z_0_town, forc_hgt_u,
                           wind_hgt_canyon, hwr, ur, canyon_u_wind);

        Real fm = 0.0;
        for (int iter = 0; iter < 3; ++iter) {
          Real ustar;
          Real temp1, temp12m;
          Real temp2, temp22m;
          FrictionVelocity(iter, forc_hgt_u, z_d_town, z_0_town, obu, ur, um,
                           ustar, temp1, temp12m, temp2, temp22m, fm);

          Real ramu = 1.0 / (ustar * ustar / um);
          Real rahu = 1.0 / (temp1 * ustar);
          Real rawu = 1.0 / (temp2 * ustar);

          Real canyon_wind_pow2 =
              std::pow(canyon_u_wind, 2.0) + std::pow(ustar, 2.0);
          Real canyon_wind = std::pow(canyon_wind_pow2, 0.5);

          computeNewTafAndQaf(c, canyon_wind, thm, rahu, rawu, taf, qaf);
          const Real dth = thm - taf;
          const Real dqh = forc_q - qaf;
          const Real tstar = temp1 * dth;
          const Real qstar = temp2 * dqh;
          const Real thvstar =
              tstar * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * qstar;
          Real zeta =
              zldis * VKC * GRAVITY * thvstar / (std::pow(ustar, 2.0) * thv);

          if (zeta > 0.0) {
            zeta = std::min(2.0, std::max(zeta, 0.01));
            um = std::max(ur, 0.1);
          } else {
            const Real beta = 1.0;
            const Real zii = 1000.0;
            const Real wc =
                beta * std::pow(-GRAVITY * ustar * thvstar * zii / thv, 0.333);
            um = std::pow(ur * ur + wc * wc, 0.5);
          }
          obu = zldis / zeta;
        }

        Taf(c) = taf;
        Qaf(c) = qaf;
      });
}

} // namespace URBANXX