#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>
#include <UrbanFluxes.hpp>
#include <cmath>
#include <iostream>

#define ALLOCATE_VIEW(viewname, type, ...)                                     \
  viewname = type(#viewname, __VA_ARGS__);

namespace URBANXX {
UrbanSurfaceFluxes::UrbanSurfaceFluxes(UrbanSharedDataBundle &bundle)
    : data_bundle(bundle) {
  ForcTemp = bundle.input.ForcTemp;
  ForcTh = bundle.input.ForcPotTemp;
  ForcRho = bundle.input.ForcRho;
  ForcQ = bundle.input.ForcSpcHumd;
  ForcPbot = bundle.input.ForcPress;
  ForcU = bundle.input.ForcWindU;
  ForcV = bundle.input.ForcWindV;

  int N_LUN = data_bundle.N_LUN;
  ALLOCATE_VIEW(Taf, Array1DR8, N_LUN);
  ALLOCATE_VIEW(Qaf, Array1DR8, N_LUN);

  Array1DR8 qafH, tafH;
  ALLOCATE_VIEW(qafH, HostArray1DR8, N_LUN);
  ALLOCATE_VIEW(tafH, HostArray1DR8, N_LUN);

  Real TAF = 283.0;
  Real QAF = 1e-4;
  for (int i = 0; i < N_LUN; ++i) {
    tafH(i) = TAF;
    qafH(i) = QAF;
  }

  Kokkos::deep_copy(Taf, tafH);
  Kokkos::deep_copy(Qaf, qafH);
}

#define GRAVITY 9.80616
#define VKC 0.4
#define ZETAM 1.574 // transition point of flux-gradient relation (wind profile)
#define ZETAT                                                                  \
  0.465 // transition point of flux-gradient relation (temperature profile)
#define CPAIR 1004.64 // specific heat of dry air [J/kg/K]

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

    canyon_u_wind = canyontop_wind * M_PI / 2.0 * factor; // eqn 3.60
  }
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

        const Real forc_t = ForcTemp(c);
        const Real forc_u = ForcU(c);
        const Real forc_v = ForcV(c);
        const Real forc_th = ForcTh(c);
        const Real forc_q = ForcQ(c);
        const Real forc_rho = ForcRho(c);

        Real taf = Taf(c);
        Real qaf = Qaf(c);

        const Real u2_plus_v2 = std::pow(forc_u, 2.0) + std::pow(forc_v, 2.0);
        const Real velocity = std::pow(u2_plus_v2, 0.5);
        const Real ur = std::max(1.0, velocity);

        const Real thm = forc_t + lapse_rate * forc_hgt_t;
        const Real thv = forc_th * (1.0 + 0.61 * forc_q);
        const Real dth = thm - taf;
        const Real dqh = forc_q - qaf;
        const Real dthv = dth * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * dqh;
        const Real zldis = forc_hgt_u - z_d_town;

        const Real hwr = data_bundle.geometry.CanyonHwr(c);

        Real um, obu;
        MoninObukIni(ur, thv, dthv, zldis, z_0_town, um, obu);

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

          Real canyon_u_wind;
          if (iter == 0) {
            ComputeCanyonUWind(ht_roof, z_d_town, z_0_town, forc_hgt_u,
                               wind_hgt_canyon, hwr, ur, canyon_u_wind);
          }

          Real canyon_wind_pow2 =
              std::pow(canyon_u_wind, 2.0) + std::pow(ustar, 2.0);
          Real canyon_wind = std::pow(canyon_wind_pow2, 0.5);
          Real canyon_resistance =
              CPAIR * forc_rho / (11.8 + 4.2 * canyon_wind);
        }
      });
}

} // namespace URBANXX