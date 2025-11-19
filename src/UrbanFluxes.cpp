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
#define ZETAM 1.574 // transition point of lux-gradient relation (wind profile)
#define ZETAT                                                                  \
  0.465 // transition point of lux-gradient relation (temperature profile)

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
      const Real term2 = -(5.0 * std::log(zeta) + zeta - 1.0);
      const Real term3 = -5.0 * (10.0 + z0m) / obu;

      u10 = um - ustar / VKC * (term1 - term2 + term3);
    }
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
}

void UrbanSurfaceFluxes::computeSurfaceFluxes() {
  std::cout << "In computeSurfaceFluxes \n";

  int N_LUN = data_bundle.N_LUN;
  Kokkos::parallel_for(
      "computeNetSurfaceFluxes", N_LUN, KOKKOS_LAMBDA(const int c) {
        const Real forc_hgt_t = 144.44377627618979; // observational height (m)
        const Real forc_hgt_u = forc_hgt_t;         // observational height (m)
        const Real z_d_town = 113.96331622200367;   // displacement height (m)
        const Real z_0_town = 0.48046005418613641;  // momentum roughness length
        const Real lapse_rate = 0.0098; // dry adiabatic lapse rate (K/m)

        const Real forc_t = ForcTemp(c);
        const Real forc_u = ForcU(c);
        const Real forc_v = ForcV(c);
        const Real forc_th = ForcTh(c);
        const Real forc_q = ForcQ(c);

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

        Real um, obu;
        MoninObukIni(ur, thv, dthv, zldis, z_0_town, um, obu);

        int iter = 1;
        Real ustar;
        Real temp1, temp12m;
        Real temp2, temp22m;
        Real fm;
        FrictionVelocity(iter, forc_hgt_u, z_d_town, z_0_town, obu, ur, um,
                         ustar, temp1, temp12m, temp2, temp22m, fm);
      });
}

} // namespace URBANXX