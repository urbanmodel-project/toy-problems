#ifndef URBAN_ALBEDO_HPP
#define URBAN_ALBEDO_HPP

#include <DataTypes.hpp>
#include <Kokkos_Core.hpp>
#include <UrbanData.hpp>

namespace URBANXX {
class UrbanAlbedo {
private:
  UrbanSharedDataBundle &data_bundle;

public:
  UrbanAlbedo(UrbanSharedDataBundle &bundle);

  void set_solar_inputs() const;
  void compute_incident_radiation() const;
  void compute_snow_albedo() const;
  void compute_combined_albedo() const;
  void compute_net_solar() const;
};
} // namespace URBANXX

#endif