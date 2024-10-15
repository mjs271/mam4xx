// mam4xx: Copyright (c) 2022,
// Battelle Memorial Institute and
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
// SPDX-License-Identifier: BSD-3-Clause

#include <mam4xx/mam4.hpp>

#include <mam4xx/aero_config.hpp>
#include <skywalker.hpp>
#include <validation.hpp>

using namespace skywalker;
using namespace mam4;
using namespace haero;
void aero_model_emissions_test(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {
    // Ensemble parameters
    // Declare array of strings for input names
    std::string input_arrays[] = {
        "dust_flux_in", "dust_dmt_vwr", "srf_temp",
        "ubot",         "vbot",         "zbot",
        "ocnfrc",       "cflx_in",      "soil_erodibility",
        "mpoly_in",     "mprot_in",     "mlip_in"};

    // Iterate over input_arrays and error if not in input
    for (std::string name : input_arrays) {
      if (!input.has_array(name.c_str())) {
        std::cerr << "Required name for array: " << name << std::endl;
        exit(1);
      }
    }

    const int dust_nflux_in_ = mam4::aero_model_emissions::dust_nflux_in;
    const int dust_nbin_ = mam4::aero_model_emissions::dust_nbin;
    const auto dust_indices_ = mam4::aero_model_emissions::dust_indices;

    const auto dust_flux_in_ = input.get_array("dust_flux_in");
    const auto dust_dmt_vwr_ = input.get_array("dust_dmt_vwr");
    const Real surface_temp_ = input.get_array("srf_temp")[0];
    const Real u_bottom_ = input.get_array("ubot")[0];
    const Real v_bottom_ = input.get_array("vbot")[0];
    const Real z_bottom_ = input.get_array("zbot")[0];
    const Real ocean_frac_ = input.get_array("ocnfrc")[0];
    const auto cflux_ = input.get_array("cflx_in");
    Real soil_erodibility_ = input.get_array("soil_erodibility")[0];
    const Real mpoly_ = input.get_array("mpoly_in")[0];
    const Real mprot_ = input.get_array("mprot_in")[0];
    const Real mlip_ = input.get_array("mlip_in")[0];

    constexpr int pcnst = mam4::pcnst;
    Real cflux[pcnst];
    for (int i = 0; i < pcnst; ++i) {
      cflux[i] = cflux_[i];
    }
    for (int i = 0; i < dust_nflux_in_; ++i) {
      cflux[dust_indices_[i]] = 0.0;
    }

    mam4::aero_model_emissions::OnlineEmissionsData online_emiss_data;
    for (int i = 0; i < dust_nflux_in_; ++i) {
      online_emiss_data.dust_flux_in[i] = dust_flux_in_[i];
    }
    online_emiss_data.surface_temp = surface_temp_;
    online_emiss_data.u_bottom = u_bottom_;
    online_emiss_data.v_bottom = v_bottom_;
    online_emiss_data.z_bottom = z_bottom_;
    online_emiss_data.ocean_frac = ocean_frac_;
    online_emiss_data.soil_erodibility = soil_erodibility_;
    mam4::aero_model_emissions::SeasaltEmissionsData seasalt_data;
    seasalt_data.mpoly = mpoly_;
    seasalt_data.mprot = mprot_;
    seasalt_data.mlip = mlip_;
    mam4::aero_model_emissions::DustEmissionsData dust_data;
    for (int i = 0; i < dust_nbin_; ++i) {
      dust_data.dust_dmt_vwr[i] = dust_dmt_vwr_[i];
    }

    mam4::aero_model_emissions::aero_model_emissions(
        online_emiss_data, seasalt_data, dust_data, cflux);

    // =========================================================================
    // NOTE: (N) indicates idx 'N' is changed in multiple function calls
    // =========================================================================
    // dust_indices = {19, 28, (22), (35)}
    // seasalt_indices = {20, 25, 29, (21), (26), (38), (22), (27), (35), (39)}
    // marine_organic_matter_indices = {(21), (26), (38), (22), (27), (39)}
    // =========================================================================
    // REARRANGED IN ORDER:
    // dust_indices = {19, (22), 28, (35)}
    // seasalt_indices = {20, (21), (22), 25, (26), (27), 29, (35), (38), (39)}
    // marine_organic_matter_indices = {(21), (22), (26), (27), (38), (39)}
    // =========================================================================
    // =========================================================================
    // calculations for 'marine_organic_matter_indices'
    // =========================================================================
    // marine_organic_matter_indices =
    //   seasalt_indices[nsalt + nsalt_om + organic_num_idx] AND
    //   seasalt_indices[nsalt + ispec]
    // WHERE organic_num_idx[organic_num_modes] = {0, 1, 3} AND
    // nsalt = 3 AND nsalt_om = 3
    // THUS marine_organic_matter_indices =
    //   seasalt_indices[3 + 3 + {0, 1, 3}] AND seasalt_indices[3 + {0, 1, 2}]
    //   = seasalt_indices[{6, 7, 9}] AND seasalt_indices[{3, 4, 5}]
    //   = seasalt_indices[{3, 4, 5, 6, 7, 9}]
    // FINALLY
    // marine_organic_matter_indices = {21, 26, 38, 22, 27, 39}
    // =========================================================================
    // =========================================================================
    // cflux_out.push_back(cflux[19]);
    // cflux_out.push_back(cflux[28]);
    // cflux_out.push_back(cflux[22]);
    // cflux_out.push_back(cflux[35]);
    std::vector<Real> cflux_out;
    for (int i = 0; i < pcnst; ++i) {
      cflux_out.push_back(cflux[i]);
      std::cout << "cflux[i] = " << cflux[i] << "\n";
    }
    output.set("cflx", cflux_out);
  });
}
