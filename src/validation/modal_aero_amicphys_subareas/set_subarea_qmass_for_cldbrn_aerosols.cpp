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
void compute_qsub_from_gcm_and_qsub_of_other_subarea(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {
    // Ensemble parameters
    // Declare array of strings for input names
    std::string input_arrays[] = {};

    // Iterate over input_arrays and error if not in input
    for (std::string name : input_arrays) {
      if (!input.has_array(name.c_str())) {
        std::cerr << "Required name for array: " << name << std::endl;
        exit(1);
      }
    }

    // using mam4::gas_chemistry::extcnt;
    // using mam4::mo_photo::PhotoTableData;
    // using mam4::mo_setext::Forcing;
    // using mam4::mo_setinv::num_tracer_cnst;

    // using View2D = DeviceType::view_2d<Real>;
    // using ConstView2D = DeviceType::view_2d<const Real>;
    // using View1D = DeviceType::view_1d<Real>;

    // const Real dt = input.get_array("delt")[0];
    // const Real rlats;
    // const View1D cnst_offline_icol[num_tracer_cnst];
    // const Forcing *forcings_in;
      // struct Forcing {
      //   // This index is in Fortran format. i.e. starts in 1
      //   int frc_ndx;
      //   bool file_alt_data;
      //   View1D fields_data[MAX_NUM_SECTIONS];
      //   int nsectors;
      // };


    const auto linoz_o3_clim_icol_ = input.get_array("linoz_o3_clim")[0];
    const Real pblh = input.get_array("pblh")[0];

    // const Atmosphere atm = validation::create_atmosphere(nlev, pblh);
    // mam4::Prognostics progs = validation::create_prognostics(nlev);
    // const mam4::mo_setsox::Config config_setsox;
    // mam4::microphysics::AmicPhysConfig config_amicphys;

    const View1D linoz_o3_clim_icol;

    auto team_policy = ThreadTeamPolicy(1u, Kokkos::AUTO);
    Kokkos::parallel_for(
        team_policy, KOKKOS_LAMBDA(const ThreadTeam &team) {
          // mam4::perform_atmospheric_chemistry_and_microphysics
        });


    std::vector<Real> cflx_out;

    output.set("cflx", cflx_out);
  });
}
