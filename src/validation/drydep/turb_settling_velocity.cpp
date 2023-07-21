// mam4xx: Copyright (c) 2022,
// Battelle Memorial Institute and
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
// SPDX-License-Identifier: BSD-3-Clause

#include <iostream>
#include <mam4xx/drydep.hpp>
#include <skywalker.hpp>
#include <validation.hpp>

using namespace skywalker;
using namespace mam4;

void turb_settling_velocity(Ensemble *ensemble) {

  // Run the ensemble.
  ensemble->process([=](const Input &input, Output &output) {
    EKAT_REQUIRE_MSG(input.has("moment"), "Required name: moment");
    EKAT_REQUIRE_MSG(input.has_array("fraction_land_use"),
                     "Required name: moment: fraction_land_use");
    EKAT_REQUIRE_MSG(input.has("radius_max"), "Required name: radius_max");
    EKAT_REQUIRE_MSG(input.has("tair"), "Required name: tair");
    EKAT_REQUIRE_MSG(input.has("pmid"), "Required name: pmid");
    EKAT_REQUIRE_MSG(input.has("radius_part"), "Required name: radius_part");
    EKAT_REQUIRE_MSG(input.has("density_part"), "Required name: density_part");
    EKAT_REQUIRE_MSG(input.has("sig_part"), "Required name: sig_part");
    EKAT_REQUIRE_MSG(input.has("fricvel"), "Required name: fricvel");
    EKAT_REQUIRE_MSG(input.has("ram1"), "Required name: ram1");
    EKAT_REQUIRE_MSG(input.has("vlc_grv"), "Required name: vlc_trb");

    auto moment = int(input.get("moment"));
    auto fraction_land_use = input.get_array("fraction_land_use");
    auto radius_max = input.get("radius_max");
    auto tair = input.get("tair");
    auto pmid = input.get("pmid");
    auto radius_part = input.get("radius_part");
    auto density_part = input.get("density_part");
    auto sig_part = input.get("sig_part");
    auto fricvel = input.get("fricvel");
    auto ram1 = input.get("ram1");
    auto vlc_grv = input.get("vlc_grv");

    Real vlc_dry = 0.0;
    Real vlc_trb = 0.0;

    drydep::modal_aero_turb_drydep_velocity(
        moment, fraction_land_use.data(), radius_max, tair, pmid, radius_part,
        density_part, sig_part, fricvel, ram1, vlc_grv, vlc_trb, vlc_dry);

    output.set("vlc_grv", vlc_grv);
    output.set("vlc_trb", vlc_trb);
  });
}