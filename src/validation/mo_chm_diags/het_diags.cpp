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
using namespace mo_chm_diags;

//constexpr const int gas_pcnst = gas_chemistry::gas_pcnst;

void het_diags(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {
  
    const auto het_rates_in = input.get_array("het_rates");
    const auto mmr_in = input.get_array("mmr");
    const auto pdel_in = input.get_array("pdel");
    const auto wght = input.get_array("wght")[0];

    Real het_rates[pver][gas_pcnst];
    Real mmr[pver][gas_pcnst];
    Real pdel[pver];
    Real wrk_wd[gas_pcnst];
    Real sox_wk;
    Real sox_species[3] = {4, -1, 3};
    Real adv_mass[gas_pcnst] = {47.998200,      34.013600,      98.078400,      64.064800,      62.132400, 
                               12.011000,     115.107340,      12.011000,      12.011000,      12.011000, 
                              135.064039,      58.442468,  250092.672000,       1.007400,     115.107340, 
                               12.011000,      58.442468,  250092.672000,       1.007400,     135.064039, 
                               58.442468,     115.107340,      12.011000,      12.011000,      12.011000, 
                           250092.672000,       1.007400,      12.011000,      12.011000,  250092.672000, 
                                1.007400};

    int counter = 0;
    for(int mm = 0; mm < gas_pcnst; mm++) {
      for(int k = 0; k < pver; k++) {
        het_rates[k][mm] = het_rates_in[counter];
        mmr[k][mm] = mmr_in[counter];
        counter++;
      }
      wrk_wd[mm] = 0;
    }
    sox_wk = 0;

    counter = 0;
    for(int k = 0; k < pver; k++) {
      pdel[k] = pdel_in[counter];
      counter++;
    }

    mo_chm_diags::het_diags(het_rates, mmr, pdel, wght, wrk_wd, sox_wk, adv_mass, sox_species); 

    std::vector<Real> wrk_wd_mm(gas_pcnst);
    counter = 0;
    for(int mm = 0; mm < gas_pcnst; mm++) {
      wrk_wd_mm[counter] = wrk_wd[mm];
    }

    output.set("wrk_wd_mm", wrk_wd_mm);
    output.set("sox_wk", sox_wk);

  });
}
