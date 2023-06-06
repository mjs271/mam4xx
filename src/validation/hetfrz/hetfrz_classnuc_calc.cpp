// mam4xx: Copyright (c) 2022,
// Battelle Memorial Institute and
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
// SPDX-License-Identifier: BSD-3-Clause

#include <iostream>
#include <mam4xx/hetfrz.hpp>
#include <skywalker.hpp>
#include <validation.hpp>

using namespace skywalker;
using namespace mam4;

void hetfrz_classnuc_calc(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {
    if (!input.has("deltat")) {
      std::cerr << "Required name: "
                << "deltat" << std::endl;
      exit(1);
    }

    if (!input.has("temperature")) {
      std::cerr << "Required name: "
                << "temperature" << std::endl;
      exit(1);
    }

    if (!input.has("pressure")) {
      std::cerr << "Required name: "
                << "pressure" << std::endl;
      exit(1);
    }

    if (!input.has("supersatice")) {
      std::cerr << "Required name: "
                << "supersatice" << std::endl;
      exit(1);
    }

    if (!input.has_array("fn")) {
      std::cerr << "Required name: "
                << "fn" << std::endl;
      exit(1);
    }

    if (!input.has("r3lx")) {
      std::cerr << "Required name: "
                << "r3lx" << std::endl;
      exit(1);
    }

    if (!input.has_array("hetraer")) {
      std::cerr << "Required name: "
                << "hetraer" << std::endl;
      exit(1);
    }

    if (!input.has_array("awcam")) {
      std::cerr << "Required name: "
                << "awcam" << std::endl;
      exit(1);
    }

    if (!input.has_array("awfacm")) {
      std::cerr << "Required name: "
                << "awfacm" << std::endl;
      exit(1);
    }

    if (!input.has_array("dstcoat")) {
      std::cerr << "Required name: "
                << "dstcoat" << std::endl;
      exit(1);
    }

    if (!input.has_array("total_aer_num")) {
      std::cerr << "Required name: "
                << "total_aer_num" << std::endl;
      exit(1);
    }

    if (!input.has_array("coated_aer_num")) {
      std::cerr << "Required name: "
                << "coated_aer_num" << std::endl;
      exit(1);
    }

    if (!input.has_array("uncoated_aer_num")) {
      std::cerr << "Required name: "
                << "uncoated_aer_num" << std::endl;
      exit(1);
    }

    if (!input.has_array("total_interstitial_aer_num")) {
      std::cerr << "Required name: "
                << "total_interstitial_aer_num" << std::endl;
      exit(1);
    }

    if (!input.has_array("total_cloudborne_aer_num")) {
      std::cerr << "Required name: "
                << "total_cloudborne_aer_num" << std::endl;
      exit(1);
    }

    // Now read in the values and arrays from the input ensemble
    auto deltat = input.get("deltat");
    auto temperature = input.get("temperature");
    auto pressure = input.get("pressure");
    auto supersatice = input.get("supersatice");
    auto fn = input.get_array("fn");
    auto r3lx = input.get("r3lx");
    auto icnlx = input.get("icnlx");
    auto hetraer = input.get_array("hetraer");
    auto awcam = input.get_array("awcam");
    auto awfacm = input.get_array("awfacm");
    auto dstcoat = input.get_array("dstcoat");
    auto total_aer_num = input.get_array("total_aer_num");
    auto coated_aer_num = input.get_array("coated_aer_num");
    auto uncoated_aer_num = input.get_array("uncoated_aer_num");
    auto total_interstitial_aer_num =
        input.get_array("total_interstitial_aer_num");
    auto total_cloudborne_aer_num = input.get_array("total_cloudborne_aer_num");

    skywalker::Real frzbcimm = 0.0;
    skywalker::Real frzdiumm = 0.0;
    skywalker::Real frzbccnt = 0.0;
    skywalker::Real frzducnt = 0.0;
    skywalker::Real frzbcdep = 0.0;
    skywalker::Real frzdudep = 0.0;

    // call hetfrz::hetfrz_classnuc_calc
    hetfrz::hetfrz_classnuc_calc(
        deltat, temperature, pressure, supersatice, fn.data(), r3lx, icnlx,
        hetraer.data(), awcam.data(), awfacm.data(), dstcoat.data(),
        total_aer_num.data(), coated_aer_num.data(), uncoated_aer_num.data(),
        total_interstitial_aer_num.data(), total_cloudborne_aer_num.data(),
        frzbcimm, frzdiumm, frzbccnt, frzducnt, frzbcdep, frzdudep);

    output.set("frzbcimm", frzbcimm);
    output.set("frzduimm", frzdiumm);
    output.set("frzbccnt", frzbccnt);
    output.set("frzducnt", frzducnt);
    output.set("frzbcdep", frzbcdep);
    output.set("frzdudep", frzdudep);
  });
}