// mam4xx: Copyright (c) 2022,
// Battelle Memorial Institute and
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
// SPDX-License-Identifier: BSD-3-Clause

#include "mam4xx/aero_modes.hpp"
#include "mam4xx/conversions.hpp"
#include "mam4xx/modal_aero_uptake.hpp"

#include <catch2/catch.hpp>

#include <haero/atmosphere.hpp>
#include <haero/constants.hpp>
#include <haero/floating_point.hpp>
#include <haero/haero.hpp>

#include "atmosphere_utils.hpp"
#include "mam4xx/conversions.hpp"

#include <ekat/logging/ekat_logger.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <cmath>
#include <sstream>

using namespace mam4;

TEST_CASE("makoh_quartic", "") {
    ekat::Comm comm;
    ekat::logger::Logger<> logger("modal aero uptake unit tests",
                                ekat::logger::LogLevel::info, comm);
    Real p0 = 0;
    Real p1 = 0;
    Real p2 = 0;
    Real p3 = 0;
    std::complex<double> x[4] = {0,0,0,0};

    makoh_quartic(p3, p2, p1, p0, x);
    for(int i = 0; i < 4; i++) {
        logger.info("x[{}] = {}", i, x[i].real());
        //REQUIRE(FloatingPoint<Real>::equiv(x[i], 0));
    }
        
    p0 = 1;
    p1 = 4;
    p2 = 6;
    p3 = 4;

    makoh_quartic(p3, p2, p1, p0, x);
    for(int i = 0; i < 4; i++) {
        logger.info("x[{}] = {}", i, x[i].real());
        //REQUIRE(FloatingPoint<Real>::equiv(x[i], 2));
    }
}