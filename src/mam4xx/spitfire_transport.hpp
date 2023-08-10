// mam4xx: Copyright (c) 2022,
// Battelle Memorial Institute and
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MAM4XX_SPITFIRE_HPP
#define MAM4XX_SPITFIRE_HPP

#include <haero/math.hpp>

namespace mam4::spitfire {

 //##############################################################################
 // The minmod function 
 //##############################################################################
 KOKKOS_INLINE_FUNCTION
 Real minmod(const Real aa, const Real bb){
    return 0.5*(Kokkos::copysign(1.0, aa) + Kokkos::copysign(1.0, bb))*haero::min(haero::abs(aa),haero::abs(bb));
 }

//##############################################################################
// The medan function
//##############################################################################
KOKKOS_INLINE_FUNCTION
Real median(const Real aa, const Real bb, const Real cc){
    return aa + minmod(bb-aa,cc-aa);
}

KOKKOS_INLINE_FUNCTION
void cfdotmc_pro(xw, ff, fdot)

KOKKOS_INLINE_FUNCTION
void get_flux(const int nk, const ColumnView &xw, const ColumnView &phi, const ColumnView &vel, const Real deltat, const ColumnView &flux){


    // Set fluxes at boundaries to zero
    flux[0] = 0.0;
    flux[nk-1] = 0.0;

    // Compute the vertical integral of phi
    // See Rasch and Lawrence (1989), Eq (3) but note we are using a pressure coordinate here
    // This needs to be parallelized.
    Real psi[nk];
    for (int k = 1; k < nk-2; ++k){
        psi[k] = phi[k-1] * (xw[k] - xw[k-1]) + psi[k-1];
    }


}

}

#endif