// mam4xx: Copyright (c) 2022,
// Battelle Memorial Institute and
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MAM4XX_MODAL_AERO_UPTAKE_HPP
#define MAM4XX_MODAL_AERO_UPTAKE_HPP

#include <mam4xx/mam4_types.hpp>
#include <haero/floating_point.hpp>
#include <haero/haero.hpp>


namespace mam4 {

//-----------------------------------------------------------------------
//     solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
//     where p0, p1, p2, p3 are real
//-----------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
void makoh_quartic(Real p3, Real p2, Real p1, Real p0, std::complex<double> cx[4]) {
    //temporary variables
    Real third, qq, rr;
    std::complex<double> cb, cb0, cb1, crad, cy, czero;

    //set complex zeros and 1/3 values
    czero = 0;
    third= 1.0 / 3.0;


    qq = (-p2 * p2) / 36.0 + (p3 * p1 - 4.0 * p0) / 12.0;
    rr = -haero::pow((p2 / 6.0), 3.0) + p2 * (p3 * p1 - 4 * p0) / 48.0 + (4.0 * p0 * p2 - p0 * p3 * p3 - p1 * p1) / 16.0;

    crad = rr*rr + qq*qq*qq;
    crad = haero::sqrt(crad);

    cb = rr-crad;
    printf("crad = %f\n", crad.real());
    printf("cb = %f\n", cb.real()); 
    if(cb == czero) {
        //insoluble particle
        cx[0] = haero::pow(-p1, (1.0/3.0)); //it won't let me fractionally exponent a negative number
        cx[1] = cx[0];
        cx[2] = cx[0];
        cx[3] = cx[0];
    } else {
        printf("crad = %f\n", crad.real());
        printf("cb = %f\n", cb.real());
        printf("qq = %f\n", qq);
        printf("third = %f\n", third);
        cb = haero::pow(cb, (1.0/3.0));
        printf("cb = %f\n", cb.real());
        cy = -cb + qq/cb + p2/6;
        printf("cb = %f\n", cb.real());
        printf("cy = %f\n", cy.real());

        cb0 = haero::sqrt(cy*cy - p0);
        cb1 = (p3*cy - p1)/(2.0*cb0);
        printf("cb0 = %f\n", cb0.real());
        printf("cb1 = %f\n", cb1.real());

        cb = p3/2.0 + cb1;
        crad = cb*cb - 4.0*(cy + cb0);
        crad = haero::sqrt(crad);

        printf("cb = %f\n", cb.real());
        printf("crad = %f\n", crad.real());

        cx[0] = (-cb+crad) / 2.0;
        cx[1] = (-cb-crad) / 2.0;

        cb = p3/2 - cb1;
        crad = cb*cb - 4.0*(cy - cb0);
        crad = haero::sqrt(crad);
        cx[2] = (-cb + crad) / 2.0;
        cx[3] = (-cb - crad) / 2.0;   
    }
}

} //namespace mam4

#endif