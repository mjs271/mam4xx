// mam4xx: Copyright (c) 2022,
// Battelle Memorial Institute and
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MAM4XX_DRYDEP_HPP
#define MAM4XX_DRYDEP_HPP

#include <haero/atmosphere.hpp>
#include <mam4xx/aero_config.hpp>

namespace mam4 {

class DryDep {

public:
  struct Config {

    Config(){};

    Config(const Config &) = default;
    ~Config() = default;
    Config &operator=(const Config &) = default;
  };

private:
  Config config_;

public:
  // name -- unique name of the process implemented by this class
  const char *name() const { return "MAM4 dry deposition"; }

  // init -- initializes the implementation with MAM4's configuration
  void init(const AeroConfig &aero_config,
            const Config &process_config = Config());

  // The value of n_land_type is taken from mozart as defined in mo_drydep.F90
  static constexpr int n_land_type = 11;

  // validate -- validates the given atmospheric state and prognostics against
  // assumptions made by this implementation, returning true if the states are
  // valid, false if not
  KOKKOS_INLINE_FUNCTION
  bool validate(const AeroConfig &config, const ThreadTeam &team,
                const Atmosphere &atm, const Prognostics &progs) const;

  // compute_tendencies -- computes tendencies and updates diagnostics
  // NOTE: that both diags and tends are const below--this means their views
  // NOTE: are fixed, but the data in those views is allowed to vary.
  KOKKOS_INLINE_FUNCTION
  void compute_tendencies(const AeroConfig &config, const ThreadTeam &team,
                          Real t, Real dt, const Atmosphere &atm,
                          const Prognostics &progs, const Diagnostics &diags,
                          const Tendencies &tends) const;
};

namespace drydep {


 const static int pver = 72;

//==============================================================================
// Calculate the radius for a moment of a lognormal size distribution
//==============================================================================
KOKKOS_INLINE_FUNCTION
Real radius_for_moment(const int moment, const Real sig_part,
                       const Real radius_part, const Real radius_max) {
  const Real lnsig = haero::log(sig_part);
  return haero::min(radius_max, radius_part) *
         haero::exp((moment - 1.5) * lnsig * lnsig);
}

//==========================================================================
// Calculate dynamic viscosity of air, unit [kg m-1 s-1]. See RoY94 p. 102
//==========================================================================
KOKKOS_INLINE_FUNCTION
Real air_dynamic_viscosity(const Real temp) {
  // (BAD CONSTANTS)
  return 1.72e-5 * (haero::pow(temp / 273.0, 1.5)) * 393.0 / (temp + 120.0);
}

//==========================================================================
// Calculate kinematic viscosity of air, unit [m2 s-1]
//==========================================================================
KOKKOS_INLINE_FUNCTION
Real air_kinematic_viscosity(const Real temp, const Real pres) {
  const Real vsc_dyn_atm = air_dynamic_viscosity(temp);
  const Real rho = pres / Constants::r_gas_dry_air / temp;
  return vsc_dyn_atm / rho;
}

//======================================================
// Slip correction factor [unitless].
// See, e.g., SeP97 p. 464 and Zhang L. et al. (2001),
// DOI: 10.1016/S1352-2310(00)00326-5, Eq. (3).
// ======================================================
KOKKOS_INLINE_FUNCTION
Real slip_correction_factor(const Real dyn_visc, const Real pres,
                            const Real temp, const Real particle_radius) {

  // [m]
  const Real mean_free_path =
      2.0 * dyn_visc /
      (pres * haero::sqrt(8.0 / (Constants::pi * Constants::r_gas_dry_air *
                                 temp))); // (BAD CONSTANTS)

  const Real slip_correction_factor =
      1.0 +
      mean_free_path *
          (1.257 + 0.4 * haero::exp(-1.1 * particle_radius / mean_free_path)) /
          particle_radius; // (BAD CONSTANTS)

  return slip_correction_factor;
}

//====================================================================
// Calculate the Schmidt number of air [unitless], see SeP97 p.972
//====================================================================
KOKKOS_INLINE_FUNCTION
Real schmidt_number(const Real temp, const Real pres, const Real radius,
                    const Real vsc_dyn_atm, const Real vsc_knm_atm) {

  //  slip correction factor [unitless]
  const Real slp_crc = slip_correction_factor(vsc_dyn_atm, pres, temp, radius);

  // Brownian diffusivity of particle [m2/s], see SeP97 p.474
  const Real dff_aer = Constants::boltzmann * temp * slp_crc /
                       (6.0 * Constants::pi * vsc_dyn_atm * radius);

  return vsc_knm_atm / dff_aer;
}

//=======================================================================================
// Calculate the bulk gravitational settling velocity [m s-1]
//  - using the terminal velocity of sphere falling in a fluid based on Stokes's
//    law and
//  - taking into account the influces of size distribution.
//=======================================================================================
KOKKOS_INLINE_FUNCTION
Real gravit_settling_velocity(const Real particle_radius,
                              const Real particle_density,
                              const Real slip_correction,
                              const Real dynamic_viscosity,
                              const Real particle_sig) {

  // Calculate terminal velocity following, e.g.,
  //  -  Seinfeld and Pandis (1997),  p. 466
  //  - Zhang L. et al. (2001), DOI: 10.1016/S1352-2310(00)00326-5, Eq. 2.

  const Real gravit_settling_velocity =
      (4.0 / 18.0) * particle_radius * particle_radius * particle_density *
      Constants::gravity * slip_correction / dynamic_viscosity;

  // Account for size distribution (i.e., we are calculating the bulk velocity
  // for a particle population instead of a single particle).

  const Real lnsig = haero::log(particle_sig);
  const Real dispersion = haero::exp(2.0 * lnsig * lnsig);

  return gravit_settling_velocity * dispersion;
}

KOKKOS_INLINE_FUNCTION
Real gamma(const int n_land_type) {
  const Real gamma_array[DryDep::n_land_type] = {
      0.56, 0.54, 0.54, 0.56, 0.56, 0.56, 0.50, 0.54, 0.54, 0.54, 0.54};
  return gamma_array[n_land_type];
}

KOKKOS_INLINE_FUNCTION
Real alpha(const int n_land_type) {
  const Real alpha_array[DryDep::n_land_type] = {
      1.50, 1.20, 1.20, 0.80, 1.00, 0.80, 100.0, 50.0, 2.0, 1.2, 50.0};
  return alpha_array[n_land_type];
}

KOKKOS_INLINE_FUNCTION
Real radius_collector(const int n_land_type) {
  const Real radius_collector_array[DryDep::n_land_type] = {
      10.0e-3, 3.5e-3, 3.5e-3,  5.1e-3, 2.0e-3, 5.0e-3,
      -1.0e0,  -1.0e0, 10.0e-3, 3.5e-3, -1.0e+0};
  return radius_collector_array[n_land_type];
}

KOKKOS_INLINE_FUNCTION
int iwet(const int n_land_type) {
  const int iwet_array[DryDep::n_land_type] = {-1, -1, -1, -1, -1, -1,
                                               1,  -1, 1,  -1, -1};
  return iwet_array[n_land_type];
}

KOKKOS_INLINE_FUNCTION
void modal_aero_turb_drydep_velocity(
    const int moment, const Real fraction_landuse[DryDep::n_land_type],
    const Real radius_max, const Real tair, const Real pmid,
    const Real radius_part, const Real density_part, const Real sig_part,
    const Real fricvel, const Real ram1, const Real vlc_grv, Real vlc_trb,
    Real vlc_dry) {

  /// TODO - figure out how/where we need to resolve

  // Calculate size-INdependent thermokinetic properties of the air
  const Real vsc_dyn_atm = air_dynamic_viscosity(tair);
  const Real vsc_knm_atm = air_kinematic_viscosity(tair, pmid);

  // Calculate the mean radius and Schmidt number of the moment
  const Real radius_moment =
      radius_for_moment(moment, sig_part, radius_part, radius_max);
  const Real shm_nbr =
      schmidt_number(tair, pmid, radius_moment, vsc_dyn_atm, vsc_knm_atm);

  // Initialize deposition velocities averages over different land surface types
  Real vlc_trb_wgtsum = 0.0;
  Real vlc_dry_wgtsum = 0.0;

  // Loop over different land surface types. Calculate deposition velocities of
  // those different surface types. The overall deposition velocity of a grid
  // cell is the area-weighted average of those land-type-specific velocities.
  for (int lt = 0; lt < DryDep::n_land_type; ++lt) {

    const Real lnd_frc = fraction_landuse[lt];

    if (lnd_frc != 0.0) {
      //----------------------------------------------------------------------
      // Collection efficiency of deposition mechanism 1 - Brownian diffusion
      //----------------------------------------------------------------------
      const Real brownian = haero::pow(shm_nbr, (-gamma(lt)));

      //----------------------------------------------------------------------
      // Collection efficiency of deposition mechanism 2 - interception
      //----------------------------------------------------------------------
      Real interception = 0.0;
      const Real rc = radius_collector(lt);
      if (rc > 0.0) {
        // vegetated surface
        interception = 2.0 * haero::square(radius_moment / rc);
      }

      //----------------------------------------------------------------------
      // Collection efficiency of deposition mechanism 3 - impaction
      //----------------------------------------------------------------------
      Real stk_nbr = 0.0;
      if (rc > 0.0) {
        // vegetated surface
        stk_nbr = vlc_grv * fricvel / (Constants::gravity * rc);
      } else {
        // non-vegetated surface
        stk_nbr = vlc_grv * fricvel * fricvel /
                  (Constants::gravity * vsc_knm_atm); //  SeP97 p.965
      }

      static constexpr Real beta =
          2.0; // (BAD CONSTANT) empirical parameter $\beta$ in Eq. (7c) of
               // Zhang L. et al. (2001)
      const Real impaction =
          haero::pow(stk_nbr / (alpha(lt) + stk_nbr),
                     beta); // Eq. (7c) of Zhang L. et al.  (2001)

      //-----------------------------------------------------
      // Stick fraction, Eq. (10) of Zhang L. et al.  (2001)
      //-----------------------------------------------------
      Real stickfrac = 1.0;
      static constexpr Real stickfrac_lowerbnd =
          1.0e-10; // (BAD CONSTANT) lower bound of stick fraction
      if (iwet(lt) < 0) {
        stickfrac =
            haero::max(stickfrac_lowerbnd, haero::exp(-haero::sqrt(stk_nbr)));
      }

      //----------------------------------------------------------------------------------
      // Using the numbers calculated above, compute the quasi-laminar layer
      // resistance following Zhang L. et al. (2001), Eq. (5)
      //----------------------------------------------------------------------------------
      static constexpr Real eps0 =
          3.0; // (BAD CONSTANT) empirical parameter $\varepsilon_0$ in Eq. (5)
               // of Zhang L. et al. (2001)
      const Real rss_lmn = 1.0 / (eps0 * fricvel * stickfrac *
                                  (brownian + interception + impaction));

      //--------------------------------------------------------------------
      // Total resistence and deposition velocity of turbulent deposition,
      // see Eq. (21) of Zender et al. (2003)
      //--------------------------------------------------------------------
      const Real rss_trb = ram1 + rss_lmn + ram1 * rss_lmn * vlc_grv;
      const Real vlc_trb_ontype = 1.0 / rss_trb;

      //--------------------------------------------------------------------------------
      // Contributions to the single-value bulk deposition velocities of the
      // grid cell
      //--------------------------------------------------------------------------------
      vlc_trb_wgtsum = vlc_trb_wgtsum + lnd_frc * (vlc_trb_ontype);
      vlc_dry_wgtsum = vlc_dry_wgtsum + lnd_frc * (vlc_trb_ontype + vlc_grv);
    }
  } // lt=0,n_land_type-1

  vlc_trb = vlc_trb_wgtsum;
  vlc_dry = vlc_dry_wgtsum;
}

//==========================================================================
// Calculate particle velocity of gravitational settling
//==========================================================================
KOKKOS_INLINE_FUNCTION
void modal_aero_gravit_settling_velocity(const int moment,
                                         const Real radius_max, const Real tair,
                                         const Real pmid,
                                         const Real radius_part,
                                         const Real density_part,
                                         const Real sig_part, Real vlc_grv) {

  const Real vsc_dyn_atm = air_dynamic_viscosity(tair);

  const Real radius_moment =
      radius_for_moment(moment, sig_part, radius_part, radius_max);

  const Real slp_crc =
      slip_correction_factor(vsc_dyn_atm, pmid, tair, radius_part);

  vlc_grv = gravit_settling_velocity(radius_moment, density_part, slp_crc,
                                     vsc_dyn_atm, sig_part);
}

//---------------------------------------------------------------------------------
// !DESCRIPTION:
//
// Calc aerodynamic resistance over oceans and sea ice from Seinfeld and Pandis,
// p.963.
//
// Author: Natalie Mahowald
// Code refactor: Hui Wan, 2023
//---------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
void calcram(const Real landfrac, const Real icefrac, const Real ocnfrac,
             const Real obklen, const Real ustar, const Real tair,
             const Real pmid, const Real pdel, const Real ram1_in,
             const Real fv_in, Real ram1_out, Real fv_out) {

  static constexpr Real lndfrc_threshold =
      0.000000001; // (BAD CONSTANT) fraction, unitless
  static constexpr Real zzocn =
      0.0001; // (BAD CONSTANT) Ocean aerodynamic roughness length
  static constexpr Real zzsice =
      0.0400; // (BAD CONSTANT) Sea ice aerodynamic roughness length
  static constexpr Real xkar = 0.4; // (BAD CONSTANT) Von Karman constant
  //---------------------------------------------------------------------------
  // Friction velocity:
  //  - If the grid cell has a land fraction larger than a threshold (~zero),
  //    then use cam_in%fv.
  //  - Otherwise, use the ustar calculated in the atmosphere.
  //---------------------------------------------------------------------------
  if (landfrac > lndfrc_threshold) {
    fv_out = fv_in;
  } else {
    fv_out = ustar;
  }

  // fvitt -- fv == 0 causes a floating point exception in
  // dry dep of sea salts and dust

  if (fv_out == 0.0) {
    fv_out = 1.e-12;
  }

  //-------------------------------------------------------------------
  // Aerodynamic resistence
  //-------------------------------------------------------------------

  if (landfrac > lndfrc_threshold) {
    // If the grid cell has a land fraction larger than a threshold (~zero),
    // simply use cam_in%ram1
    ram1_out = ram1_in;
  } else {
    // If the grid cell has a land fraction smaller than the threshold,
    // calculate aerodynamic resistence

    // calculate psi, psi0, temp
    const Real zz =
        pdel * Constants::r_gas_dry_air * tair / pmid / Constants::gravity /
        2.0; // use half the layer height like Ganzefeld and Lelieveld, 1995

    Real psi = 0.0;
    Real psi0 = 0.0;
    if (obklen != 0.0) {
      psi = haero::min(haero::max(zz / obklen, -1.0), 1.0);
      psi0 = haero::min(haero::max(zzocn / obklen, -1.0), 1.0);
    }
    Real temp = zz / zzocn;

    // special treatment for ice-dominant cells
    if (icefrac > 0.5) {
      psi = 0.0;
      if (obklen > 0.0) {
        psi0 = haero::min(haero::max(zzsice / obklen, -1.0), 1.0);
      }
      temp = zz / zzsice;
    }

    // calculate aerodynamic resistence
    if (psi > 0.0) {
      ram1_out = 1.0 / xkar / ustar * (haero::log(temp) + 4.7 * (psi - psi0));
    } else {
      const Real nu = haero::pow((1.00 - 15.000 * psi), 0.25);
      const Real nu0 = haero::pow(1.000 - 15.000 * psi0, 0.25);

      if (ustar != 0.0) {
        ram1_out =
            1.0 / xkar / ustar *
            (haero::log(temp) +
             haero::log(
                 ((haero::square(nu0) + 1.0) * haero::square(nu0 + 1.0)) /
                 ((haero::square(nu) + 1.0) * haero::square(nu + 1.0))) +
             2.0 * (haero::atan(nu) - haero::atan(nu0)));
      } else {
        ram1_out = 0.0;
      }
    }
  }
}

//==========================================================================================
// Calculate deposition velocities caused by turbulent dry deposition and
// gravitational settling of aerosol particles
//
// Reference:
//  L. Zhang, S. Gong, J. Padro, and L. Barrie:
//  A size-seggregated particle dry deposition scheme for an atmospheric aerosol
//  module Atmospheric Environment, 35, 549-560, 2001.
//
// History:
//  - Original version by X. Liu.
//  - Calculations for gravitational and turbulent dry deposition separated into
//    different subroutines by Hui Wan, 2023.
//==========================================================================================
KOKKOS_INLINE_FUNCTION
void model_aero_depvel_part(const Real fraction_landuse[DryDep::n_land_type],
                            const Real tair, const Real pmid, const Real ram1,
                            const Real fricvel, const Real radius_part,
                            const Real density_part, const Real sig_part,
                            const int moment, Real vlc_dry, Real vlc_trb,
                            Real vlc_grv) {

  static constexpr Real radius_max = 50.0e-6;

  //------------------------------------------------------------------------------------
  // Calculate deposition velocity of gravitational settling in all grid layers
  //------------------------------------------------------------------------------------

  modal_aero_gravit_settling_velocity(moment, radius_max, tair, pmid,
                                      radius_part, density_part, sig_part,
                                      vlc_grv);

  // vlc_dry is just the gravitational settling velocity for now.
  vlc_dry = vlc_grv;

  //------------------------------------------------------------------------------------
  // For the lowest model layer:
  //  - Calculate turbulent dry deposition velocity, vlc_trb.
  //  - Add vlc_trb to vlc_grv to give the total deposition velocity, vlc_dry.
  //------------------------------------------------------------------------------------
  modal_aero_turb_drydep_velocity(moment, fraction_landuse, radius_max, tair,
                                  pmid, radius_part, density_part, sig_part,
                                  fricvel, ram1, vlc_grv, vlc_trb, vlc_dry);
}

//==========================================================================================
// Calculate deposition velocities caused by turbulent dry deposition and
// gravitational settling of aerosol particles
//
// Reference:
//  L. Zhang, S. Gong, J. Padro, and L. Barrie:
//  A size-seggregated particle dry deposition scheme for an atmospheric aerosol
//  module Atmospheric Environment, 35, 549-560, 2001.
//
// History:
//  - Original version by X. Liu.
//  - Calculations for gravitational and turbulent dry deposition separated into
//   different subroutines by Hui Wan, 2023.
//==========================================================================================
KOKKOS_INLINE_FUNCTION
void modal_aero_depvel_part(const Real fraction_landuse[DryDep::n_land_type],
                            const Real tair, const Real pmid, const Real ram1,
                            const Real fricvel, const Real radius_part,
                            const Real density_part, const Real sig_part,
                            const int moment, Real vlc_dry, Real vlc_trb,
                            Real vlc_grv) {

  // use a maximum radius of 50 microns when calculating deposition velocity
  static constexpr Real radius_max = 50.0e-6; //(BAD CONSTANT)

  //------------------------------------------------------------------------------------
  // Calculate deposition velocity of gravitational settling in all grid layers
  //------------------------------------------------------------------------------------
  modal_aero_gravit_settling_velocity(moment, radius_max, tair, pmid,
                                      radius_part, density_part, sig_part,
                                      vlc_grv);

  // vlc_dry is just the gravitational settling velocity for now.
  vlc_dry = vlc_grv;

  //------------------------------------------------------------------------------------
  // For the lowest model layer:
  //  - Calculate turbulent dry deposition velocity, vlc_trb.
  //  - Add vlc_trb to vlc_grv to give the total deposition velocity, vlc_dry.
  //------------------------------------------------------------------------------------

  modal_aero_turb_drydep_velocity(moment, fraction_landuse, radius_max, tair,
                                  pmid, radius_part, density_part, sig_part,
                                  fricvel, ram1, vlc_grv, vlc_trb, vlc_dry);
}

} // namespace drydep

KOKKOS_INLINE_FUNCTION
void drydep_diags_for_1_tracer(const Real vlc_dry, const Real vlc_trb,
                               const Real vlc_grv, Real sflx, Real dep_trb,
                               Real dep_grv) {

  // apportion dry deposition into turb and gravitational settling for tapes
  if (vlc_dry != 0.0) {
    // turbulent dry deposition portion of sflx [kg/m2/s] or [1/m2/s]
    dep_trb = sflx * vlc_trb / vlc_dry;
    // gravitational settling   portion of slfx [kg/m2/s] or [1/m2/s]
    dep_grv = sflx * vlc_grv / vlc_dry;
  }
}

} // namespace mam4

#endif