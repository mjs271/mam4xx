#include <mam4xx/mam4.hpp>

#include <mam4xx/aero_config.hpp>
#include <mam4xx/rename.hpp>
#include <skywalker.hpp>
#include <validation.hpp>

using namespace skywalker;
using namespace mam4;

void mam_rename_1subarea(Ensemble *ensemble) {

  ensemble->process([=](const Input &input, Output &output) {
    int nlev = 1;
    Real pblh = 1000;
    Atmosphere atm(nlev, pblh);
    mam4::Prognostics progs(nlev);
    mam4::Diagnostics diags(nlev);
    mam4::Tendencies tends(nlev);

    mam4::AeroConfig mam4_config;
    mam4::RenameProcess process(mam4_config);

    const int nmodes = AeroConfig::num_modes();
    const int naerosol_species = AeroConfig::num_aerosol_ids();

    const Real zero = 0;
    const Real smallest_dryvol_value = 1.0e-25; // FIXME: BAD_CONSTANT

    auto iscloudy_vector = input.get_array("iscldy");
    bool iscloudy = iscloudy_vector[0];
    // int dest_mode_of_mode = input.get_array("dest_mode_of_mode");
    int dest_mode_of_mode[nmodes] = {0, 1, 0, 0};

    auto set_mode_values = [](const std::vector<Real> &vector_in,
                              Real values[nmodes]) {
      for (int m = 0; m < nmodes; ++m) {
        values[m] = vector_in[m];
      }
    };


    auto set_mode_aerosol_values = [](const std::vector<Real> &vector_in,
                                       Real values[nmodes][naerosol_species]) {
      int count = 0;
      for (int m = 0; m < nmodes; ++m) {
        for (int ispec = 0; ispec < naerosol_species; ++ispec) {
          values[m][ispec] = vector_in[count];
          count++;
        }
      }
    };

    auto qnum_cur_vector = input.get_array("qnum_cur");
    Real qnum_cur[nmodes];
    set_mode_values(qnum_cur_vector, qnum_cur);

    Real qnumcw_cur[nmodes]={zero};
    Real qaercw_cur[nmodes][naerosol_species]={{zero}};
    Real qaercw_del_grow4rnam[nmodes][naerosol_species]={{zero}};

    if (iscloudy){
        auto qnumcw_cur_vector = input.get_array("qnumcw_cur");
        set_mode_values(qnumcw_cur_vector, qnumcw_cur);

        auto qaercw_cur_vector = input.get_array("qaercw_cur");
        set_mode_aerosol_values(qaercw_cur_vector, qaercw_cur );

        auto qaercw_del_grow4rnam_vector = input.get_array("qaercw_del_grow4rnam");
        set_mode_aerosol_values(qaercw_del_grow4rnam_vector, qaercw_del_grow4rnam );
    }


    Real qaer_del_grow4rnam[nmodes][naerosol_species];
    auto qaer_del_grow4rnam_vector = input.get_array("qaer_del_grow4rnam");
    set_mode_aerosol_values(qaer_del_grow4rnam_vector, qaer_del_grow4rnam );

    auto qaer_cur_vector = input.get_array("qaer_cur");
    Real qaer_cur[nmodes][naerosol_species];
    set_mode_aerosol_values(qaer_cur_vector, qaer_cur );

    Rename this_rename;

    Real sz_factor[nmodes];
    Real fmode_dist_tail_fac[nmodes];
    Real v2n_lo_rlx[nmodes];
    Real v2n_hi_rlx[nmodes];
    Real ln_diameter_tail_fac[nmodes];
    int num_pairs = 0;
    Real diameter_cutoff[nmodes];
    Real ln_dia_cutoff[nmodes];
    Real diameter_threshold[nmodes];
    Real mass_2_vol[naerosol_species];

    rename::find_renaming_pairs(dest_mode_of_mode,    // in
                                sz_factor,            // out
                                fmode_dist_tail_fac,  // out
                                v2n_lo_rlx,           // out
                                v2n_hi_rlx,           // out
                                ln_diameter_tail_fac, // out
                                num_pairs,            // out
                                diameter_cutoff,      // out
                                ln_dia_cutoff,
                                diameter_threshold);

    Real dgnum_amode[nmodes];
    for (int m = 0; m < nmodes; ++m) {
      dgnum_amode[m] = modes(m).nom_diameter;
    }

    Real molecular_weight_rename[naerosol_species] = {
        150, 115, 150, 12, 58.5, 135, 250092}; // [kg/kmol]
    for (int iaero = 0; iaero < naerosol_species; ++iaero) {
      mass_2_vol[iaero] = molecular_weight_rename[iaero] / aero_species(iaero).density;
    }

    this_rename.mam_rename_1subarea_(iscloudy,
                                     smallest_dryvol_value,
                                     dest_mode_of_mode,    // in
                                     sz_factor,            // in
                                     fmode_dist_tail_fac,  // in
                                     v2n_lo_rlx,           // in
                                     v2n_hi_rlx,           // in
                                     ln_diameter_tail_fac, // in
                                     num_pairs,            // in
                                     diameter_cutoff,      // in
                                     ln_dia_cutoff,        // in
                                     diameter_threshold,   // in
                                     mass_2_vol,
                                     dgnum_amode, // in
                                     qnum_cur, 
                                     qaer_cur,
                                     qaer_del_grow4rnam,
                                     qnumcw_cur, qaercw_cur,
                                     qaercw_del_grow4rnam);

    auto save_mode_values = [](const Real values[nmodes],
                               std::vector<Real> &values_vector) {
      for (int i = 0; i < nmodes; ++i)
        values_vector[i] = values[i];
    };

    auto save_mode_aerosol_values =
        [](const Real values[nmodes][naerosol_species],
           std::vector<Real> &values_vector) {
          int count = 0;
          for (int m = 0; m < nmodes; ++m) {
            for (int ispec = 0; ispec < naerosol_species; ++ispec) {
              values_vector[count] = values[m][ispec];
              count++;
            }
          }
        };

    std::vector<Real> qnum_cur_out(nmodes+1, 0);
    save_mode_values(qnum_cur, qnum_cur_out);
    output.set("qnum_cur", qnum_cur_out);

    std::vector<Real> qaer_cur_out((nmodes+1)*naerosol_species,0);
    save_mode_aerosol_values(qaer_cur, qaer_cur_out);
    output.set("qaer_cur", qaer_cur_out);
    
    // std::vector<Real> qaer_del_grow4rnam_out((nmodes+1)*naerosol_species,0);
    // save_mode_aerosol_values(qaer_del_grow4rnam, qaer_del_grow4rnam_out);
    // output.set("qaer_del_grow4rnam", qaer_del_grow4rnam_out);

    if (iscloudy) {
      std::vector<Real> qnumcw_cur_out(nmodes + 1, 0);
      save_mode_values(qnumcw_cur, qnumcw_cur_out);
      output.set("qnumcw_cur", qnumcw_cur_out);

      std::vector<Real> qaercw_cur_out((nmodes + 1) * naerosol_species, 0);
      save_mode_aerosol_values(qaercw_cur, qaercw_cur_out);
      output.set("qaercw_cur", qaercw_cur_out);

      // std::vector<Real> qaercw_del_grow4rnam_out(
      //     (nmodes + 1) * naerosol_species, 0);
      // save_mode_aerosol_values(qaercw_del_grow4rnam, qaercw_del_grow4rnam_out);
      // output.set("qaercw_del_grow4rnam", qaercw_del_grow4rnam_out);
    }    

  });
}
