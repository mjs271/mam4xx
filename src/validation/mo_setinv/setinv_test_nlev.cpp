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
void setinv_test_nlev(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {
    using View1D = DeviceType::view_1d<Real>;
    using View2D = DeviceType::view_2d<Real>;
    using View2DHost = typename HostType::view_2d<Real>;
    using View1DHost = typename HostType::view_1d<Real>;
    // Ensemble parameters
    // Declare array of strings for input names
    std::string input_arrays[] = {
        "tfld",      "h2ovmr",  "vmr",       "pmid",
        "ncol",      "lchnk",   "pcols",     "pver",
        "gas_pcnst", "nfs",     "boltz_cgs", "num_tracer_cnst",
        "has_n2",    "m_ndx",   "n2_ndx",    "has_o2",
        "o2_ndx",    "has_h2o", "h2o_ndx",   "cnst_offline_yaml"};

    // Iterate over input_arrays and error if not in input
    for (std::string name : input_arrays) {
      if (!input.has_array(name.c_str())) {
        std::cerr << "Required name for array: " << name << std::endl;
        exit(1);
      }
    }

    const Real tfld_in = input.get_array("tfld")[0];
    const Real h2ovmr_in = input.get_array("h2ovmr")[0];
    const Real pmid_in = input.get_array("pmid")[0];
    auto c_off_in = input.get_array("cnst_offline_yaml");

    const int nlev = mam4::nlev;
    const int nfs = mam4::mo_setinv::nfs;
    const int num_tracer_cnst = mam4::mo_setinv::num_tracer_cnst;

    ColumnView tfld = haero::testing::create_column_view(nlev);
    auto tfld_h = Kokkos::create_mirror_view(tfld);
    ColumnView h2ovmr = haero::testing::create_column_view(nlev);
    auto h2ovmr_h = Kokkos::create_mirror_view(h2ovmr);
    ColumnView pmid = haero::testing::create_column_view(nlev);
    auto pmid_h = Kokkos::create_mirror_view(pmid);

    View2D invariants("invariants", nlev, nfs);
    View1DHost invariants_h("invariants_h", nlev);

    View1D c_off[num_tracer_cnst];
    for (int i = 0; i < num_tracer_cnst; ++i) {
      c_off[i] = View1D("c_off", nlev);
    } //

    View2DHost c_off_h("c_off_h", nlev, num_tracer_cnst);

    for (int k = 0; k < nlev; ++k) {
      tfld_h(k) = tfld_in;
      h2ovmr_h(k) = h2ovmr_in;
      pmid_h(k) = pmid_in;
      for (int i = 0; i < num_tracer_cnst; ++i) {
        c_off_h(k, i) = c_off_in[i];
      }
    }

    Kokkos::deep_copy(tfld, tfld_h);
    Kokkos::deep_copy(h2ovmr, h2ovmr_h);
    Kokkos::deep_copy(pmid, pmid_h);
    for (int i = 0; i < num_tracer_cnst; ++i) {
      const auto c_off_h_at_i = Kokkos::subview(c_off_h, Kokkos::ALL, i);
      Kokkos::deep_copy(c_off[i], c_off_h_at_i);
    }

    // Single-column dispatch.
    auto team_policy = ThreadTeamPolicy(1u, Kokkos::AUTO);
    Kokkos::parallel_for(
        team_policy, KOKKOS_LAMBDA(const ThreadTeam &team) {
          mam4::mo_setinv::setinv(team, invariants, tfld, h2ovmr, c_off, pmid);
        });

    std::vector<Real> invariants_out(nfs);
    for (int i = 0; i < nfs; ++i) {
      const auto invariants_at_i = Kokkos::subview(invariants, Kokkos::ALL, i);
      Kokkos::deep_copy(invariants_h, invariants_at_i);
      invariants_out[i] = invariants_h[i];
      for (int k = 0; k < nlev; ++k) {
        // make sure every level got the same answer
        EKAT_ASSERT(invariants_h(0) == invariants_h(k));
      }
    }

    output.set("invariants", invariants_out);
  });
}
