#include <mam4xx/mam4.hpp>
#include <catch2/catch.hpp>
#include <ekat/logging/ekat_logger.hpp>
#include <ekat/mpi/ekat_comm.hpp>

TEST_CASE("test_constructor", "mam4_aging_process") {
  mam4::AeroConfig mam4_config;
  mam4::AgingProcess process(mam4_config);
  REQUIRE(process.name() == "MAM4 aging");
  REQUIRE(process.aero_config() == mam4_config);
}

TEST_CASE("test_compute_tendencies", "mam4_aging_process") {


  ekat::Comm comm;
  ekat::logger::Logger<> logger("aging unit tests",
                            ekat::logger::LogLevel::debug, comm);
  std::ostringstream ss;


  mam4::AeroConfig mam4_config;
  mam4::AgingProcess process(mam4_config);


    ss << "\n aging compute tendencies";
    logger.debug(ss.str());
    ss.str("");  


}

TEST_CASE("test_cond_coag_mass_to_accum", "mam4_aging_process"){

ekat::Comm comm;
ekat::logger::Logger<> logger("aging unit tests",
                          ekat::logger::LogLevel::debug, comm);
std::ostringstream ss;




//const auto naero =  mam4::AeroConfig::num_aerosol_ids();
//const auto nmodes = mam4::AeroConfig::num_modes();



}

TEST_CASE("transfer_aged_pcarbon_to_accum", "mam4_aging_process"){




}

TEST_CASE("mam4_pcarbond_aging_1subarea", "mam4_aging_process"){


}