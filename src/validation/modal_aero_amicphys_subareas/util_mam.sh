#! /bin/bash

# amicphys_subareas_driver.cpp

# compute_qsub_from_gcm_and_qsub_of_other_subarea.cpp
# form_gcm_of_gases_and_aerosols_from_subareas.cpp
# get_partition_factors.cpp
# set_subarea_gases_and_aerosols.cpp
# set_subarea_qmass_for_cldbrn_aerosols.cpp
# set_subarea_qmass_for_intrst_aerosols.cpp
# set_subarea_qnumb_for_cldbrn_aerosols.cpp
# set_subarea_qnumb_for_intrst_aerosols.cpp
# set_subarea_rh.cpp
# setup_subareas.cpp

# compute_qsub_from_gcm_and_qsub_of_other_subarea.yaml
# form_gcm_of_gases_and_aerosols_from_subareas.yaml
# get_partition_factors.yaml
# mam_compute_qsub_from_gcm_and_qsub_of_other_subarea.py
# mam_form_gcm_of_gases_and_aerosols_from_subareas.py
# mam_get_partition_factors.py
# mam_set_subarea_gases_and_aerosols.py
# mam_set_subarea_qmass_for_cldbrn_aerosols.py
# mam_set_subarea_qmass_for_intrst_aerosols.py
# mam_set_subarea_qnumb_for_cldbrn_aerosols.py
# mam_set_subarea_qnumb_for_intrst_aerosols.py
# mam_set_subarea_rh.py
# mam_setup_subareas.py
# set_subarea_gases_and_aerosols.yaml
# set_subarea_qmass_for_cldbrn_aerosols.yaml
# set_subarea_qmass_for_intrst_aerosols.yaml
# set_subarea_qnumb_for_cldbrn_aerosols.yaml
# set_subarea_qnumb_for_intrst_aerosols.yaml
# set_subarea_rh.yaml
# setup_subareas.yaml

files=(
  form_gcm_of_gases_and_aerosols_from_subareas.cpp
  get_partition_factors.cpp
  set_subarea_gases_and_aerosols.cpp
  set_subarea_qmass_for_cldbrn_aerosols.cpp
  set_subarea_qmass_for_intrst_aerosols.cpp
  set_subarea_qnumb_for_cldbrn_aerosols.cpp
  set_subarea_qnumb_for_intrst_aerosols.cpp
  set_subarea_rh.cpp
  setup_subareas.cpp
)

for f in ${files[@]}; do
  cp compute_qsub_from_gcm_and_qsub_of_other_subarea.cpp "${f}"
done
