#!/bin/bash

target_dir="HPC_package_HBV_2nov2026"

mkdir HPC_packages/$target_dir
mkdir HPC_packages/$target_dir/src
mkdir HPC_packages/$target_dir/resources
## Copy MatLab scripts:
cp main_script.m HPC_packages/$target_dir/src
cp country_level_analyses.m HPC_packages/$target_dir/src
cp assign_PAP_VL_params.m HPC_packages/$target_dir/src
cp make_daly_mat.m HPC_packages/$target_dir/src
cp HBVmodel.m HPC_packages/$target_dir/src
cp get_range_of_FracSPosHighVL.m HPC_packages/$target_dir/src
cp safe_greater_or_equal_to.m HPC_packages/$target_dir/src
cp get_treatment_eligible_ageindices.m HPC_packages/$target_dir/src

## Copy parameter files: pre-existing (MJdV) .mat files, and new (.csv) files.
cp ../resources/*.mat ../src/HPC_packages/$target_dir/resources
cp ../resources/GHO_in_facility_births_summary.csv HPC_packages/$target_dir/resources/
cp ../resources/ANC_coverage_data.csv HPC_packages/$target_dir/resources/
cp ../resources/p_ANCHIVtest_women_by_age.csv HPC_packages/$target_dir/resources/
cp ../resources/treatment_cascade.csv HPC_packages/$target_dir/resources/
## Copy latest WUENIC vaccination coverage:
cp ../resources/WUENIC2025_HepB* HPC_packages/$target_dir/resources/
## Copy the "Intervention_parameters" file:
cp /mnt/c/Users/mpickles/OneDrive\ -\ Imperial\ College\ London/Dropbox_copy/Hepatitis\ B/Modelling/Intervention_parameters.xlsx HPC_packages/$target_dir/resources/
## Copy list of ISOs (countries):
cp ../resources/ListOfISOs.txt ../src/HPC_packages/$target_dir/resources
