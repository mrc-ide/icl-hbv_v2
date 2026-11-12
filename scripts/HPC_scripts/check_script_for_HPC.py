## Script makes sure the code should run on the HPC:
## Check all files in place.
## Makes sure the key MatLab code script constants are set correctly.
## Check that "N_SCENARIOS = " in process_results_HPCv4.py matches the number of scenarios in the MatLab code.
## Check for job_HBVmodel    scripts - and check number of countries.
## Make outputs/ folder:

import sys, os
import glob

## Make sure the "outputs" folder exists:
if not os.path.exists("outputs"):
    os.makedirs("outputs")


####################################################################
## Check I have all the necessary files:
needed_src_files = ["HBVmodel.m","country_level_analyses.m","get_treatment_eligible_ageindices.m","make_daly_mat.m","assign_PAP_VL_params.m","get_range_of_FracSPosHighVL.m","main_script.m","safe_greater_or_equal_to.m"]
needed_resources_files = ["ANC_coverage_data.csv","WUENIC2025_HepB3_vaccination_coverage.xlsx","stochastic_parameters_mat.mat","GHO_in_facility_births_summary.csv","WUENIC2025_HepB_BD_vaccination_coverage.xlsx","treatment_2016_map.mat","Intervention_parameters.xlsx","p_ANCHIVtest_women_by_age.csv","treatment_cascade.csv","ListOfISOs.mat","params_map.mat","treatment_rates_map.mat","ListOfISOs.txt","particle_index_countries.mat","vaccination_coverages.mat","WHO_region_map.mat","reference_prevalences.mat"]

current_resources_files = os.listdir(path='resources/')
for f in needed_resources_files:
    if not(f in current_resources_files):
        print("***ERROR***: Missing file",f," in resources/")
        ##sys.exit(1)
current_src_files = os.listdir(path='src/')
for f in needed_src_files:
    if not(f in current_src_files):
        print("***ERROR***: Missing file",f," in src/")
        ##sys.exit(1)


#########################################################################
## Check main_script.m for the following:
## 1. RUN_ON_CLUSTER = 1;
## 2. num_stochas_runs (normally 200 - print a warning if less)

main_script_file = open("src/main_script.m","r")
main_script = main_script_file.readlines()
main_script_file.close()

RUN_ON_CLUSTER = 0
num_stochas_runs = 0
for line in main_script:
    if line[0:16]== "RUN_ON_CLUSTER =":
        if line[0:19]== "RUN_ON_CLUSTER = 1;":
            if(RUN_ON_CLUSTER==0):
                RUN_ON_CLUSTER=1
            else:
                print("***ERROR*** - more than one RUN_ON_CLUSTER statement. Exiting")
                sys.exit(1)
        else:
            print("***ERROR*** - unexpected statement for RUN_ON_CLUSTER:",line,"Exiting.")
            sys.exit(1)
    if line[0:18]== "num_stochas_runs =":
        if line[0:24]== "num_stochas_runs = 200;":
            if(num_stochas_runs==0):
                num_stochas_runs=1
            else:
                print("***ERROR*** - more than one RUN_ON_CLUSTER statement. Exiting")
                sys.exit(1)
        else:
            print("Warning: num_stochas_runs =",int(line.split(" = ")[1].split(";")[0]))

if(RUN_ON_CLUSTER==0):
    print("***ERROR*** - No RUN_ON_CLUSTER statement. Exiting\n")
    sys.exit(1)

########################################################
## Check country_level_analyses.m file:
country_level_analyses_file = open("src/country_level_analyses.m","r")
country_level_analyses = country_level_analyses_file.readlines()
country_level_analyses_file.close()

NUM_SCENARIOS_SET = 0
# num_scenarios = 104;
# %start_scenario = 16;
# start_scenario = 100


for line in country_level_analyses:
    ## Remove annoying whitespace on the left:
    line = line.lstrip() 
    if line[0:5]== "num_s":
        if line[0:15]== "num_scenarios =":
            try:
                n_scenarios = int(line[15:].split(";")[0].lstrip())
            except:
                print("Error - can't get number of scenarios. Exiting")
                sys.exit(1)
            if(NUM_SCENARIOS_SET==0):
                if(n_scenarios!=37):
                    print("Warning: n_scenarios =",n_scenarios)
                NUM_SCENARIOS_SET=1
            else:
                print("***ERROR*** - more than one num_scenarios statement. Exiting")
                sys.exit(1)
        else:
            print("***ERROR*** - unexpected statement for RUN_ON_CLUSTER: %s.\nExiting",line)
            sys.exit(1)
    elif line[0:14]== "start_scenario":
        start_scenario = int(line.split("=")[1].split(";")[0].lstrip())
        if(start_scenario!=1):
            print("***WARNING*** - start_scenario=",start_scenario)
if(NUM_SCENARIOS_SET==0):
    print("***ERROR*** - No num_scenarios statement. Exiting\n")
    sys.exit(1)

###########################################################################
HBVmodel_file = open("src/HBVmodel.m","r")
HBVmodel_script = HBVmodel_file.readlines()
HBVmodel_file.close()


for line in HBVmodel_script:
    if line[0:15]== "PRINT_VERBOSE =":
        if line[0:18]!= "PRINT_VERBOSE = 0;":
            print("Warning - PRINT_VERBOSE non-zero",line[0:18])
    if line[0:15]== "WARN_NEGATIVE =":
        if line[0:18]!= "WARN_NEGATIVE = 0;":
            print("Warning - WARN_NEGATIVE non-zero",line[0:18])
    if line[0:12]== "RUN_CHECKS =":
        if not(line[0:15]== "RUN_CHECKS = 0;"):
            print("Warning - RUN_CHECKS non-zero",line[0:15])
sys.exit(1)



################################################################################
## Look in job script:
job_script_files = glob.glob("./job*.sh")
if(len(job_script_files)==0):
    print("***ERROR*** - no job script present. Exiting\n")
elif(len(job_script_files)>1):
    print("Warning: multiple python scripts. Currently assuming process_results_HPCv4.py is the correct one")

job_script_FILE = open(job_script_files[0])
job_script = job_script_FILE.readlines()
job_script_FILE.close()
PBS_settings = [ele for ele in job_script if ele.startswith("#PBS")]
N_JOBS = [ele for ele in PBS_settings if ele.startswith("#PBS -J ")]
if(len(N_JOBS)==1):
    N_JOBS = N_JOBS[0].split("-J")[1].rstrip().lstrip()
    if(not(N_JOBS=="1-110")):
        print("Warning - running for N_JOBS=",N_JOBS,"countries")


## Tidy up job_script to remove comments and blank lines
job_script = [ele for ele in job_script if not (ele.startswith("#") or ele=="\n")]
#job_script = [ele for ele in job_script if not ele=="\n"]

python_filename = [ele for ele in job_script if ele.startswith("python3")]
if(len(python_filename)==1):
    python_filename = python_filename[0].split()[1]
else:
    print("***ERROR*** - no/multiple python scripts (",python_filename,") listed in ",job_script_file[0])

#################################################################################
## Check the script process_results_HPCv4.py is expecting the same number of scenarios as the MatLab code is generating:

process_results_files = glob.glob("./process_results*.py")
if(len(process_results_files)>1):
    if(python_filename in process_results_files):
        print("Warning: multiple python scripts (",process_results_files,"). Currently assuming ",python_filename," is the correct one")

with open(python_filename,'r') as f:
    targets = [line for line in f if "N_SCENARIOS =" in line]
f.close()
if(len(targets)>1):
    print("***ERROR*** - multiple definitions of N_SCENARIOS in ",python_filename,": ",targets)
N_SCENARIOS_PYTHON = int(targets[0].split("=")[1].lstrip().rstrip())
if(N_SCENARIOS_PYTHON!=n_scenarios):
    print("***ERROR*** - python script ",python_filename,"does not match the number of scenarios in the matlab code country_level_analyses.m: ",N_SCENARIOS_PYTHON,n_scenarios,"Exiting")
    sys.exit(1)
print("=================Finished - check warnings===============")
#print(country_level_analyses[0:4])

