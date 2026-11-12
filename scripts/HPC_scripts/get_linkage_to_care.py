import sys, os
import glob, re


job_script_files = glob.glob("./job*.sh.o*")

## Check if all the

check = ".".join(job_script_files[0].split(".")[:-1])
list_of_isos = []
for j in job_script_files:
    j_split = j.split(".")
    job_id_to_check = ".".join(j_split[:-1])
    if(job_id_to_check!=check):
        print("Error")
        sys.exit(1)
    country_iso = int(j_split[-1])
    list_of_isos += [country_iso]
list_of_isos.sort()


#############################################################
## Functions:
#############################################################
def get_model_est(this_line,estimate_type, data_line):
    model_estimate = float(this_line.split(" ")[0])
    target_est = float(this_line.split("=")[-1])
    if(abs(model_estimate-target_est)>0):
        if(abs(model_estimate-target_est)/(target_est+1e-10)>0.01 and abs(model_estimate-target_est)>0.01):
            if((model_estimate==0 and target_est<1e-2) or (target_est==0 and model_estimate<1e-2)):
                pass
            else:
                print("Error - model and target",estimate_type," in 2016 differ: ",model_estimate,target_est,data_line)
            
            ##sys.exit(1)

    return model_estimate

########################
iso = 1


def read_this_job_script(iso, file_tag):
    outstring_this_job = ""
    job_output_file = open(file_tag+"."+str(iso),"r")

    job_output = job_output_file.readlines()
    job_output_file.close()

    for l in job_output:
        res_run_scenario = re.search(r'The "[0-9]+" scenario \(run number [0-9]+\) started at', l)
        if res_run_scenario:
            i_scenario = int(l.split('"')[1])
            i_run = int(l.split(")")[0].split(" ")[-1])
        else:
            res_2016linkage = re.search(r'CHECKS2016:',l)
            if res_2016linkage:
                ##print(i_run,i_scenario,l)
                ## CHECKS2016: care=0.0300 vs target=0.0300; Tx (of elig)=0.0000 vs target=0.0000; Tx (of sAg+)=0.0000 vs target=0.0000 prop_in_care=0.0000
                this_line_care = l.split("care=")[1].split(";")[0]
                ##
                model_care = get_model_est(this_line_care,"care",l)
                this_line_Txofelig = l.split("Tx (of elig)=")[1].split(";")[0]
                model_Txofelig = get_model_est(this_line_Txofelig,"Txofelig",l)
                ##
                this_line_TxofsAg = l.split("Tx (of sAg+)=")[1].split("prop_in_care=")[0]
                model_TxofsAg = get_model_est(this_line_TxofsAg,"TxofsAg",l)
                prop_in_care2016 = float(l.split("prop_in_care=")[1])

            res_2026linkage = re.search(r'PROPORTION_REMAIN_IN_CARE_2026',l)
            if res_2026linkage:
                prop_linking_2026 = l.split("=")[-1]
                outstring_this_job += str(iso)+","+str(i_run)+","+str(i_scenario)+","+str(model_care)+","+str(model_Txofelig)+","+str(model_TxofsAg)+","+str(prop_in_care2016)+","+str(prop_linking_2026)

    return outstring_this_job
#print(job_output[0:3])

outstring = "ISO,i_run,i_scenario,,prop_eligible_treatedby2016,prop_sAgpos_treatedby2016,prop_link_care2016\n"
for iso in list_of_isos:
    outstring += read_this_job_script(iso,check)
print(outstring)
