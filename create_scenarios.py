## Code creates scenarios...
import sys, os
import glob

## Put BD first as it's the only one that has complicated transitions.
header = "BD,HepB3,PAP,ANCscreen,BirthCohortScreen,CommunityScreen,InFacilityOptOutScreen,DecentralisedServices,TxEligibility,LongActingTx,FunctionalCure"
n_interventions = len(header.split(",")) 

## SQ scenario is dealt with separately:
## "HepB3_WUENIC2025","BD_WUENIC2025""PAP_SQ","ANCScreen_No","BirthCohortScreen_No","CommunityScreen_No","InFacilityOptOutScreen_No","DecentralisedServices_No","TxEligibility_Current","TxTDF","FunctCure_No"
SQ_scenario = ','.join(["-1"]*n_interventions)

continued_improvement_string = "0,0,0,0,0,0,0,0,0,0,0"
base_scenario_string = continued_improvement_string

########################################################################
## These are to keep track of which intervention is which.
## Currently not running universal PAP (PAP for all pregnant women).

possible_values = {"HepB3":["HepB3_ContImp","HepB3WHOtarget"],
                   "PAP":["PAP_SQ","PAP_HVL_ContImp","PAP_HVL_targeted","PAP_PoC"],
                   "ANCscreen":["ANCScreen_No","ANCScreen_Yes"],
                   "BirthCohortScreen":["BirthCohortScreen_No","BirthCohortScreen_Yes"],
                   "CommunityScreen":["CommunityScreen_No","CommunityScreen_Yes"],
                   "InFacilityOptOutScreen":["InFacilityOptOutScreen_No","InFacilityOptOutScreen_Yes"],
                   "DecentralisedServices":["DecentralisedServices_No","DecentralisedServices_Yes"],
                   "TxEligibility":["TxEligibility_Current","TxEligibility_PoC","TxEligibility_Universal"],
                   "LongActingTx":["TxTDF","TxLA"],
                   "FunctionalCure":["FunctCure_No","FunctCure_BepiLike","FunctCure_Future"]}

BD_values = ["BD_ContImp","BD_IFexpansion","BD_OOFexpansion","BD_IF_OOFexpansion"]

BD_allowed_transitions = {"BD_ContImp":["BD_IFexpansion","BD_OOFexpansion"],
                          "BD_IFexpansion":["BD_IF_OOFexpansion"],
                          "BD_OOFexpansion":["BD_IF_OOFexpansion"],
                          "BD_IF_OOFexpansion":[]}




#####################################################################
## Functions:
def check_base_scenario_string_format(scenario_string,header):
    split_string = scenario_string.split(",")
    if(len(split_string)!=len(header.split(","))):
       print("Error - base_scenario_string length is ",len(split_string),". Expected length is ",len(header.split(",")))
       sys.exit(1)
    for c in split_string:
        if not(c.isdigit()):
            print("Error - in scenario_string there is a non-numeric character",c)
            sys.exit(1)
    
    
### Part-written code - want to check we have the right ISOs.
def check_ISOs_to_run():
    main_script_file = open("src/main_script.m","r")
    main_script = main_script_file.readlines()
    main_script_file.close()

    run_on_cluster_line = [line for line in main_script if line.startswith('RUN_ON_CLUSTER =')]
    if(len(run_on_cluster_line)!=1):
        print("Error - wrong number of lines in main_script.m beginning with 'RUN_ON_CLUSTER = '. Exiting\n")
        sys.exit(1)
    
    if run_on_cluster_line[0][0:19]== "RUN_ON_CLUSTER = 1;":
        print("CLUSTER")
    elif run_on_cluster_line[0][0:19]== "RUN_ON_CLUSTER = 0;":        
        print("NOT CLUSTER")
        countries_to_run_line = [line for line in main_script if line.lstrip().startswith('countries_to_run = ')]
        
        print(countries_to_run_line)
    else:
        print("Error - RUN_ON_CLUSTER takes value ",run_on_cluster_line[0][17:18])
        sys.exit(1)
        




##############################################################################
##############################################################################
## Main code:
##############################################################################
##############################################################################


current_scenario_files = glob.glob("current_scenario_???.csv")


if current_scenario_files:
    if(len(current_scenario_files)>1):
        print("Error - code currently can only handle one file called current_scenario_[ISO].csv. The following are in the directory:",current_scenario_files)
        sys.exit(1)
    else:
        current_scenario_file = open(current_scenario_files[0],"r")
        base_scenario_string = current_scenario_file.read().rstrip()
        current_scenario_file.close()
        check_base_scenario_string_format(base_scenario_string,header)
else:
    base_scenario_string = continued_improvement_string
    ##if(len(sys.argv)>1):
    ##base_scenario_string = "0,0,0,0,0,0,0,0,0,0,0"

##print("Current scenario",base_scenario_string)
#check_ISOs_to_run()




previous_iteration = ""

base_scenario = base_scenario_string.rstrip().split(",")

# for i in range(n_interventions):    
#     ## BD
#     if(i==0):
        
if(base_scenario_string==continued_improvement_string):
    new_scenarios = continued_improvement_string+"\n"
else:
    new_scenarios = ""
    
list_of_interventions = header.split(",")
for i in range(n_interventions):

    ## Deal with BD separately (currently it is the only intervention where different combinations are allowed). PAP has multiple states but you just pick the best individual one; BD can be IF, OOF or IF+OOF.
    if list_of_interventions[i]=="BD":
        current_bd_index = int(base_scenario[i])
        current_bd = BD_values[current_bd_index]
        new_bd_names = BD_allowed_transitions[current_bd]
        new_bd_values = [BD_values.index(n) for n in new_bd_names]
        for j in new_bd_values:
            this_scenario_list = base_scenario[:]
            this_scenario_list[i] = str(j)
            this_scenario = ",".join(this_scenario_list) + "\n"
            new_scenarios += this_scenario
    else:
        
        current_intervention_value = int(base_scenario[i])
        if(current_intervention_value==0):
            current_intervention = list_of_interventions[i]
            n_new_interventions = len(possible_values[current_intervention])
            for j in range(1,n_new_interventions):
                this_scenario_list = base_scenario[:]
                this_scenario_list[i] = str(j)
                this_scenario = ",".join(this_scenario_list) + "\n"
                new_scenarios += this_scenario
#print("New scenarios:")                
#print(new_scenarios)

new_scenario_file = open("new_scenarios.csv","w")
new_scenario_file.write(new_scenarios)
new_scenario_file.close()

    ##"BD,HepB3,PAP,ANCscreen,BirthCohortScreen,CommunityScreen,InFacilityOptOutScreen,DecentralisedServices,TxEligibility,LongActingTx,FunctionalCure"


















# i_PAP_SQ = "300"           ## No PAP unless already present.
# i_PAP_HVL_CI = "301"  ## Increase PAP to 20% of HVL unless already present.
# i_PAP_HVL_targeted = "302" ## Increase PAP to X% of HVL unless already present (X% placeholder but will be based on VL testing availability).
# i_PAP_PoC = "303"          ## Eligibility based on PoC test with given sensitivity, specificity and (globally constant) coverage.
# ##i_PAP_all = "304"        ## All eligible (no VL criteria). Uptake to national ANC-1 value.
# allowed_transitions_PAP = {i_PAP_HVL_CI:[i_PAP_HVL_targeted,i_PAP_PoC]}
                           
# i_treat_SQ = "400"          ## Current treatment (Capped at most recent treatemnt data - currently Polaris 2025).
# i_treat_CI = "401"
# i_treat_IFscreening = "402"
# i_treat_PoC = "403"
# ##i_treat_universal
# i_treat_LA = "404"
# i_treat_decentralised = "405"
# i_treat_IFscreening_PoC = "406"
# i_treat_IFscreening_LA = "407"
# i_treat_IFscreening_decentralised = "408"
# i_treat_PoC_LA = "409"
# i_treat_PoC_decentralised = "410"
# i_treat_LA_decentralised = "411"

# i_treat_IFscreening_PoC_LA = "412"
# i_treat_IFscreening_PoC_decentralised = "413"

# allowed_transitions_Treatment = {i_treat_CI:[i_treat_IFscreening,i_treat_PoC,i_treat_LA,i_treat_decentralised],
#                                  i_treat_IFscreening:[i_treat_IFscreening_PoC,i_treat_IFscreening_LA,i_treat_IFscreening_decentralised],
#                                  i_treat_PoC:[i_treat_IFscreening_PoC,i_treat_PoC_LA,i_treat_PoC_decentralised],
#                                  i_treat_LA:[i_treat_IFscreening_LA,i_treat_PoC_LA,i_treat_LA_decentralised],
#                                  i_treat_decentralised:[i_treat_IFscreening_decentralised,i_treat_PoC_decentralised,i_treat_LA_decentralised],
                                 
#                                  }
                                 

# i_functcure_NoCure = "501"
# i_functcure_Bepilike = "502"
# i_functcure_future_cure = "503"
# allowed_transitions_FunctCure = {i_functcure_NoCure:[i_functcure_Bepilike,i_functcure_future_cure]}



# ##allowed_transitions_Screening

# def print_error(tool,current_state,newvalue,allowed_transitions):
#     if(allowed_transitions=="None"):
#         print("Error: for",tool,"in current state",current_state," no transitions are allowed.",current_state,"->",newvalue,"is not allowed. Exiting")
#     else:    
#         print("Error: for",tool,"in current state",current_state," the allowed new intervention states are",allowed_transitions,".",current_state,"->",newvalue,"is not allowed. Exiting")
#     sys.exit(1)
          

# ## Given a current portfolio string add a specified new intervention and return an updated portfolio string:
# def add_one_intervention_to_portfolio(current_portfolio,intervention_to_add,newvalue):
#     current_intervention_list = current_portfolio.split(",")
#     if intervention_to_add=="HepB3":
#         ## Pull out the "0" state (HepB3)
#         current_state = current_intervention_list[0]
#         if(current_state in allowed_transitions_HepB3.keys()):
#             allowed_transitions = allowed_transitions_HepB3[current_state]
#             if(newvalue in allowed_transitions):
#                 current_intervention_list[0]=newvalue
#                 new_portfolio = ",".join(current_intervention_list)
#             else:
#                 print_error("HepB3",current_state,newvalue,allowed_transitions)
#         else:
#             print_error("HepB3",current_state,newvalue,"None")
            
#             return new_portfolio
#     ####
#     elif intervention_to_add=="BD":
#         ## Pull out the "1" state (BD)
#         current_state = current_intervention_list[1]
#         if(current_state in allowed_transitions_BD.keys()):
#             allowed_transitions = allowed_transitions_BD[current_state]
#             if(newvalue in allowed_transitions):
#                 current_intervention_list[0]=newvalue
#                 new_portfolio = ",".join(current_intervention_list)
#             else:
#                 print_error("BD",current_state,newvalue,allowed_transitions)
#         else:
#             print_error("BD",current_state,newvalue,"None")
#     ####
#     elif intervention_to_add=="PAP":
#         ## Pull out the "2" state (PAP)
#         current_state = current_intervention_list[2]
#         if(current_state in allowed_transitions_PAP.keys()):
#             allowed_transitions = allowed_transitions_PAP[current_state]
#             if(newvalue in allowed_transitions):
#                 current_intervention_list[0]=newvalue
#                 new_portfolio = ",".join(current_intervention_list)
#             else:
#                 print_error("PAP",current_state,newvalue,allowed_transitions)
#         else:
#             print_error("PAP",current_state,newvalue,"None")
#     ####
#     elif intervention_to_add=="Treatment":
#         ## Pull out the "3" state (Treatment)
#         current_state = current_intervention_list[3]
#         if(current_state in allowed_transitions_Treatment.keys()):
#             allowed_transitions = allowed_transitions_Treatment[current_state]
#             if(newvalue in allowed_transitions):
#                 current_intervention_list[0]=newvalue
#                 new_portfolio = ",".join(current_intervention_list)
#             else:
#                 print_error("Treatment",current_state,newvalue,allowed_transitions)
#         else:
#             print_error("Treatment",current_state,newvalue,"None")
#     ####
#     elif intervention_to_add=="FunctCure":
#         ## Pull out the "4" state (FunctCure)
#         current_state = current_intervention_list[4]
#         if(current_state in allowed_transitions_FunctCure.keys()):
#             allowed_transitions = allowed_transitions_FunctCure[current_state]
#             if(newvalue in allowed_transitions):
#                 current_intervention_list[0]=newvalue
#                 new_portfolio = ",".join(current_intervention_list)
#             else:
#                 print_error("FunctCure",current_state,newvalue,allowed_transitions)
#         else:
#             print_error("FunctCure",current_state,newvalue,"None")
#     ####
#     elif intervention_to_add=="Screening":
#         ## Pull out the "5" state (Screening)
#         current_state = current_intervention_list[5]
#         if(current_state in allowed_transitions_Screening.keys()):
#             allowed_transitions = allowed_transitions_Screening[current_state]
#             if(newvalue in allowed_transitions):
#                 current_intervention_list[0]=newvalue
#                 new_portfolio = ",".join(current_intervention_list)
#             else:
#                 print_error("Screening",current_state,newvalue,allowed_transitions)
#         else:
#             print_error("Screening",current_state,newvalue,"None")
            
#             return new_portfolio

#     sys.exit(1)



# s = add_one_intervention_to_portfolio("1,1,1,1,1,No additional screening\n","BD",i_BD_IF)
# print(s)
# sys.exit(1)







# ###################################################################    
# ## Order of outputs:
# ## HepB3, BD, PAP, Treatment, FunctCure, AdditionalScreening
# if previous_iteration=="START":
#     ## Status quo: Post-2025 all coverage (BD,HepB3,PAP,diagnosis,treatment) kept fixed at 2025 levels. No functional cure or additional screening.
#     i_SQ = "0,0,0,0,0,No additional screening\n"
#     i_CI = "0,0,1,1,1,No additional screening\n" ## Continued improvement - HepB3 and BD remain at current levels. PAP increases to 5% (2026-2030) of HVL. Treatment and diagnosis increase annually by region-specific rate
#     i_CI_plusB3 = "0,0,1,1,1,No additional screening\n"  ## ContImp + HepB3 to 90% (or current coverage if higher) increasing T_INTERVENTION_START_HepB3-T_INTERVENTION_END_HepB3
#     i_CI_plusBD_IF = "0,0,1,1,1,No additional screening\n"  ## ContImp + in-facility (IF) BD optimization: BD to IF ceiling in countries with BD at present
#     i_CI_plusBD_OOF = "0,0,1,1,1,No additional screening\n"  ## ContImp + BD increases to also reach to 60% of non-facility births (BD introduced in non-BD countries).
#     i_CI_plusBD_IF_OOF = "0,0,1,1,1,No additional screening\n"  ## ContImp + combined IF/OOF BD.
#     i_CI_plusPAP_HVL_targeted = "0,0,1,1,1,No additional screening\n"  ## ContImp + X% PAP coverage of high VL (X%=access to VL monitoring). 
#     i_CI_plusPAP_PoC = "0,0,1,1,1,No additional screening\n"  ## ContImp + PoC test to determine PAP eligibility. Assume 90% sensitivity, and 90% specificity.
#     i_CI_plusPAP_all = "0,0,1,1,1,No additional screening\n"  ## ContImp + PAP ANC-1 coverage of all pregnant women (HVL+LVL)
#     i_CI_plusDx_ANCscreening = "0,0,1,1,1,No additional screening\n"  ## ContImp with additional diagnosis among pregnant women during ANC (ANC-1 capped at HIV screening %)
#     i_CI_plusDx_BirthCohort = 11;  ## ContImp with additional diagnosis through birth cohort screening (those born within 5 years of country introduction of HepB3)
#     i_CI_plusDx_IFscreening = 12; ## ContImp with diagnosis through in-facility screening (e.g. acute care). Use WHO 
#     i_CI_plusDx_CommunityScreening = 13; ## ContImp with additional diagnosis through (realistic) community screening similar to PROLIFICA.
#     ##i_CI_plusDx_CommunityScreening_perfect = 14; ## ContImp with additional diagnosis through perfect community screening.
#     i_CI_plusDx_IntegratedServices = 14;  ## CompInt with integrated services (TBD?)
#     i_CI_plusTx_PoCeligibility = 15; ## CompInt with PoC treatment eligibility 
#     i_CI_plusTx_treatall = 16;
#     i_CI_plusTx_LA = 17;
#     i_CI_plusDecentralisedDxTx = 18;
#     i_CI_plusTx_cure_Bepi = 19;
#     i_CI_plusTx_cure_improved = 20;


