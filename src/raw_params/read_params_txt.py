## Script will take the parameters for each country
import sys
import numpy as np

## List of files to read:
# country_s_e_HCCdeaths.txt
# stochas_params_mat.csv
# BD_table.csv
# params.txt
# treatment_rates.txt
# HepB3_table.csv
# dwvec.csv
# ListOfISOs.txt
# num_in_treatment_2016.txt

## TO DO:
# country_s_e_HCCdeaths.txt
# stochas_params_mat.csv
# params.txt

len_sexratio_params = 212
len_migration_params = 212
n_age_gps = 101 ## (0--0, 1--1, 2--2,..., 98--98, 99--99, 100--100)
n_years_totalpop = 152 ##  (1950 to 2101)
#############################################################
ListOfISOs_file = open("ListOfISOs.txt","r")
ListOfISOs_raw = ListOfISOs_file.read().rstrip().splitlines()
ListOfISOs = []
for raw_ISO in ListOfISOs_raw:
    ListOfISOs += [raw_ISO.split(":")[1].rstrip().lstrip()]
##ListOfISOs = ListOfISOs[:2]    
##print(ListOfISOs)

#############################################################
num_in_treatment_2016_file = open("num_in_treatment_2016.txt","r")
num_in_treatment_2016_raw = num_in_treatment_2016_file.read().rstrip().splitlines()
num_in_treatment_2016_dict = {}
for line in num_in_treatment_2016_raw:
    line_split = line.rstrip().lstrip().split(":")
    num_in_treatment_2016_dict[line_split[0]] = line_split[1].lstrip()
if(not(all(key in num_in_treatment_2016_dict.keys() for key in ListOfISOs))):
    print("No match between ListOfISOs and countries in num_in_treatment_2016.txt. Exiting")
    sys.exit(1)
##print(num_in_treatment_2016_dict)
##sys.exit(1)


#############################################################
# BD_table.csv
BD_table_file = open("BD_table.csv","r")
BD_table_raw = BD_table_file.read().rstrip().splitlines()
BD_table_dict = {}
BD_table_header = BD_table_raw[0].rstrip().replace("x","BirthDoseCoverage_").split(",")[1:] ## First column is ISO code:
#print(BD_table_header)
for line in BD_table_raw[1:]:
    BD_table_country_data = line.rstrip().split(",")
    if(BD_table_country_data[0]) in ListOfISOs:
        BD_table_dict[BD_table_country_data[0]] = BD_table_country_data[1:]
    ##else:    ## BD_table.csv has 194 countries so ignore this
        ##print("Unknown ISO ",BD_table_country_data[0])
        ##sys.exit(1)

if(not(all(key in BD_table_dict.keys() for key in ListOfISOs))):
    print("No match between ListOfISOs and countries in BD_table.csv. Exiting")
    sys.exit(1)
    

#############################################################
# HepB3_table.csv
HepB3_table_file = open("HepB3_table.csv","r")
HepB3_table_raw = HepB3_table_file.read().rstrip().splitlines()
HepB3_table_dict = {}
HepB3_table_header = HepB3_table_raw[0].rstrip().replace("x","B3Coverage_").split(",")[1:] ## First column is ISO code:
##print(HepB3_table_header)
for line in HepB3_table_raw[1:]:
    HepB3_table_country_data = line.rstrip().split(",")
    if(HepB3_table_country_data[0]) in ListOfISOs:
        HepB3_table_dict[HepB3_table_country_data[0]] = HepB3_table_country_data[1:]
    ##else:   ## HepB3_table.csv has 194 countries so ignore this:
        ##print("Unknown ISO ",HepB3_table_country_data[0])
        ##sys.exit(1)

if(not(all(key in HepB3_table_dict.keys() for key in ListOfISOs))):
    print("No match between ListOfISOs and countries in HepB3_table.csv. Exiting")
    sys.exit(1)
    ##print(BD_table_country_data)
##print(HepB3_table_dict)

    


#############################################################
treatment_rates_file = open("treatment_rates.txt","r")
treatment_rates_raw = treatment_rates_file.read().rstrip().splitlines()

treatment_rates_dict = {}


i=0
while(i<len(treatment_rates_raw)):
    line = treatment_rates_raw[i].rstrip()    
    country_iso = line[:4]
    treatment_rates_dict[country_iso] = []
    ##print(country_iso)
    check = line[3:]
    if(not(check==":[")):
        print("Error - not an ISO",line)
        sys.exit(1)
    i += 1
    while not(treatment_rates_raw[i]=="]"):
        line = treatment_rates_raw[i].rstrip()
        ##print(line,i)
        num  = line.rstrip(",")
        treatment_rates_dict[country_iso] += [num]
        i += 1
    i += 1
    #print(i,treatment_rates_dict)



    ## treat_coverage_in_2016 (in HBVmodel.m) is HBsAg_treat_cov_all_ages (in country_level_analysis.m) and pop_size_HBsAg_treatment_2016_vec(4) which is the country-level of pop_size_HBsAg_treatment_map(ISO).
    

# AFG:[
#   2016,
#   0,
#   0.036911305717558185,
#   1,
#   0.13968763412904533,
#   1
# ]

########################################################################################
## params.txt:
## This is a pain because it's so big! Here are the fields:
# "country_name": "Afghanistan",
# "Pop_byAgeGroups_1950"
# "total_pop_female"
# "total_pop_male"
# "fert"
# "sex_ratios"
# "net_migration"
# "MortalityRate_Women"
# "MortalityRate_Men"
# "CancerDeathRate": 0.5,
# "ClearanceRateWomenCoFactor": 1,
# "Efficacy_BirthDoseVacc_HbSAg": 0.95,
# "Efficacy_BirthDoseVacc_HbEAg": 0.83,
# "Efficacy_InfantVacc": 0.95,
# "p_VerticalTransmission_HbSAg_NoIntv": 0.17050716633937812,
# "p_VerticalTransmission_HbEAg_NoIntv": 0.9,
# "beta_1to15"
# "beta_5plus"
# "ReducInTransmission": 0,
# "YearReducInTransmission": 2100

params_file = open("params.txt","r")
##params_file = open("temp1.txt","r")
## Creates an array where each entry is a country:
params_raw = params_file.read().rstrip().split(" COUNTRY: ")[1:]
## Note that we need to remove the first entry of the above after splitting (as it's the bit before the first " COUNTRY: ", which is nothing)


# def cut_country_params(country_params_raw,s):
#     l = len(s)
#     i_temp = country_params_raw.index(s)
#     ## Add 1 to remove newline:
#     country_params_raw = country_params_raw[(i_temp+l+1):]
#     return country_params_raw

def cut_country_params(country_params_raw, start_string, end_string):
    l_start = len(start_string)
    l_end = len(end_string)
    i_start = country_params_raw.index(start_string)
    i_end = country_params_raw.index(end_string)
    ## Add 1 to remove newline:
    this_param_block = country_params_raw[(i_start+l_start+1):i_end]
    country_params_raw = country_params_raw[(i_end+l_end+1):]
    return [this_param_block, country_params_raw]

def process_pop1950(pop1950_param_block):
    age_gp_sep = '''
    ],
    [
      '''
    ## Each array entry of pop1950_param_byage_eachsex is F + M pops for a given age group:
    pop1950_param_byage_eachsex = pop1950_param_block.split(age_gp_sep)
    pop1950_byage_F = []
    pop1950_byage_M = []
    for x in pop1950_param_byage_eachsex:
        pop_thisagegp_by_sex = [y.lstrip().rstrip() for y in x.split(",")]
        pop1950_byage_F += [pop_thisagegp_by_sex[0]]
        pop1950_byage_M += [pop_thisagegp_by_sex[1]]
    return [pop1950_byage_F, pop1950_byage_M]

def process_total_pop(totalpop_param_block):
    ## total_pop_F/M are two 101 x 152 matrix of 152 years (1950 to 2101) for 101 age groups (0--0, 1--1, 2--2,..., 98--98, 99--99, 100--100)
    time_sep = '''
    ],
    [
'''

    totalpop_param_block = totalpop_param_block.lstrip("[")
    totalpop_param_bytime = [x.lstrip().rstrip() for x in totalpop_param_block.split(time_sep)]
    if(not(len(totalpop_param_bytime)==n_age_gps)):
        print("Error in process_total_pop - number of age groups wrong. Exiting")
        sys.exit(1)
    
    # The "-1" is so it's obvious if I don't fill in everything.
    totalpop_matrix = -1 * np.ones((n_age_gps, n_years_totalpop))

    for i in range(len(totalpop_param_bytime)):
        totalpop_overtime_forthisage = totalpop_param_bytime[i]
        totalpop_overtime_forthisage_numeric = [int(float(x.lstrip().rstrip())) for x in totalpop_overtime_forthisage.split(",")]
        for j in range(len(totalpop_overtime_forthisage_numeric)):
            totalpop_matrix[i][j] = totalpop_overtime_forthisage_numeric[j]
        
    return totalpop_matrix

## Can process total_pop_female_param_block, total_pop_male_param_block
## which are two 101 x 152 matrix of 152 years (1950 to 2101) for 101 age groups (0--0, 1--1, 2--2,..., 98--98, 99--99, 100--100)
## and fert (which is a 212 x 100 (I think))
## totalpop_matrix = -1 * np.ones((n_age_gps, n_years_totalpop))
def process_block_overtime(param_block, size_x, size_y, vartype):
    time_sep = '''
    ],
    [
'''

    param_block = param_block.lstrip("[")
    param_block_bytime = [x.lstrip().rstrip() for x in param_block.split(time_sep)]
    if(not(len(param_block_bytime)==size_x)):
        print("Error in process_block_overtime() - number of x groups",len(param_block_bytime)," wrong. Exiting")
        sys.exit(1)
    
    # The "-1" is so it's obvious if I don't fill in everything.
    processed_matrix = -1 * np.ones((size_x, size_y))

    for i in range(len(param_block_bytime)):
        param_overy_forthisx = param_block_bytime[i]
        if(vartype=="int"):
            param_overy_forthisx_numeric = [int(float(x.lstrip().rstrip())) for x in param_overy_forthisx.split(",")]
        elif(vartype=="float"):
            param_overy_forthisx_numeric = [float(x.lstrip().rstrip()) for x in param_overy_forthisx.split(",")]
        else:
            print("Unknown var type in process_block_overtime(). Exiting")
            sys.exit(1)
        
        for j in range(len(param_overy_forthisx_numeric)):
            processed_matrix[i][j] = param_overy_forthisx_numeric[j]

    ## Check dimensions of array:
    ##print(processed_matrix.shape)
    ## Check no elements are still "-1":
    check = np.argwhere(processed_matrix == -1)
    if len(check)>0:
        print("Error in process_block_overtime(): some elements of matrix are still -1",check)
        sys.exit(1)

    return processed_matrix


# def process_fertility(fert_param_block):
#     if not(len(fert_param_block) == 140):
#         print("Fertility block length is ",len(fert_param_block),". Exiting\n")
#         ##print(fert_param_block)
#         sys.exit(1)

#     return 1

## Fertility is really big, so this will make it smaller:
## Take a 1000 x 212 matrix and turn it into a 20 x 212 matrix
def rationalise_fertility_params(fert_alltimesteps):

    fertility_byageband_alltimesteps = -1 * np.ones((20, 212))
    for j in range(212):
        ## Extract one row and reshape into a 20x50 matrix (each row is one 5-year age band)
        fert_thisyear = fert_alltimesteps[:,j].reshape(20,50)
        ## Check that all fertility rates in a given 5 year age band are the same:
        for k in range(20):
            if not(len(set(fert_thisyear[k,:]))):
                print("Error - cannot compress fertility age bands in rationalise_fertility_params(). Exiting")
                sys.exit(1)
            fertility_byageband_alltimesteps[k,j] = fert_thisyear[k,0]
        #print(fert_thisyear)
    
    ## Check - I think UNPD only have rates >1950 so the first 60 rows should be duplicates:
    for j in range(1,60):
        if(not(all(fertility_byageband_alltimesteps[:,j]==fertility_byageband_alltimesteps[:,0]))):
            print("Error in rationalise_fertility_params() - the fertility params from 1890-1949 differ unexpectedly. Exiting")
            sys.exit(1)
        ##print(j+1890,all(fertility_byageband_alltimesteps[:,j]==fertility_byageband_alltimesteps[:,0]))
    fertility_byageband_alltimesteps = fertility_byageband_alltimesteps[:,60:]
    ##print("Fert shape:",fertility_byageband_alltimesteps.shape)
    return fertility_byageband_alltimesteps


##
## Mortality has a lot of repetition (1890-1950) so can remove to make it smaller:
## Take a 1000 x 212 matrix and turn it into a 20 x 212 matrix
def rationalise_mortality_params(mortality_alltimesteps):
    ## UNPD only have rates >1950 so the first 60 rows should be duplicates:

    for j in range(1,60):
        if(not(all(mortality_alltimesteps[j,:]==mortality_alltimesteps[0,:]))):
            print("Error in rationalise_mortality_params() - the mortality params from 1890-1949 differ unexpectedly. Exiting", j,mortality_alltimesteps[j,:],mortality_alltimesteps[0,:])
            sys.exit(1)
    mortality_alltimesteps = mortality_alltimesteps[60:,:]

    # Mortality rates change every 5th year:
    ##30
    mortality_alltimesteps_final = -1 * np.ones((31, 21))
    for j in range(30):
        ## 
        if(not(all(mortality_alltimesteps[(j*5),:]==mortality_alltimesteps[(j*5+1),:]))
           or not(all(mortality_alltimesteps[(j*5),:]==mortality_alltimesteps[(j*5+2),:]))
           or not(all(mortality_alltimesteps[(j*5),:]==mortality_alltimesteps[(j*5+3),:]))
           or not(all(mortality_alltimesteps[(j*5),:]==mortality_alltimesteps[(j*5+4),:]))):
            print("Error in rationalise_mortality_params() - the mortality params over the five year period",str(1950+j),"-",str(1954+j)," differ unexpectedly. Exiting", j,mortality_alltimesteps[j,:],mortality_alltimesteps[0,:])
            sys.exit(1)
        mortality_alltimesteps_final[j,:] = mortality_alltimesteps[(j*5),:]

    mortality_alltimesteps_final[30,:] = mortality_alltimesteps[151,:]
    ##mortality_alltimesteps = mortality_alltimesteps[0:5:152,:]

    return mortality_alltimesteps_final


def process_sex_ratio(sexratio_param_block):
    sexratio_params = [s.lstrip().rstrip() for s in sexratio_param_block.split(",")]

    if(len(sexratio_params)!=len_sexratio_params):
        print("Error - length of sex ratio params is not 212. Exiting")
        sys.exit(1)
    return sexratio_params

## Process a 1D parameter block (sexratio, net_migration)
def process_vector(param_block, len_params):
    params_vect = [s.lstrip().rstrip() for s in param_block.split(",")]

    if(len(params_vect)!=len_params):
        print("Error - length of params vector is not ",len_params,". Exiting")
        sys.exit(1)
    return params_vect

    ## This marks the end of a given parameter set (e.g. fertility):

def extract_last_params(params_block):
    ## Remove the trailing whitespace and "}" from params_block; then split into lines and remove whitespace and trailing comma from each line, then stick in an array:
    last_params = [x.lstrip().rstrip(",") for x in params_block.rstrip().rstrip("}").splitlines()]

    list_of_params = ['"CancerDeathRate": ', '"ClearanceRateWomenCoFactor": ', '"Efficacy_BirthDoseVacc_HbSAg": ', '"Efficacy_BirthDoseVacc_HbEAg": ', '"Efficacy_InfantVacc": ', '"p_VerticalTransmission_HbSAg_NoIntv": ', '"p_VerticalTransmission_HbEAg_NoIntv": ', '"beta_1to15": ', '"beta_5plus": ', '"ReducInTransmission": ', '"YearReducInTransmission": ']

    output_string = ""
    for i in range(len(list_of_params)):
        param_line = last_params[i]
        
        if(last_params[i].find(list_of_params[i])==0):
            this_param = last_params[i].replace(list_of_params[i],"")
            output_string += this_param + ","
        else:
            print("Error: Missing parameter",last_params[i])
            sys.exit(1)
    output_string = output_string.rstrip(",")
    return [list_of_params,output_string]
    ##print(last_params)
    
    
##################################################################################
## Main code:
##################################################################################
## This marks the end of a given parameter set (e.g. fertility):
end_str = '''    ]
  ],'''

end_str_sexratio = '''
  ],'''


                                          

sex_ratio_dict = {}
pop1950_byage_F = {}
pop1950_byage_M = {}
total_pop_F = {}
total_pop_M = {}
fert = {}
migration_dict = {}

mortalityF = {}
mortalityM = {}
last_params = {}

for country_params_raw in params_raw:
    country_ISO = country_params_raw[0:3]
    print("Parsing ",country_ISO)
    [pop1950_param_block, country_params_raw] = cut_country_params(country_params_raw,
                                                                   '''"Pop_byAgeGroups_1950": [
    [''',
    end_str)
    [pop1950_byage_F[country_ISO], pop1950_byage_M[country_ISO]] = process_pop1950(pop1950_param_block)
    ##print("pop1950_byage_F:",len(pop1950_byage_F[country_ISO]))

    ### Total pop F (101 x 152 matrix):
    [total_pop_female_param_block, country_params_raw] = cut_country_params(country_params_raw,
                                                                            '''"total_pop_female": [
    [''',end_str)
    ##total_pop_F[country_ISO] = process_total_pop(total_pop_female_param_block)
    total_pop_F[country_ISO] = process_block_overtime(total_pop_female_param_block, n_age_gps, n_years_totalpop, "int")

    
    ### Total pop M (101 x 152 matrix):
    [total_pop_male_param_block, country_params_raw] = cut_country_params(country_params_raw,
                                                                          '''"total_pop_male": [
    [''',end_str)
    ##total_pop_M[country_ISO] = process_total_pop(total_pop_male_param_block)
    total_pop_M[country_ISO] = process_block_overtime(total_pop_male_param_block, n_age_gps, n_years_totalpop, "int")
    ##print("Total_pop size:",total_pop_M[country_ISO].shape)

    ### Fertility (1000 x 212 matrix):    
    [fert_param_block, country_params_raw] = cut_country_params(country_params_raw, '''"fert": [
    [''',end_str)
    fert_alltimesteps = process_block_overtime(fert_param_block, 1000, 212, "float")
    ## The fertility parameters are 1000x212 - the 1000 means "every timestep for all age groups".
    ## UN WPP gives fertility for every 5 year age group, so this is a lot of repetition.
    fert[country_ISO] = rationalise_fertility_params(fert_alltimesteps)

    ##
    
    ##print("Fertility:",fert[country_ISO].shape)

    
    ### Sex ratio (1D vector 212 long - number of years 1890-2101 inclusive):
    [sexratio_param_block, country_params_raw] = cut_country_params(country_params_raw, '"sex_ratios": [',end_str_sexratio)
    sex_ratio_dict[country_ISO] = process_sex_ratio(sexratio_param_block)

    ### Migration (1D vector 212 long - number of years 1890-2101 inclusive):
    [migration_param_block, country_params_raw] = cut_country_params(country_params_raw, '"net_migration": [',end_str_sexratio)
    
    migration_dict[country_ISO] = process_vector(migration_param_block, len_migration_params)
    ##process_vector(param_block, len_params)
    ##print(migration_dict[country_ISO])
    
    ##print("new country_params_raw: ",country_params_raw[0:100])

    ### MortalityRate_Women (21 x 212 matrix):    
    [mortalityF_param_block, country_params_raw] = cut_country_params(country_params_raw, '''"MortalityRate_Women": [
    [''',end_str)
    mortality_alltimesteps = process_block_overtime(mortalityF_param_block, 212, 21, "float")
    mortalityF[country_ISO] = rationalise_mortality_params(mortality_alltimesteps)
    
    ### MortalityRate_Men (21 x 212 matrix):    
    [mortalityM_param_block, country_params_raw] = cut_country_params(country_params_raw, '''"MortalityRate_Men": [
    [''',end_str)
    mortality_alltimesteps = process_block_overtime(mortalityM_param_block, 212, 21, "float")
    mortalityM[country_ISO] = rationalise_mortality_params(mortality_alltimesteps)
    ##print("MortalityM:",mortalityM[country_ISO].shape)

    [list_of_params, last_params[country_ISO]] = extract_last_params(country_params_raw)
    #print(country_params_raw)
    



##params_dict = {}
##print(sex_ratio_dict)
##print(pop1950_byage_F)


# outstring = "Year,"+",".join(ListOfISOs)+"\n"
# for i in range(len_sexratio_params):
#     outstring += str(i+1890)+","
#     for country_ISO in ListOfISOs:
#         outstring += sex_ratio_dict[country_ISO][i]+","
#     outstring = outstring.rstrip(",") ## remove trailing comma.
#     outstring += "\n"
# print(outstring)

##############################################################################################
## Generate output files:
##############################################################################################




def make_outputfile_allcountries_byyear(output_dict, start_year, n_years, outfilename):
    outstring = "Year,"+",".join(ListOfISOs)+"\n"
    for i in range(n_years):
        outstring += str(i+start_year)+","
        for country_ISO in ListOfISOs:
            outstring += output_dict[country_ISO][i]+","
        outstring = outstring.rstrip(",") ## remove trailing comma.
        outstring += "\n"
    outfile = open(outfilename,"w")
    outfile.write(outstring)
    outfile.close()


## For fert:
def make_outputfile_onecountry(countryISO, output_dict, age_gp_labels, start_year, n_years, outfilename):
    outstring = "Year,"+",".join(age_gp_labels)+"\n"
    ##print(output_dict[country_ISO].shape)
    ##print(len(age_gp_labels))
    for j in range(n_years):
        outstring += str(j+start_year)+","
        for i in range(len(age_gp_labels)):
            outstring += str(output_dict[country_ISO][i][j])+","
        outstring = outstring.rstrip(",") ## remove trailing comma.
        outstring += "\n"
    outfile = open(outfilename+"_"+countryISO+".csv","w")
    outfile.write(outstring)
    outfile.close()
    ##return(outstring)


def make_outputfile_onecountry_combinesexes(countryISO, output_dictF, output_dictM, age_gp_labels, year_labels, outfilename):
    outstring = "Year,"+",".join([str(x)+"F" for x in age_gp_labels])+",".join([str(x)+"M" for x in age_gp_labels])+"\n"
    for j in range(len(year_labels)):
        outstring += str(year_labels[j])+","
        for i in range(len(age_gp_labels)):
            outstring += str(output_dictF[country_ISO][i][j])+","
        for i in range(len(age_gp_labels)):
            outstring += str(output_dictM[country_ISO][i][j])+","
        outstring = outstring.rstrip(",") ## remove trailing comma.
        outstring += "\n"
    outfile = open(outfilename+"_"+countryISO+".csv","w")
    outfile.write(outstring)
    outfile.close()
    ##return(outstring)


## Version for mortality as transposed:
def make_outputfile_onecountry_combinesexes_transpose(countryISO, output_dictF, output_dictM, age_gp_labels, year_labels, outfilename):
    outstring = "Year,"+",".join([str(x)+"F" for x in age_gp_labels])+",".join([str(x)+"M" for x in age_gp_labels])+"\n"
    for i in range(len(year_labels)):
        outstring += str(year_labels[i])+","
        for j in range(len(age_gp_labels)):
            outstring += str(output_dictF[country_ISO][i][j])+","
        for j in range(len(age_gp_labels)):
            outstring += str(output_dictM[country_ISO][i][j])+","
        outstring = outstring.rstrip(",") ## remove trailing comma.
        outstring += "\n"
    outfile = open(outfilename+"_"+countryISO+".csv","w")
    outfile.write(outstring)
    outfile.close()
    ##return(outstring)

##############################################################################
## Now write the files:
##############################################################################


## "Pop_byAgeGroups_1950" - needs to be written 1 file (combine F+M)
outstring = "Year,"+",".join(ListOfISOs)+"\n"
## 21 age groups:
age_gps = ["0","1-4"]+[str(x)+"-"+str(x+4) for x in range(5,100,5)]
for i in range(len(age_gps)):
    outstring += age_gps[i]+"F,"
    for country_ISO in ListOfISOs:
        outstring += pop1950_byage_F[country_ISO][i]+","
    outstring = outstring.rstrip(",") ## remove trailing comma.
    outstring += "\n"
for i in range(len(age_gps)):
    outstring += age_gps[i]+"M,"
    for country_ISO in ListOfISOs:
        outstring += pop1950_byage_M[country_ISO][i]+","
    outstring = outstring.rstrip(",") ## remove trailing comma.
    outstring += "\n"
outfile = open("processed_params/params_pop1950.csv","w")
outfile.write(outstring)
outfile.close()

    


year_labels = [str(x) for x in range(1950, (1950+152))]

total_pop_age_labels = range(101)
for country_ISO in ListOfISOs:
    make_outputfile_onecountry_combinesexes(country_ISO, total_pop_F, total_pop_M, total_pop_age_labels, year_labels, "processed_params/TotalPop")


fert_age_range_labels = [str(x)+"-"+str(x+4) for x in range(0,100,5)]
for country_ISO in ListOfISOs:
    ##make_outputfile_onecountry(country_ISO, fert, fert_age_range_labels, 1890, 212, "Fertility")
    ## rationalise_fertility_params() reduces the years from 1890-2101 to 1950-2101
    make_outputfile_onecountry(country_ISO, fert, fert_age_range_labels, 1950, 152, "processed_params/Fertility")

make_outputfile_allcountries_byyear(sex_ratio_dict, 1890, len_sexratio_params, "processed_params/param_sexratio.csv")
make_outputfile_allcountries_byyear(migration_dict, 1890, 212, "processed_params/param_netmigration.csv")
    
    
### MortalityRate_Women (21 x 212 matrix):    
##mortalityF[country_ISO]
mortality_age_labels = ["0","1-4"]+[str(x)+"-"+str(x+4) for x in range(5,100,5)]
mortality_year_labels = [str(x)+"-"+str(x+4) for x in range(1950,(1950+152),5)]

for country_ISO in ListOfISOs:
    make_outputfile_onecountry_combinesexes_transpose(country_ISO, mortalityF, mortalityM, mortality_age_labels, mortality_year_labels, "processed_params/Mortality")

list_of_params = ",".join(list_of_params).replace('"','').replace(':','').replace(' ','')
outstring_lastparams = list_of_params + "\n"
for country_ISO in ListOfISOs:
    outstring_lastparams += country_ISO + "," + last_params[country_ISO] + "\n"
outstring_lastparams = outstring_lastparams.rstrip()
outfile = open("processed_params/params_bycountry.csv","w")
outfile.write(outstring_lastparams)
outfile.close()
# "Pop_byAgeGroups_1950"
# "total_pop_female"
# "total_pop_male"
# "fert"
# "sex_ratios"
# "net_migration"
# "MortalityRate_Women"
# "MortalityRate_Men"
# "CancerDeathRate": 0.5,
# "ClearanceRateWomenCoFactor": 1,
# "Efficacy_BirthDoseVacc_HbSAg": 0.95,
# "Efficacy_BirthDoseVacc_HbEAg": 0.83,
# "Efficacy_InfantVacc": 0.95,
# "p_VerticalTransmission_HbSAg_NoIntv": 0.17050716633937812,
# "p_VerticalTransmission_HbEAg_NoIntv": 0.9,
# "beta_1to15"
# "beta_5plus"
# "ReducInTransmission": 0,
# "YearReducInTransmission": 2100
