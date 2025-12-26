


load(fullfile(basedir,'resources','ListOfISOs.mat')) % contains ListOfISOs


load("../outputs/results_countries_default_stochastic_run_1.mat")
a = outMap("Status quo infant & BD")

%% a.keys()
a_AFG = a("AFG")
a_AGO = a("AGO")
a_ALB = a("ALB")

##load("../../../../../mpickles/Documents/Hepatitis_B/icl-hbv/outputs/a.mat")
load("../../../../OneDrive - Imperial College London/Dropbox_copy/Hepatits B/Modelling/icl-hbv/outputs/results_countries_default_stochastic_run_1.mat")
b = outMap("Status quo infant & BD")
b_AFG = b("AFG")
b_AGO = b("AGO")
b_ALB = b("ALB")


isequal(a_AFG,b_AFG)
isequal(a_AGO,b_AGO)
isequal(a_ALB,b_ALB)

isequal(a_AFG,b_AFG) & isequal(a_AGO,b_AGO) &  isequal(a_ALB,b_ALB)
