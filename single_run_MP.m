%%

clear


%%num_stochas_runs = 200;
num_stochas_runs = 1;

start_year = 1890;
T0 = 2101-start_year; % must go up to year 2101 to get incidence readings up to year 2100
end_year = start_year + T0;
num_year_divisions = 10;
assert(rem(log10(num_year_divisions),1)==0) % ensure that num_year_divisions is a multiple of 10
% log10(1) = 0, log10(10) = 1, log10(100) = 2, etc.
% rem(x,1) gives 0 if x is an integer
dt = 1/num_year_divisions; % time step



%% Make additions to path as files in different places:
currentFolder = pwd;
pat = fullfile('icl-hbv','src');
if ~endsWith(currentFolder,pat)
    warning(['Please run this script from within the folder ' pat])
end
basedir = fileparts(fileparts(currentFolder)); % path for folder two levels up from script
addpath(fullfile(basedir,'src')) % path for Mike's testing code (why do we have multiple folders?)
addpath(fullfile(basedir,'src','analysis')) % path for function country_level_analyses
addpath(fullfile(basedir,'src','model')) % path for function HBVmodel
%% MatLab suggests addpath is neater (because avoids clutter of multiple paths). But it only makes sense (to me)
%% if I make some standard library functions which I call from many scripts. So I will remove it here.

%sensitivity_analysis_list = {'default','infant_100','treat_medium','treat_high'};
sensitivity_analysis_list = {'default'};
num_sensitivity_analyses = length(sensitivity_analysis_list);



load(fullfile(basedir,'resources','ListOfISOs.mat')) % contains ListOfISOs, a 1x110 cell
%% MP: ListOfISOs doesn't really do anything useful. 
allcountry_ISOs = "";
for country_num = 1:110
    ISO = ListOfISOs{country_num};
    %%thiscountry_ISO = jsonencode(ListOfISOs(ISO),PrettyPrint=true);
    thiscountry_ISO = string(country_num)+ ": "+ISO+newline
    allcountry_ISOs = strcat(allcountry_ISOs,thiscountry_ISO);
end
writelines(allcountry_ISOs,"ListOfISOs.txt");


load(fullfile(basedir,'resources','vaccination_coverages.mat')) % contains BD_table and HepB3_table
writetable(BD_table,"BD_table.csv",'WriteRowNames',true); 
writetable(HepB3_table,"HepB3_table.csv",'WriteRowNames',true); 
%% BD_table and HepB3_table are 194x40 tables. 
%% The 40 refers to years (I think). 194 is countries (including ones not modelled).
%% Note that the rows have the country names but this is suppressed by default when writing using writetable.

load(fullfile(basedir,'resources','reference_prevalences.mat')) % contains country_s_e_HCCdeaths_map
% country_s_e_HCCdeaths_map is a container map with 110 entries indexed by
% country name (e.g. ZWE)
% country_s_e_HCCdeaths_map("ZWE") is a struct with the following fields:
%   source_HBsAg: 'CDA'
%   HBsAg_prevs_year_1: 2016
%   country_HBsAg_prevalences_by_ages_mid_1_young_old: [0.0120 0.0850]
% I believe it's more complicated for a few countries (specifically China):
%   source_HBsAg: 'Cui'
%   HBsAg_prevs_year_1: 1992
%   country_HBsAg_prevalences_by_ages_mid_1_young_old: [0.0967 0.1005]
%   HBsAg_prevs_middle_year_1: [18×2 double]

allcountry_hccdeaths = "";
for country_num = 1:110
    ISO = ListOfISOs{country_num};
    thiscountry_hccdeaths = jsonencode(country_s_e_HCCdeaths_map(ISO),PrettyPrint=true);
    allcountry_hccdeaths = strcat(allcountry_hccdeaths," COUNTRY: ",ISO,thiscountry_hccdeaths);
end
writelines(allcountry_hccdeaths,"country_s_e_HCCdeaths.txt");




load(fullfile(basedir,'resources','params_map.mat')) % contains params_map and dwvec

% dwvec: 1x15 vector, giving DALY weights (1 per model box) - note that one box is death.
% dwvec gets stuck in params in country_level_analyses.m
% in HBVmodel.m dwvec is *only* used in the following line:
% yld_1yr(k, ag, OutputEventNum-1) = sum( state_prev_vec .* demog.dwvec' );
writematrix(dwvec, "dwvec.csv");

% params_map:
% container map with 110 entries indexed by country name (e.g. ZWE). Each
% entry is a structure, e.g. params_map("ZWE") has the following fields:
%                       country_name: 'Zimbabwe'
%                Pop_byAgeGroups_1950: [100×2 double]    
%                    total_pop_female: [101×152 double]
%                      total_pop_male: [101×152 double]
%                                fert: [1000×212 double]
%                          sex_ratios: [1.0300 1.0300 1.0300 1.0300 1.0300 1.0300 1.0300 1.0300 1.0300 1.0300 1.0300 1.0300 1.0300 1.0300 … ] (1×212 double)
%                       net_migration: [-3.0200e-04 -3.0200e-04 -3.0200e-04 -3.0200e-04 -3.0200e-04 -3.0200e-04 -3.0200e-04 -3.0200e-04 … ] (1×212 double)
%                 MortalityRate_Women: [212×21 double]
%                   MortalityRate_Men: [212×21 double]
%                     CancerDeathRate: 0.5000
%          ClearanceRateWomenCoFactor: 1
%        Efficacy_BirthDoseVacc_HbSAg: 0.9500
%        Efficacy_BirthDoseVacc_HbEAg: 0.8300
%                 Efficacy_InfantVacc: 0.9500
% p_VerticalTransmission_HbSAg_NoIntv: 0.0762
% p_VerticalTransmission_HbEAg_NoIntv: 0.9000
%                          beta_1to15: 1.0000e-03
%                          beta_5plus: 1.0000e-03
%                 ReducInTransmission: 0
%             YearReducInTransmission: 2100
allcountry_params = "";
for country_num = 1:110
    ISO = ListOfISOs{country_num};
    thiscountry_params = jsonencode(params_map(ISO),PrettyPrint=true);
    allcountry_params = allcountry_params + " COUNTRY: " + ISO + thiscountry_params + newline;
end
writelines(allcountry_params,"params.txt");





load(fullfile(basedir,'resources','stochastic_parameters_mat.mat')) % contains stochas_params_mat

%% stochas_params_mat is a 200x883 matrix. I think the 200 refers to the number of stochastic runs.
%% MP question: what is the 883? 883 is prime.
writematrix(stochas_params_mat, "stochas_params_mat.csv");

load(fullfile(basedir,'resources','treatment_2016_map.mat')) % contains num_in_treatment_2016_map and pop_size_HBsAg_treatment_map

%% num_in_treatment_2016_map:
%  A container map with 110 entries indexed by country name (e.g. ZWE). Each
% entry is a double (one number!)
allcountry_num_in_treatment_2016 = "";
for country_num = 1:110
    ISO = ListOfISOs{country_num};
    thiscountry_num_in_treatment_2016 = jsonencode(num_in_treatment_2016_map(ISO),PrettyPrint=true);
    thiscountry_num_in_treatment_2016 = ISO + ":" + thiscountry_num_in_treatment_2016 + newline
    allcountry_num_in_treatment_2016 = strcat(allcountry_num_in_treatment_2016,thiscountry_num_in_treatment_2016);
end
writelines(allcountry_num_in_treatment_2016,"raw_params/num_in_treatment_2016.txt");



load(fullfile(basedir,'resources','treatment_rates_map.mat')) % contains treatment_rates_map
allcountry_treatment_rates = "";
for country_num = 1:110
    ISO = ListOfISOs{country_num};
    thiscountry_treatment_rates = jsonencode(treatment_rates_map(ISO),PrettyPrint=true);
    thiscountry_treatment_rates = ISO + ":" + thiscountry_treatment_rates + newline
    allcountry_treatment_rates = strcat(allcountry_treatment_rates,thiscountry_treatment_rates);
end
writelines(allcountry_treatment_rates,"raw_params/treatment_rates.txt");

filename_diaries = 'diary_countries.out';
diary(fullfile(basedir,'outputs',filename_diaries))
% writes a copy of all subsequent keyboard input and the resulting output (except it does not include graphics) to the named file
% if the file already exists, output is appended to the end of the file


% ListOfScenarios = {...
%     'Status quo infant & BD',...
%     'Status quo infant & BD expansion to 25%',...
%     'Status quo infant & BD expansion to 50%',...
%     'Status quo infant & BD expansion to 75%',...
%     'Status quo infant & BD expansion to 90%',...
%     'Status quo infant & BD drop 5 2020',...
%     'Status quo infant & BD drop 10 2020',...
%     'Status quo infant & BD drop 15 2020',...
%     'Status quo infant & BD drop 20 2020',...
%     'Status quo infant & BD delayed expansion 2023 to 2030',...
%     'Status quo infant & BD delayed expansion 2023 to 2033',...
%     'Status quo infant & BD delayed expansion 2025 to 2040'};
% num_scenarios = length(ListOfScenarios);
% assert(num_scenarios==12)


ListOfScenarios = {'Status quo infant & BD'};
num_scenarios = length(ListOfScenarios);



ListOfAllISOs = ListOfISOs;
%%num_all_countries = length(ListOfAllISOs); %% MP: repetition of num_countries
%%assert(num_all_countries==110)
%%num_countries = length(ListOfISOs);
num_countries = 1;



% variables for importing values for the 7 fitted parameters for a particular country-particle combination
country_start_cols = 2:8:(1+num_countries*8); % there are 1 + num_countries*8 + 2 columns in stochas_params_mat


begin_time_run = datetime('now');
disp(['This analysis started at ' datestr(begin_time_run)])


disp(['Number of sensitivity analyses: ' num2str(num_sensitivity_analyses)])
disp(['Number of scenarios: ' num2str(num_scenarios)])
disp(['Number of countries: ' num2str(num_countries)])
disp(['Number of runs per country: ' num2str(num_stochas_runs)])



ages = 0:dt:(100 - dt); % The age group at the beginning of each compartment; 1 x 1000 double; [0 0.1 0.2 ... 99.8 99.9]
num_age_steps = length(ages); % Number of age-categories (assuming max age is 100 years (i.e. oldest person we see at the start of an iteration is aged 100-dt)

theta = 0.01; % proportion of persons newly infected that develop severe symptoms during acute infection
CFR_Acute = 0.05; % proportion of people with severe acute infection that die from severe acute infection
rate_6months = 2; % the rate necessary to move everyone out of States 14 and 15 into State 2 in 6 months

ECofactor = 15; % Multiple for rate of transmission for horizontal transmission with HbEAg vs HbSAg
% epsilon_2 and epsilon_3 in Table S2, p. 6 of the appendix to the publication


Snames = {
    'Susceptible', ... % 1
    'HBV: Immune Tolerant', ... % 2
    'HBV: Immune Reactive', ... % 3
    'HBV: Asymptomatic Carrier', ... % 4
    'HBV: Chronic Hep B', ... % 5
    'HBV: Comp Cirrhosis', ... % 6
    'HBV: Decomp Cirrhosis', ... % 7
    'HBV: Liver Cancer', ... % 8
    'HBV: Immune (Rec. or vacc.)', ... % 9
    'HBV: TDF-Treatment', ... % 10
    'Prematurely dead due to HBV', ... % 11
    '3TC-Treatment', ... % 12
    'Failed 3TC-Treatment', ...  % 13
    'Non-severe acute', ...  % 14
    'Severe acute' ...  % 15
    }; % 1 x 15 cell array
num_states = length(Snames);
assert(num_states==15) % 15 disease states


% risk of chronic carriage due to horizontal and vertical transmission
risk_ages_edmunds_func = @(age) exp(-0.645*age^(0.455));
p_ChronicCarriage = 0.885 * ones(1, num_age_steps, 2, 2);
% does not differ over disease stages
% 4 dimensions necessary so that it can be multiplied by X
risk_p_ChronicCarriage = arrayfun(@(xx) risk_ages_edmunds_func(xx), ages(:,ages>0.5));
risk_p_ChronicCarriage = repmat(risk_p_ChronicCarriage,1,1,2,2);
p_ChronicCarriage(:, ages>0.5, :, :) = risk_p_ChronicCarriage;


sensitivity_analysis_hours_vec = repmat(duration(0,0,0),1,num_sensitivity_analyses);

for sensitivity_analysis_num=1:num_sensitivity_analyses


    sensitivity_analysis = sensitivity_analysis_list{sensitivity_analysis_num};
    disp(['Sensitivity analysis: ' sensitivity_analysis])


    begin_time_sensitivity_analysis = datetime('now');
    disp(['The "' sensitivity_analysis '" sensitivity analysis started at ' datestr(begin_time_sensitivity_analysis)])


    for stochas_run=1:num_stochas_runs

        country_level_analyses(sensitivity_analysis,...
            num2str(stochas_run),num_stochas_runs,...
            ListOfScenarios,num_scenarios,...
            ListOfAllISOs,ListOfISOs,num_countries,...
            BD_table,HepB3_table,...
            num_in_treatment_2016_map,pop_size_HBsAg_treatment_map,treatment_rates_map,...
            country_s_e_HCCdeaths_map,...
            params_map,dwvec,stochas_params_mat,country_start_cols,...
            basedir,filename_diaries,...
            num_states,num_year_divisions,dt,ages,num_age_steps,start_year,T0,end_year,...
            theta,CFR_Acute,rate_6months,ECofactor,p_ChronicCarriage)

    end

    time_taken_for_sensitivity_analysis = datetime('now') - begin_time_sensitivity_analysis;
    disp(['The duration of this sensitivity analysis (' sensitivity_analysis ') was ' char(time_taken_for_sensitivity_analysis) ' hours.\n\n'])
end






diary off


