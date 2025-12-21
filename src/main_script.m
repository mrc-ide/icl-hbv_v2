%% Version of icl-hbv model modified for Cor-HepB projects by MP, December 2025.

clear

%% To remove to allow full run:
%% num_stochas_runs = 2;
%% num_scenarios = 1;
%% num_all_countries = 3;
%% num_countries = 3;
%% sensitivity_analysis_list = {'default'};


% This script creates 4 x 200 = 800 results files, each of which is less than 20 MB in size.
% A single sensitivity analysis, which contains 110 countries, each containing 200 particles, will take about 600 hours (25 days) to run on a desktop computer.
% This script, which contains 4 sensitivity analyses, will therefore take about 100 days to run on a desktop computer.


currentFolder = pwd;
% MP: I have modified this to put all MatLab code in src/ as this is simpler.
pat = fullfile('icl-hbv','src');
if ~endsWith(currentFolder,pat)
    warning(['Please run this script from within the folder ' pat])
end
% fileparts() gives a file path one directory level up.
basedir = fileparts(currentFolder); % path for folder one level up from script.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets up number of scenarios, sensitivity analyses and stochastic runs: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MP: Number of stochastic runs per country x scenario x sensitivity analysis:
% MP: There are 200 runs from calibration so this is default. Use less for testing.
%num_stochas_runs = 200;
num_stochas_runs = 2;


sensitivity_analysis_list = {'default'};
%sensitivity_analysis_list = {'default','infant_100','treat_medium','treat_high'};
num_sensitivity_analyses = length(sensitivity_analysis_list);


ListOfScenarios = {...
    'Status quo infant & BD',...
    'Status quo infant & BD expansion to 25%',...
    'Status quo infant & BD expansion to 50%',...
    'Status quo infant & BD expansion to 75%',...
    'Status quo infant & BD expansion to 90%',...
    'Status quo infant & BD drop 5 2020',...
    'Status quo infant & BD drop 10 2020',...
    'Status quo infant & BD drop 15 2020',...
    'Status quo infant & BD drop 20 2020',...
    'Status quo infant & BD delayed expansion 2023 to 2030',...
    'Status quo infant & BD delayed expansion 2023 to 2033',...
    'Status quo infant & BD delayed expansion 2025 to 2040'};
num_scenarios = length(ListOfScenarios);
assert(num_scenarios==12)
num_scenarios = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up parameters - first set parameters that are set by hand in code:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_year = 1890;
T0 = 2101-start_year; % must go up to year 2101 to get incidence readings up to year 2100
end_year = start_year + T0;
num_year_divisions = 10;
assert(rem(log10(num_year_divisions),1)==0) % ensure that num_year_divisions is a multiple of 10
% log10(1) = 0, log10(10) = 1, log10(100) = 2, etc.
% rem(x,1) gives 0 if x is an integer
dt = 1/num_year_divisions; % time step

ages = 0:dt:(100 - dt); % The age group at the beginning of each compartment; 1 x 1000 double; [0 0.1 0.2 ... 99.8 99.9]
num_age_steps = length(ages); % Number of age-categories (assuming max age is 100 years (i.e. oldest person we see at the start of an iteration is aged 100-dt)

theta = 0.01; % proportion of persons newly infected that develop severe symptoms during acute infection
CFR_Acute = 0.05; % proportion of people with severe acute infection that die from severe acute infection
rate_6months = 2; % the rate necessary to move everyone out of States 14 and 15 into State 2 in 6 months

ECofactor = 15; % Multiple for rate of transmission for horizontal transmission with HbEAg vs HbSAg
% epsilon_2 and epsilon_3 in Table S2, p. 6 of the appendix to the publication



% risk of chronic carriage due to horizontal and vertical transmission
% This is an "anonymous" function - similar to inline function in C (matlab
% has inline functions but recommends anonymous functions as more secure and faster).
risk_ages_edmunds_func = @(age) exp(-0.645*age^(0.455));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SERNIK - stuff to do: move these into here.
% params.dwvec = dwvec;
% params.SpeedUpELoss_F = 9.5;
% params.CancerRate_WomenCoFactor = 1;
% params.CirrhosisRate_WomenCoFactor = 1;    
% Progression Parameters
% Prog gives basic (identical across age, gender and treatment) progression rates between disease states
% Transactions then gives one finer control by allowing one to adjust rates
% according to age, gender and treatment
% Transactions.Values is then used in the disease progression part of the
% model

% Prog = zeros(num_states, num_states); % Non-Age Specific Prog parameters stored as (from, to)
% Fill-in transitions from Immune Tolerant
% Prog(2, 3) = 0.1;
% Fill-in transitions from Immune Reactive
% Prog(3, 4) = 0.05;
% Prog(3, 5) = 0.005;
% Prog(3, 6) = 0.028;
% Fill-in constant transitions from Asymptomatic Carrier
% Prog(4, 5) = 0.01;
% Prog(4, 9) = 0.01;
% Fill-in transitions from Chronic Hep B.
% Prog(5, 6) = 0.04;
% Fill-in transitions from Compensated Cirrhosis
% Prog(6, 7) = 0.04;
% Prog(6, 11) = 0.04;
% Fill-in transitions from Decompensated Cirrhsis
% Prog(7, 8) = 0.04;
% Prog(7, 11) = 0.30;
% Fill-in transitions from Liver Cancer
% Prog(8, 11) = params.CancerDeathRate;
% Fill-in transition from TDF-Treatment
% Prog(10, 11) = 0.001;
% Fill-in transition from 3TC-Treatment
% Prog(12, 13) = 0.2;
% Fill-in transition from Failed 3TC-Treatment
% Prog(13, 8) = 0.04;
% Prog(13, 11) = 0.3;
% Fill-in transition from Severe acute
% Prog(15, 11) = CFR_Acute * rate_6months;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up parameters - now load parameters set in calibration:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load(fullfile(basedir,'resources','ListOfISOs.mat')) % contains ListOfISOs
load(fullfile(basedir,'resources','vaccination_coverages.mat')) % contains BD_table and HepB3_table
load(fullfile(basedir,'resources','reference_prevalences.mat')) % contains country_s_e_HCCdeaths_map
load(fullfile(basedir,'resources','params_map.mat')) % contains params_map and dwvec
load(fullfile(basedir,'resources','stochastic_parameters_mat.mat')) % contains stochas_params_mat
load(fullfile(basedir,'resources','treatment_2016_map.mat')) % contains num_in_treatment_2016_map and pop_size_HBsAg_treatment_map
load(fullfile(basedir,'resources','treatment_rates_map.mat')) % contains treatment_rates_map



% treatment_boundaries_vec contains the following (NOTE - this is *not* a complete list, but only treatment_boundaries_vec([1 2 3 5]) are used in the code:
% (1) treatment year, 
% (2) rate to keep number of people in treatment constant,
% (3) rate to have 40% of eligible people in treatment by 2030, whether the constant rate is less than the 40% rate
% (5) rate to have 80% of eligible people in treatment by 2030, whether the 40% rate is less than the 80% rate



filename_diaries = 'diary_countries.out';
diary(fullfile(basedir,'outputs',filename_diaries))
% writes a copy of all subsequent keyboard input and the resulting output (except it does not include graphics) to the named file
% if the file already exists, output is appended to the end of the file




ListOfAllISOs = ListOfISOs;
num_all_countries = length(ListOfAllISOs);
assert(num_all_countries==110)
num_countries = length(ListOfISOs);
assert(num_countries==110)

num_all_countries = 3;
num_countries = 3;

% variables for importing values for the 7 fitted parameters for a particular country-particle combination
country_start_cols = 2:8:(1+num_all_countries*8); % there are 1 + num_all_countries*8 + 2 columns in stochas_params_mat


begin_time_run = datetime('now');
disp(['This analysis started at ' datestr(begin_time_run)])


disp(['Number of sensitivity analyses: ' num2str(num_sensitivity_analyses)])
disp(['Number of scenarios: ' num2str(num_scenarios)])
disp(['Number of countries: ' num2str(num_countries)])
disp(['Number of runs per country: ' num2str(num_stochas_runs)])




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


% does not differ over disease stages
% 4 dimensions necessary so that it can be multiplied by X
risk_p_ChronicCarriage = arrayfun(@(xx) risk_ages_edmunds_func(xx), ages(:,ages>0.5));
risk_p_ChronicCarriage = repmat(risk_p_ChronicCarriage,1,1,2,2);

% Set up the probability that an individual becomes a chronic carrier given
% they get infected at a specific age, their sex, and *accessible* (whether
% can be reached by treatment progs, 1=no, 2=yes).
p_ChronicCarriage = 0.885 * ones(1, num_age_steps, 2, 2);
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
    sensitivity_analysis_hours_vec(sensitivity_analysis_num) = time_taken_for_sensitivity_analysis;
    assert(all(sensitivity_analysis_hours_vec(1:sensitivity_analysis_num)>0))
    average_time_per_sensitivity_analysis = mean(sensitivity_analysis_hours_vec(1:sensitivity_analysis_num));
    min_time_per_sensitivity_analysis = min(sensitivity_analysis_hours_vec(1:sensitivity_analysis_num));
    max_time_per_sensitivity_analysis = max(sensitivity_analysis_hours_vec(1:sensitivity_analysis_num));
    num_sensitivity_analyses_left = num_sensitivity_analyses - sensitivity_analysis_num;
    mean_time_left = num_sensitivity_analyses_left * average_time_per_sensitivity_analysis;
    min_time_left = num_sensitivity_analyses_left * min_time_per_sensitivity_analysis;
    max_time_left = num_sensitivity_analyses_left * max_time_per_sensitivity_analysis;
    if num_sensitivity_analyses_left>0
        disp(['There are ' num2str(num_sensitivity_analyses_left) ' sensitivity analyses left, which will take about ' char(mean_time_left) ' hours.'])
    end


end


end_time_run = datetime('now');
disp(end_time_run)
time_taken_for_run = end_time_run - begin_time_run;
disp(['The duration of this analysis was ' char(time_taken_for_run) ' hours.\n\n'])


diary off


