function output = HBVmodel(source_HBsAg,...
    num_disease_states,num_year_divisions,dt,ages,num_age_steps,start_year,num_years_simul,...
    theta,ECofactor, ...
    treatment_rate_params, treat_start_year,treat_coverage_in_2016, ...
    params, PAP_VL_params, PAP_cov_params, ...
    p_ChronicCarriage,Prog,Transactions, ...
    scenario_BDcoverage, scenario_BDcoverage_fromMAP, ...
    scenario_BDcoverage_fromCPAD, scenario_HepB3coverage, ...
    HAS_TREATMENT,...
    ISO, scenario_num, scenario_CohortTesting, ...
    num_year_1980_2100, life_expectancy, ...
    stochas_run_str, sensitivity_analysis, basedir, store_results_as_text)


% X-stocks are (infection_state, age, sex(1=women, 2=men), accessible*)   {*accessible
% specifies whether this person can be reached by treatment progs, 1=no, 2=yes}
% So X() is 4D array of size num_disease_states *  num_age_steps * num_sexes * num_treat_blocks
num_sexes = 2; % F=1, M=2.
i_female = 1;
i_male = 2;
num_treat_blocks = 2;
i_notreat = 1;
%% i_yestreat = 2; %% not currently needed.

%% Note that these are also defined in country_level_analyses.m (need to match).
I_NO_COHORT_TEST = 1;
I_COHORT_TEST = 2;

i_Susc = 1;         % 'Susceptible', 
i_ImmTol = 2;       % 'HBV: Immune Tolerant' : HBeAg+ with very high HBV DNA (>1e6IU/ml), normal ALT
i_ImmReact = 3;     % 'HBV: Immune Reactive' :  HBeAg+ with high HBV DNA (>20000IU/ml), elevated ALT
i_AsymptCarr = 4;   % 'HBV: Asymptomatic Carrier': HBeAg- low/undetectable HBV DNA, normal ALT
i_Chronic = 5;      % 'HBV: Chronic Hep B': HBeAg-, moderate to high HBV DNA levels, fluctuating/persistently elevated ALT 
i_CompCirr = 6;     % 'HBV: Comp Cirrhosis',
i_DecompCirr = 7;   % 'HBV: Decomp Cirrhosis',
i_HCC = 8;          % 'HBV: Liver Cancer',
i_Immune = 9;       % 'HBV: Immune (Rec. or vacc.)',
i_TDFtreat = 10;    % 'HBV: TDF-Treatment',
i_HBVdeath = 11;    % 'Prematurely dead due to HBV', ... % 11
i_3TCtreat = 12;    % '3TC-Treatment', ... % 12
i_3TCfailed = 13;   % 'Failed 3TC-Treatment', ...  % 13
i_NonSevAcute = 14; % 'Non-severe acute', ...  % 14
i_SevereAcute = 15; % 'Severe acute' ...  % 15
 
i_alive = [1:10 12:15];
i_eAgpos_chronic = 2:3;  %% Immune Tolerant, Immune Reactive, Non-severe + severe acute.
i_eAgpos = [2:3 14:15];  %% Immune Tolerant, Immune Reactive, Non-severe + severe acute.
i_sAgpos_notEagpos_notreat = [4:8 13];     %% Asymptomatic carrier, Chronic, Comp+Decom Cirr, HCC, failed 3TC.
i_treateligible = [3 5 6 7]; %% Immune Reactive, Chronic, Comp+Decomp Cirr
i_sAgpos = [2:8 10 12:15];   %% Includes 10 (TDFtreat) and 12 (3TCtreat) states
i_sAgpos_chronic = [2:8 10 12:13]; 

DUMMY_VALUE = -99;  % Used in initialising arrays to a dummy value (-99 should be easy to spot).

%% Establish basic simulation parameters
agegroups_5yr = 1 + floor(ages / 5); % categorises the ages into age-groups of 5 year width; 1 x 1000 double; [1 1 ... 20 20], each number present 50 times
agegroups_1yr = 1 + floor(ages); % categorises the ages into age-groups of 1 year width; 1 x 1000 double; [1 1 ... 100 100], each number present 10 times
					
i6mo = find(ages >= 0.5, 1); % markers for key age boundaries
i1y = find(ages >= 1, 1);
i5y = find(ages >= 5, 1);
i15y = find(ages >= 15, 1);

end_year = start_year + num_years_simul; % 2101
TimeSteps = start_year:dt:end_year; % 1 x 2101 double; [1890 1890.1 1890.2 ... 2099.8 2099.9 2100 2100.1 ... 2100.8 2100.9 2101]


%% Establish intervention parameters
% Definition of natural history model (number of stages)
% Definition of stocks (first dimension is age, second dimension is status)


% ---- Intervention Parameters ----

%% MP: note - I am using p_VerticalTransmission_HbSAg_NoBD instead of "p_VerticalTransmission_HbSAg_NoIntv" (similarly for EAg).
%% The "Intv" refers to birth-dose vaccination (either normal vaccination, microarray patches (MAPs), or compact prefilled auto-disable devices (CPAD)).
%% WLASNY - cut 11 April 2026 (as we move transmission to be by high/low VL):
%% p_VerticalTransmission_HbSAg_NoBD = params.p_VerticalTransmission_HbSAg_NoIntv; % probability of transmission from an HBeAg-, HBsAg+ mother to her baby without intervention
%% p_VerticalTransmission_HbEAg_NoBD = params.p_VerticalTransmission_HbEAg_NoIntv;

%% p_VerticalTransmission_HbSAg_BD = p_VerticalTransmission_HbSAg_NoBD * (1 - params.Efficacy_BirthDoseVacc_HbSAg);
%% p_VerticalTransmission_HbSAg_BirthDose_MAP_CPAD = p_VerticalTransmission_HbSAg_NoBD * (1 - efficacy_MAP_CPAD_HbSAg);
%% assert(p_VerticalTransmission_HbSAg_NoBD>=0 && p_VerticalTransmission_HbSAg_NoBD<=1)
%% mustBeBetween(p_VerticalTransmission_HbSAg_BD, 0, 1)
%% mustBeBetween(p_VerticalTransmission_HbSAg_BirthDose_MAP_CPAD, 0, 1)
% probability of transmission from an HBeAg+ mother to her baby after the baby is given BD vaccination
%%% KAWA:
%%%p_VerticalTransmission_HbEAg_BD = p_VerticalTransmission_HbEAg_NoBD * (1 - params.Efficacy_BirthDoseVacc_HbEAg); 
%%%p_VerticalTransmission_HbEAg_BirthDose_MAP_CPAD = p_VerticalTransmission_HbEAg_NoBD * (1 - efficacy_MAP_CPAD_HbEAg); 
%%%p_VerticalTransmission_Tr_NoBD = p_VerticalTransmission_HbEAg_NoBD * (1 - Efficacy_Treatment_MTCT); % probability of transmission from an HBeAg+ mother on treatment to her baby without intervention


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TAM: PAP chunk 1:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PAP in addition to Birth dose (Fraction of those who have BD that also get PAP)
% cov_BirthDoseAndTDF_EAgHighVL = 0;
% cov_BirthDoseAndTDF_SAgHighVL = 0;
% cov_BirthDoseAndTDF_EAgLowVL = 0;
% cov_BirthDoseAndTDF_SAgLowVL = 0;
% % PAP instead of BD (Fraction of those who do not get BD that do get PAP)
% cov_TDFOnly_EAgHighVL = 0;
% cov_TDFOnly_SAgHighVL = 0;
% cov_TDFOnly_EAgLowVL = 0;
% cov_TDFOnly_SAgLowVL = 0;
% % This scale-up parameter pertains to both types of PAP usage (instantaneously)
% TScaleup_PAP = 2025;

% PMTCT Asssumptions
% IMPORT THE PAP_VL_params (Now mandatory to give a full set of parameters)
%%FracEPosHighVL = PAP_VL_params.FracEPosHighVL;
%%FracSPosHighVL = PAP_VL_params.FracSPosHighVL;
%%p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL = PAP_VL_params.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL;
%%p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL = PAP_VL_params.PAP_VL_params.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL;
%%pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc = PAP_VL_params.pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc;
%%pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP = PAP_VL_params.pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP;
%%pr_VerticalTransmission_HbSAgLowVL_PAP = PAP_VL_params.pr_VerticalTransmission_HbSAgLowVL_PAP;
%%pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc = PAP_VL_params.pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc;
%%pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP = PAP_VL_params.pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP;
%%pr_VerticalTransmission_HbSAgHighVL_PAP = PAP_VL_params.pr_VerticalTransmission_HbSAgHighVL_PAP;
%%pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc = PAP_VL_params.pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc;
%%pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP = PAP_VL_params.pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP;
%%pr_VerticalTransmission_HbEAgLowVL_PAP = PAP_VL_params.pr_VerticalTransmission_HbEAgLowVL_PAP;
%%pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc = PAP_VL_params.pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc;
%%pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP = PAP_VL_params.pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP;
%%pr_VerticalTransmission_HbEAgHighVL_PAP = PAP_VL_params.pr_VerticalTransmission_HbEAgHighVL_PAP;


% % Confirm that the correct PAP_VL_params have been entered:
% PAP_VL_params_required_params = {'FracEPosHighVL',...
%                             'FracSPosHighVL',...
%                             'p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL',...
%                             'p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL',...
%                             'pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc',...
%                             'pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP',...
%                             'pr_VerticalTransmission_HbSAgLowVL_PAP',...
%                             'pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc',...
%                             'pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP',...
%                             'pr_VerticalTransmission_HbSAgHighVL_PAP',...
%                             'pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc',...
%                             'pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP',...
%                             'pr_VerticalTransmission_HbEAgLowVL_PAP',...
%                             'pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc',...
%                             'pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP',...
%                             'pr_VerticalTransmission_HbEAgHighVL_PAP'};
% 
% for i=1:length(PAP_VL_params_required_params)
%     assert(exist(PAP_VL_params_required_params{i},'var')>0)
% end

% Compute p_VertTrans_HbSAgHighVL_NoIntv and p_VertTrans_HbSAgLowVL_NoIntv
p_HbSAg_av = params.p_VerticalTransmission_HbSAg_NoIntv; % From fitting (this average rate to be preserved)
%% These are now calculated in assign_PAP_VL_params:
p_VertTrans_HbSAgLowVL_NoIntv = PAP_VL_params.p_VertTrans_HbSAgLowVL_NoIntv;
p_VertTrans_HbSAgHighVL_NoIntv = PAP_VL_params.p_VertTrans_HbSAgHighVL_NoIntv;

% Check close to original p_HbSAg_av value:
assert(abs(p_HbSAg_av - (p_VertTrans_HbSAgLowVL_NoIntv*(1-PAP_VL_params.FracSPosHighVL) + p_VertTrans_HbSAgHighVL_NoIntv*PAP_VL_params.FracSPosHighVL))<0.001)
%mustBeBetween(p_VertTrans_HbSAgLowVL_NoIntv, 0, 1)
%mustBeBetween(p_VertTrans_HbSAgHighVL_NoIntv, 0, 1)
assert(p_VertTrans_HbSAgLowVL_NoIntv>=0 && p_VertTrans_HbSAgLowVL_NoIntv<=1)
assert(p_VertTrans_HbSAgHighVL_NoIntv>=0 && p_VertTrans_HbSAgHighVL_NoIntv<=1)


p_HbEAg_av = params.p_VerticalTransmission_HbEAg_NoIntv; % From assumption in the version of the model used in fitting  (this average rate to be preserved) [NB. if this varied in different region then specifiy that here!)]
%% These are now calculated in assign_PAP_VL_params:
p_VertTrans_HbEAgLowVL_NoIntv = PAP_VL_params.p_VertTrans_HbEAgLowVL_NoIntv;
p_VertTrans_HbEAgHighVL_NoIntv = PAP_VL_params.p_VertTrans_HbEAgHighVL_NoIntv;

% Check close to original p_HbEAg_av value:
assert(abs(p_HbEAg_av - (p_VertTrans_HbEAgLowVL_NoIntv*(1-PAP_VL_params.FracEPosHighVL) + p_VertTrans_HbEAgHighVL_NoIntv*PAP_VL_params.FracEPosHighVL))<0.001)
%mustBeBetween(p_VertTrans_HbEAgLowVL_NoIntv, 0, 1)
%mustBeBetween(p_VertTrans_HbEAgHighVL_NoIntv, 0, 1)
assert(p_VertTrans_HbEAgLowVL_NoIntv>=0 && p_VertTrans_HbEAgLowVL_NoIntv<=1)
assert(p_VertTrans_HbEAgHighVL_NoIntv>=0 && p_VertTrans_HbEAgHighVL_NoIntv<=1)

% Compute the transmission probabilities using probability ratios (i.e.
% (1-effectiveness)). Note that for now we do not assume independence
% between the different interventions.
p_VertTrans_HbSAgLowVL_BD         = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_BD * p_VertTrans_HbSAgLowVL_NoIntv;
p_VertTrans_HbSAgLowVL_PAP        = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_PAP * p_VertTrans_HbSAgLowVL_NoIntv;
p_VertTrans_HbSAgLowVL_MAP        = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_MAP * p_VertTrans_HbSAgLowVL_NoIntv;
p_VertTrans_HbSAgLowVL_CPAD       = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_CPAD * p_VertTrans_HbSAgLowVL_NoIntv;
p_VertTrans_HbSAgLowVL_BD_PAP     = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_BD_PAP * p_VertTrans_HbSAgLowVL_NoIntv;
p_VertTrans_HbSAgLowVL_MAP_PAP    = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_MAP_PAP * p_VertTrans_HbSAgLowVL_NoIntv;
p_VertTrans_HbSAgLowVL_CPAD_PAP   = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_CPAD_PAP * p_VertTrans_HbSAgLowVL_NoIntv;
%%p_VertTrans_HbSAgLowVL_Treat      = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_Treat * p_VertTrans_HbSAgLowVL_NoIntv;
%%p_VertTrans_HbSAgLowVL_BD_Treat   = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_BD_Treat * p_VertTrans_HbSAgLowVL_NoIntv;
%%p_VertTrans_HbSAgLowVL_MAP_Treat  = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_MAP_Treat * p_VertTrans_HbSAgLowVL_NoIntv;
%%p_VertTrans_HbSAgLowVL_CPAD_Treat = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_CPAD_Treat * p_VertTrans_HbSAgLowVL_NoIntv;

p_VertTrans_HbSAgHighVL_BD         = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_BD * p_VertTrans_HbSAgHighVL_NoIntv;
p_VertTrans_HbSAgHighVL_PAP        = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_PAP * p_VertTrans_HbSAgHighVL_NoIntv;
p_VertTrans_HbSAgHighVL_MAP        = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_MAP * p_VertTrans_HbSAgHighVL_NoIntv;
p_VertTrans_HbSAgHighVL_CPAD       = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_CPAD * p_VertTrans_HbSAgHighVL_NoIntv;
p_VertTrans_HbSAgHighVL_BD_PAP     = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_BD_PAP * p_VertTrans_HbSAgHighVL_NoIntv;
p_VertTrans_HbSAgHighVL_MAP_PAP    = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_MAP_PAP * p_VertTrans_HbSAgHighVL_NoIntv;
p_VertTrans_HbSAgHighVL_CPAD_PAP   = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_CPAD_PAP * p_VertTrans_HbSAgHighVL_NoIntv;
%%p_VertTrans_HbSAgHighVL_Treat      = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_Treat * p_VertTrans_HbSAgHighVL_NoIntv;
%%p_VertTrans_HbSAgHighVL_BD_Treat   = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_BD_Treat * p_VertTrans_HbSAgHighVL_NoIntv;
%%p_VertTrans_HbSAgHighVL_MAP_Treat  = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_MAP_Treat * p_VertTrans_HbSAgHighVL_NoIntv;
%%p_VertTrans_HbSAgHighVL_CPAD_Treat = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_CPAD_Treat * p_VertTrans_HbSAgHighVL_NoIntv;

p_VertTrans_HbEAgLowVL_BD         = PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_BD * p_VertTrans_HbEAgLowVL_NoIntv;
p_VertTrans_HbEAgLowVL_PAP        = PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_PAP * p_VertTrans_HbEAgLowVL_NoIntv;
p_VertTrans_HbEAgLowVL_MAP        = PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_MAP * p_VertTrans_HbEAgLowVL_NoIntv;
p_VertTrans_HbEAgLowVL_CPAD       = PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_CPAD * p_VertTrans_HbEAgLowVL_NoIntv;
p_VertTrans_HbEAgLowVL_BD_PAP     = PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_BD_PAP * p_VertTrans_HbEAgLowVL_NoIntv;
p_VertTrans_HbEAgLowVL_MAP_PAP    = PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_MAP_PAP * p_VertTrans_HbEAgLowVL_NoIntv;
p_VertTrans_HbEAgLowVL_CPAD_PAP   = PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_CPAD_PAP * p_VertTrans_HbEAgLowVL_NoIntv;
%%p_VertTrans_HbEAgLowVL_Treat      = PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_Treat * p_VertTrans_HbEAgLowVL_NoIntv;
%%p_VertTrans_HbEAgLowVL_BD_Treat   = PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_BD_Treat * p_VertTrans_HbEAgLowVL_NoIntv;
%%p_VertTrans_HbEAgLowVL_MAP_Treat  = PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_MAP_Treat * p_VertTrans_HbEAgLowVL_NoIntv;
%%p_VertTrans_HbEAgLowVL_CPAD_Treat = PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_CPAD_Treat * p_VertTrans_HbEAgLowVL_NoIntv;


p_VertTrans_HbEAgHighVL_BD         = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_BD * p_VertTrans_HbEAgHighVL_NoIntv;
p_VertTrans_HbEAgHighVL_PAP        = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_PAP * p_VertTrans_HbEAgHighVL_NoIntv;
%%p_VertTrans_HbEAgHighVL_Treat      = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_Treat * p_VertTrans_HbEAgHighVL_NoIntv;
p_VertTrans_HbEAgHighVL_MAP        = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_MAP * p_VertTrans_HbEAgHighVL_NoIntv;
p_VertTrans_HbEAgHighVL_CPAD       = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_CPAD * p_VertTrans_HbEAgHighVL_NoIntv;
p_VertTrans_HbEAgHighVL_BD_PAP     = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_BD_PAP * p_VertTrans_HbEAgHighVL_NoIntv;
%%p_VertTrans_HbEAgHighVL_BD_Treat   = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_BD_Treat * p_VertTrans_HbEAgHighVL_NoIntv;
p_VertTrans_HbEAgHighVL_MAP_PAP    = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_MAP_PAP * p_VertTrans_HbEAgHighVL_NoIntv;
%%p_VertTrans_HbEAgHighVL_MAP_Treat  = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_MAP_Treat * p_VertTrans_HbEAgHighVL_NoIntv;
p_VertTrans_HbEAgHighVL_CPAD_PAP   = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_CPAD_PAP * p_VertTrans_HbEAgHighVL_NoIntv;
%%p_VertTrans_HbEAgHighVL_CPAD_Treat = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_CPAD_Treat * p_VertTrans_HbEAgHighVL_NoIntv;

%% Treatment:
%% By default, treatment is only for EAg+ individuals:
p_VertTrans_HbEAg_Treat      = PAP_VL_params.pRatio_VertTrans_Treat * params.p_VerticalTransmission_HbEAg_NoIntv;
p_VertTrans_HbEAg_Treat_BD   = PAP_VL_params.pRatio_VertTrans_Treat_BD * params.p_VerticalTransmission_HbEAg_NoIntv;
p_VertTrans_HbEAg_Treat_MAP  = PAP_VL_params.pRatio_VertTrans_Treat_MAP * params.p_VerticalTransmission_HbEAg_NoIntv;
p_VertTrans_HbEAg_Treat_CPAD = PAP_VL_params.pRatio_VertTrans_Treat_CPAD * params.p_VerticalTransmission_HbEAg_NoIntv;

%% Allow treatment for EAg- SAg+ individuals if needed (not currently used):
%%p_VertTrans_HbSAg_Treat      = PAP_VL_params.pRatio_VertTrans_Treat * params.p_VerticalTransmission_HbSAg_NoIntv;
%%p_VertTrans_HbSAg_Treat_BD   = PAP_VL_params.pRatio_VertTrans_Treat_BD * params.p_VerticalTransmission_HbSAg_NoIntv;
%%p_VertTrans_HbSAg_Treat_MAP  = PAP_VL_params.pRatio_VertTrans_Treat_MAP * params.p_VerticalTransmission_HbSAg_NoIntv;
%%p_VertTrans_HbSAg_Treat_CPAD = PAP_VL_params.pRatio_VertTrans_Treat_CPAD * params.p_VerticalTransmission_HbSAg_NoIntv;


% Check that all the the PAP VL parameters lie in expected ranges:
p_ratioSAg_HVL_LVL_noint = p_VertTrans_HbSAgHighVL_NoIntv/p_VertTrans_HbSAgLowVL_NoIntv;
p_ratioEAg_HVL_LVL_noint = p_VertTrans_HbEAgHighVL_NoIntv/p_VertTrans_HbEAgLowVL_NoIntv;
assert(p_ratioSAg_HVL_LVL_noint>=1.0)
assert(p_ratioEAg_HVL_LVL_noint>=1.0)
Snames_PAP_VL_params = fieldnames(PAP_VL_params);
% Many of these are probability ratios, but we expect them to be <=1 (because they relate to interventions that reduce transmission):
%%for i = [1, 2, 5:numel(Snames_PAP_VL_params)]
for i = 1:numel(Snames_PAP_VL_params)
    %mustBeBetween(PAP_VL_params.(Snames_PAP_VL_params{i}),0,1)
    assert(PAP_VL_params.(Snames_PAP_VL_params{i})>=0 && PAP_VL_params.(Snames_PAP_VL_params{i})<=1)
end
    
%   Within each e/vl cat, confirm appropriate ordering
assert(safe_greater_or_equal_to(p_VertTrans_HbSAgLowVL_NoIntv,  p_VertTrans_HbSAgLowVL_BD))
assert(safe_greater_or_equal_to(p_VertTrans_HbSAgLowVL_NoIntv,  p_VertTrans_HbSAgLowVL_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbSAgLowVL_BD,      p_VertTrans_HbSAgLowVL_BD_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbSAgHighVL_NoIntv, p_VertTrans_HbSAgHighVL_BD))
assert(safe_greater_or_equal_to(p_VertTrans_HbSAgHighVL_NoIntv, p_VertTrans_HbSAgHighVL_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbSAgHighVL_BD,     p_VertTrans_HbSAgHighVL_BD_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgLowVL_NoIntv,  p_VertTrans_HbEAgLowVL_BD))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgLowVL_NoIntv,  p_VertTrans_HbEAgLowVL_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgLowVL_BD,      p_VertTrans_HbEAgLowVL_BD_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgHighVL_NoIntv, p_VertTrans_HbEAgHighVL_BD))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgHighVL_NoIntv, p_VertTrans_HbEAgHighVL_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgHighVL_BD,     p_VertTrans_HbEAgHighVL_BD_PAP))
    
%   That the high VL cat is always more transmissive than the low VL cat
assert(safe_greater_or_equal_to(p_VertTrans_HbSAgHighVL_NoIntv, p_VertTrans_HbSAgLowVL_NoIntv))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgHighVL_NoIntv, p_VertTrans_HbEAgLowVL_NoIntv))
assert(safe_greater_or_equal_to(p_VertTrans_HbSAgHighVL_BD,     p_VertTrans_HbSAgLowVL_BD))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgHighVL_BD,     p_VertTrans_HbEAgLowVL_BD))
assert(safe_greater_or_equal_to(p_VertTrans_HbSAgHighVL_BD_PAP, p_VertTrans_HbSAgLowVL_BD_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgHighVL_BD_PAP, p_VertTrans_HbEAgLowVL_BD_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbSAgHighVL_PAP,    p_VertTrans_HbSAgLowVL_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgHighVL_PAP,    p_VertTrans_HbEAgLowVL_PAP))

%   That the 'e cat' is always the same or more transmissive the the s cat
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgHighVL_NoIntv, p_VertTrans_HbSAgHighVL_NoIntv))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgLowVL_NoIntv,  p_VertTrans_HbSAgLowVL_NoIntv))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgHighVL_BD,     p_VertTrans_HbSAgHighVL_BD))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgLowVL_BD,      p_VertTrans_HbSAgLowVL_BD))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgHighVL_BD_PAP, p_VertTrans_HbSAgHighVL_BD_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgLowVL_BD_PAP,  p_VertTrans_HbSAgLowVL_BD_PAP)  )  
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgHighVL_PAP,    p_VertTrans_HbSAgHighVL_PAP))
assert(safe_greater_or_equal_to(p_VertTrans_HbEAgLowVL_PAP,     p_VertTrans_HbSAgLowVL_PAP))
%% End of checks.
% ------------------------------------------------------------------------



    %% MP: PAP coverage is now done on a scenario-by-scenario basis in country_level_analyses.m
    % %% If no PAP, then set these to zero:
    %% TAM delete me.
    %% Now replaced with e.g. PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL
    % % Coverage of PAP among those with BD
    % [cov_BirthDoseAndTDF_EAgHighVL_itt, ...
    %     cov_BirthDoseAndTDF_EAgLowVL_itt, ...
    %     cov_BirthDoseAndTDF_SAgHighVL_itt, ...
    %     cov_BirthDoseAndTDF_SAgLowVL_itt] = ...
    %     deal(zeros(size(TimeSteps)));

    %% TAM - delete me.
    % % Coverage of PAP among those not with BD
    %% Now e.g. PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgHighVL
    % [cov_TDFOnly_EAgHighVL_itt, ...
    %     cov_TDFOnly_EAgLowVL_itt, ...
    %     cov_TDFOnly_SAgHighVL_itt, ...
    %     cov_TDFOnly_SAgLowVL_itt] = ...
    %     deal(zeros(size(TimeSteps)));


%% MP: This is a little bit spaghetti code.
%% MP: Note - this deal(0) line is actually necessary. 
%% This code chunk initialises these variables to be zero.
%% These are then stored, and then finally updated for the next step.
%% This isn't ideal - not sure if there is an underlying logic making this order necessary 
%% but it means we are storing the previous timestep's value at each timestep.
[births_toHbEAgWomenHighVL, births_toHbEAgWomenLowVL, births_toHbSAgWomenHighVL, ...
    births_toHbSAgWomenLowVL, births_Total, ...
    babiesChronic_from_HbEAgWomenHighVL, babiesChronic_from_HbEAgWomenLowVL, ...
    babiesChronic_from_HbSAgWomenHighVL, babiesChronic_from_HbSAgWomenLowVL,...
    ratebirthdoses, ratebirthdoses_MAP, ratebirthdoses_CPAD,...
    pregnantWomenNeedToScreen,...
    num_mothers_PAP_HbEAg_HighVL, num_mothers_PAP_HbEAg_LowVL, num_mothers_PAP_HbSAg_HighVL, num_mothers_PAP_HbSAg_LowVL,...
    RateOfPAPInitiation, HBVPositivePregnantWomenAtANC] = deal(0);

  %% For now let's just count the number of people starting treatment (unstratified by age/sex):
  number_starting_treatment_to_print = 0;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TAM: End of PAP chunk 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ---- Load Epidemiological Parameters from fitting procedure ----

% beta_AGE are the (per-year) rate of horizontal transmission from an infected person (age gp AGE) who is sAg+ but 
% eAG- (for eAg+, this is multiplied by ECofactor, and capped at 1/yr), and a susceptible person
% (also in age gp AGE). 
% Note that for some age groups they will appear twice (e.g. 1-4 are both in the "U5" and "1to15" groups).
% This corresponds to them having 2 separate types of interactions (with other U5s, and with 1to15s).
% Finally, the code allows beta to be made time-dependent (so declines in future) via beta_scaler below - this is 
% currently not implemented.
beta_U5 = params.beta_U5; % Rate of horizontal transmission between susceptible and infected persons - UNDER FIVE
beta_1to15 = params.beta_1to15;                                % Rate of generation transmission between susceptible and infected persons - All Ages
beta_5plus = params.beta_5plus;




% ----- Infection-relate parameters -----

%% The following would make the horiontal transmission probability time-dependent (decreasing by some fraction beta_scaler.
ReducInTransmission = params.ReducInTransmission;              % Fractional reduction in transmission. Currently set to 0
YearReducInTransmission = params.YearReducInTransmission;      % Turning point year for reduction. Currently set to 2100
DurReducInTransmission = 15;                                  % Time taken to complete change - MP: note that this is not actually that. 
%% MP: beta_scaler seems to be legacy code. Currently ReducInTransmission is 0. Otherwise (even with YearReducInTransmission=2100)
%% we still get some reduction in beta, and quite a large reduction after 2080 (reaching 50% reduction in 2100).
beta_scaler = ReducInTransmission ./ (1 + exp( (TimeSteps - (YearReducInTransmission)) ./ (DurReducInTransmission / 10) ));

beta_U5_SAg = beta_U5 * (1 - ReducInTransmission) + beta_U5 * zeros(size(beta_scaler));  % NOT TIME DEPENDENT as no beta_scaler.
beta_U5_EAg = min(1.0, beta_U5_SAg * ECofactor);

beta_1to15_SAg = beta_1to15 * (1 - ReducInTransmission) + beta_1to15 * beta_scaler;  % TIME DEPENDENT as beta_scaler is.
beta_1to15_EAg = min(1.0, beta_1to15_SAg * ECofactor);

beta_5plus_SAg = beta_5plus * (1 - ReducInTransmission) + beta_5plus * beta_scaler;  % TIME DEPENDENT as beta_scaler is.
beta_5plus_EAg = min(1.0, beta_5plus_SAg * ECofactor);



% X-stocks are (infec, age, sex(1=women, 2=men), accessible*)   {*accessible
% specifies whether this person can be reached by treatment progs, 1=no, 2=yes}



% Demography
StartPop = params.Pop_byAgeGroups_1950(agegroups_1yr, :) * dt;
% params.Pop_byAgeGroups_1950 is a 100 x 2 double of 100 age groups and 2 genders
% agegroups_1yr is a 1 x 1000 double; [1 1 ... 100 100], each number present (1/dt) times
% start population size like 1950 population
% expanding params.Pop_byAgeGroups_1950 from 1 year age steps to 0.1 year age steps
% each age group repeated (10=1/dt) times therefore multiply each entry by
% dt.

X = zeros(num_disease_states, num_age_steps, num_sexes, num_treat_blocks);
% dimensions: disease states, age, gender, accessible to treatment


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise prevalence using HBsAg data:
if strcmp(source_HBsAg,'Cui')
    StartPrev_byAgeGroups = params.HBsAg_prevs_middle_year_1;
    %% MP: Magic number 18
    assert(isequal(size(StartPrev_byAgeGroups),[18 num_sexes]))
    if any(isnan(StartPrev_byAgeGroups))
       nan_positions = isnan(StartPrev_byAgeGroups);
       %% MP: Magic number 1 (also line below)
       first_non_nan_pos = find(~nan_positions(:,1), 1);
       StartPrev_byAgeGroups(1:first_non_nan_pos,:) = repmat(StartPrev_byAgeGroups(first_non_nan_pos,:),first_non_nan_pos,1);
       last_non_nan_pos = find(~nan_positions(:,1), 1, 'last'); % find the last non NaN position
       StartPrev_byAgeGroups(last_non_nan_pos:end,:) = repmat(StartPrev_byAgeGroups(last_non_nan_pos,:),size(nan_positions,1)-last_non_nan_pos+1,1);
    end
    %% MP: Magic numbers 3 1
    StartPrev_byAgeGroups = [StartPrev_byAgeGroups(1:end-1,:); repmat(StartPrev_byAgeGroups(end,:),3,1)];
    % params.HBsAg_prevs_middle_year_1 is a 18 x 2 double of age group by gender
    % params.HBsAg_prevs_middle_year_1 age groups: 0--4 5--9 10--14 15--19 20--24 25--29 30--34 35--39 40--44 45--49 50--54 55--59 60--64 65--69 70--74 75--79 80--84 85+
    assert(isequal(size(StartPrev_byAgeGroups),[20 num_sexes]))

    NumSAg = StartPrev_byAgeGroups(agegroups_5yr, :) .* StartPop;
    % agegroups_5yr is a 1 x 1000 double; [1 1 ... 20 20], each number present 50 times
    % expanding StartPrev_byAgeGroups from 5 year age steps to 0.1 year age steps
    NumNotSAg = (1 - StartPrev_byAgeGroups(agegroups_5yr, :)) .* StartPop;
elseif strcmp(source_HBsAg,'CDA')
    %% MP: Magic numbers 99.9, 6, 1, 2
    StartPrev_byAgeGroups = [repmat(params.country_HBsAg_prevalences_by_ages_mid_1_young_old(1),num_year_divisions*(5.9-0.0)+1,2); ...
        repmat(params.country_HBsAg_prevalences_by_ages_mid_1_young_old(2),num_year_divisions*(99.9-6.0)+1,2)];
    % apply prevalence in 5-year-olds to 0 to 6 year olds; apply prevalence in all ages to 6 to 99 year olds
    assert(isequal(size(StartPrev_byAgeGroups),[num_age_steps num_sexes]))

    NumSAg = StartPrev_byAgeGroups .* StartPop;
    NumNotSAg = (1 - StartPrev_byAgeGroups) .* StartPop;
elseif strcmp(source_HBsAg,'WHO')
    under_5_pos_vec_len = length(find(ages<=5.0));
    over_5_pos_vec_len = length(find(ages>5.0));
    assert(under_5_pos_vec_len+over_5_pos_vec_len==num_age_steps)
    %% MP: Magic numbers 1, 2, 2, 2
    StartPrev_byAgeGroups = [ ...
        repmat(params.country_HBsAg_prevalences_by_ages_prevacc_young_old(1),under_5_pos_vec_len,2); ...
        repmat(params.country_HBsAg_prevalences_by_ages_prevacc_young_old(2),over_5_pos_vec_len,2) ...
        ];
    assert(isequal(size(StartPrev_byAgeGroups),[num_age_steps num_sexes]))

    NumSAg = StartPrev_byAgeGroups .* StartPop;
    NumNotSAg = (1 - StartPrev_byAgeGroups) .* StartPop;
end    

X(i_Susc, :, :, i_notreat) = NumNotSAg;
%% MP: Magic number 0.5 (and 1-0.5) - put in main_script.m
X(i_ImmReact, :, :, i_notreat) = 0.5 * NumSAg;
X(i_AsymptCarr, :, :, i_notreat) = 0.5 * NumSAg;

% Demography
% Prepare an index that will allow quick population of the mu vector from
% the demographic data input (uneven age-groupings)

%% MP: Magic numbers 2:21, 5. There are 21 age groups (0-0, 1-4, 5-9, 10-14,... 95-99). The "2" is because we first pretend
%% the 0-0 and 1-4 age groups are a single age group (index 2 as it will correspond to 1-4). We later set age gp 0-0 by hand.
%% The 5 is so that overall we cover the 1000 timesteps (dt=0.1) from 0-99.9 (num_year_divisions=1/dt; when pretending the 0-0 
%% and 1-4 age groups are a single group, we have 20 of these groups, so need a multiplier of 5=1000/(10*20). 
MappingFromDataToParam = repmat(2:21,5*num_year_divisions,1);
MappingFromDataToParam = MappingFromDataToParam(:);
MappingFromDataToParam(1:num_year_divisions) = 1; %% MP: now set age group 0-0 by hand.
% MappingFromDataToParam gives the value in the mortality vectors (21 values
% corresponding to age groups 0--0, 1--4, 5--9, 10--14, ..., 80--84, 85--89, 90--94, 95--99) that should be
% used for the age groups 0, 0.1, 0.2, ..., 99.9

 
%% MP: removed %% cov_InfantVacc_itt = params.InfantVacc;
%% MP: removed %% cov_BirthDose_itt = params.scenario_BirthDose_coverage;
%% MP: dead code - can remove these as they are now checked in country_level_analyses.m:
%%assert(all(params.InfantVacc >= 0) && all(params.InfantVacc <= 1), "HepB3 coverage needs to be 0-1")
%%assert(all(params.scenario_BirthDose_coverage >= 0) && all(params.scenario_BirthDose_coverage <= 1), "BD coverage needs to be 0-1")
%%assert(isequal(size(params.InfantVacc),size(TimeSteps)))
%%assert(isequal(size(params.scenario_BirthDose_coverage),size(TimeSteps)))


assert(isequal(size(Prog),size(zeros(num_disease_states, num_disease_states)))); % Non-Age Specific Prog parameters stored as (from, to)



% Prepare for simulation

% prepare storage containers, for outputs once per year
% breakdowns by age/sex

%% TAM: mini-chunk
[NumSAg_5yr, PrevEAg_of_SAg_5yr] = deal(-99 * ones(2, max(agegroups_5yr), (num_years_simul+1))); 

%% Note that this excludes Vertical transmission
Incid_chronic_all_5yr_approx_no_VertTrans = zeros(2, max(agegroups_5yr), num_years_simul+1);
%% end of mini-chunk

[...
    Tot_Pop_1yr, Prev_Immune_Reactive_1yr, Prev_Chronic_Hep_B_1yr, Prev_Comp_Cirr_1yr, Prev_Decomp_Cirr_1yr, ...
    Prev_Liver_Cancer_1yr, Prev_TDF_treat_1yr, NumSAg_1yr, NumSAg_chronic_1yr, yld_1yr, Prev_Deaths_1yr...
    ] = deal(DUMMY_VALUE * ones(num_sexes, max(agegroups_1yr), num_years_simul+1));
[...
    Incid_chronic_all_1yr_approx,...
    Incid_Deaths_1yr_approx...
    ] = deal(DUMMY_VALUE * ones(num_sexes, max(agegroups_1yr), num_years_simul + 1));

[Prev_HCC_1yr,  NumEAg_chronic_1yr, NumEAg_chronic_acute_1yr] ...
    = deal(DUMMY_VALUE * ones(num_sexes, max(agegroups_1yr), num_years_simul + 1));


%% MP: used as a store 
if(store_results_as_text==1)
    %% Store the following:
    %% - state variables each year (max(agegroups_5yr) * num_disease_states * num_sexes* num_treat_blocks)
    %% - new cases of chronic carriage/yr (neonates, plus 5 yr age gps) = 21
    %% - deaths per year - 5 yr age groups = 20
    ncol_X_to_print_byage = num_disease_states*  num_sexes* num_treat_blocks;
    ncol_results_to_print = max(agegroups_5yr) * ncol_X_to_print_byage + 21 + 20 + 10;
    results_to_print = DUMMY_VALUE * ones(ncol_results_to_print, num_years_simul + 1);
    
    %%X_to_print = DUMMY_VALUE * ones(max(agegroups_5yr)*ncol_X_to_print, num_years_simul + 1);
end

%% Mini TAM: 
%% Get HepB3:
%% Infection stage susceptible (so x1), Age gp - age 6m (so x1), by sex (so x2) and by whether treatment accessible (so x2):
transfer_to_vacc = zeros(1, 1, 2, 2);


%% MP: Magic number 1 is because arrays NewChronicCarriage, moving_btw_states should have same number of dimensions as X()
%% but the first index is single (because not indexing over natural history states).
[...
    NewChronicCarriage, moving_btw_states, ...
    ] = deal(zeros(1, num_age_steps, num_sexes, num_treat_blocks));



% single output per year
%% TAM: extra PAP-model-specific outputs included here:
[Time, RateInfantVacc, RateBirthDoseVacc, RatePeripartumTreatment, ...
    num_births_1yr, NumDecompCirr, NumLiverCancer, ...
    PregnantWomenNeedToScreen, HBVPregnantWomenNeedToEvaluate] = deal(DUMMY_VALUE * ones(1, num_years_simul+1));
 
[num_births_toHbEAgWomenHVL_1yr_approx, num_births_toHbEAgWomenLVL_1yr_approx, num_births_toHbSAgWomenHVL_1yr_approx, num_births_toHbSAgWomenLVL_1yr_approx, ... 
    num_births_1yr_approx, ...
    num_births_chronic_HbEAgWomenHVL_1yr_approx, num_births_chronic_HbEAgWomenLVL_1yr_approx, num_births_chronic_HbSAgWomenHVL_1yr_approx, num_births_chronic_HbSAgWomenLVL_1yr_approx, ... 
    Incid_babies_chronic_1yr_approx, ...
    PeripartumTreatment_HbEAg_HighVL_approx, PeripartumTreatment_HbEAg_LowVL_approx, PeripartumTreatment_HbSAg_HighVL_approx, PeripartumTreatment_HbSAg_LowVL_approx...
    ] = deal(-99 * ones(1, num_years_simul+1));


% ----- Simulation -----

i_dt = 1; % i_dt increase every time i.e. every 0.1 years; goes from 1 to 2101 (length of TimeSteps)
OutputEventNum = 1; % OutputEventNum increase every year; goes from 1 to 212
moving_to_treatment = zeros(size(X));
initiated_treatment = false;
num_babies = 0;

%% MP: This is a little bit spaghetti code.
%% Note that female_multiplier and male_multiplier are both updated later on. They depend on 
%% sex_ratio which is defined below. It would make more sense to initialise them straight after
%% sex_ratio is initialised.
babies_ChronicCarriage = 0;
female_multiplier = 0;
male_multiplier = 0;


moving_to_treatment_by_birthcohort_testing = zeros(size(X));
moving_to_treatment_by_birthcohort_testing_this_timestep = zeros(size(X));


for time = TimeSteps 
       
    
    % Update mortality and fertility rates
    mu = zeros(num_disease_states, num_age_steps, num_sexes, num_treat_blocks);
    % The "1"s below represent the one gender we are considering at a time
    % (we need it so that mu() has the correct dimensions),
    mu(:, :, i_female, :) = repmat(params.MortalityRate_Women(OutputEventNum, MappingFromDataToParam), [num_disease_states 1 num_treat_blocks]);
    mu(:, :, i_male, :) = repmat(params.MortalityRate_Men(OutputEventNum, MappingFromDataToParam), [num_disease_states 1 num_treat_blocks]);
    % params.MortalityRate_Women is a (num_years_simul+1=212) x 21 matrix of mortality rates of 21 age groups (0--0, 1--4, 5--9, 10--14, ..., 90--94, 95--99) for every year from 1890 to 2101
    % OutputEventNum ranges from 1 to (num_years_simul+1)
    % agegroups_1yr = [1 1 1 ... 100 100 100], each number 10 times
    % selected vector copied across disease states and treatments
    % params.fert is a 1000 x (num_years_simul+1) matrix; ages in 0.1 year jumps versus 212 years
    %% MP: Magic number 1:10:end
    fert = params.fert(1:10:end, OutputEventNum);
    %% MP: Magic numbers 100 1
    assert(isequal(size(fert), [100 1]))
    %% MP: Magic number 1
    fert = repmat(fert',num_year_divisions,1);
    fert = fert(:);  % Reshape fert into a 1D vector from a matrix
    %% MP: Magic number 1
    assert(isequal(size(fert),[num_age_steps 1]))
    assert(length(params.net_migration)==(num_years_simul+1))
    net_migration = params.net_migration(OutputEventNum);
    sex_ratio = params.sex_ratios(OutputEventNum);
    assert(isscalar(net_migration))
    net_migration = repmat(net_migration, [num_disease_states num_age_steps num_sexes num_treat_blocks]);


    
    % Compute Outputs once per year
    if rem(time, 1) == 0 % only saves variables in this "for" loop every 10 time steps (or once a year, since dt=0.1)
        Time(OutputEventNum) = time;
 
        
        % Rescale population sizes of each age group and gender
        if (time >= 1950) 
            ModelPop = squeeze(sum(sum(X(i_alive,:,:,:), 1), 4));
            assert(isequal(size(ModelPop),[num_age_steps num_sexes]))
            % sum over disease state of alive people and treatment; ModelPop is 1000 x 2 i.e. age groups versus gender
            assert(isequal(size(params.total_pop_female),[101 152]))
            % params.total_pop_female is a 101 x 152 matrix of 152 years (1950 to 2101) for 101 age groups (0--0, 1--1, 2--2,..., 98--98, 99--99, 100--100)
            assert(isequal(size(params.total_pop_male),[101 152]))
            col_index = time - 1949;
            %% MP: magic numbers 1:100 represent indexes in params.total_pop for ages 0-99
            MontaguPopFemale = params.total_pop_female(1:100,col_index);
            % only want ages 0--99
            MontaguPopMale = params.total_pop_male(1:100,col_index);
            MontaguPop = [MontaguPopFemale MontaguPopMale];
            assert(isequal(size(MontaguPop),[100 num_sexes]))    %% MP: Magic number 100 is number of 1-year age gps 0-99.
            MontaguPopExpand = MontaguPop(agegroups_1yr, :) * dt;
            % agegroups_1yr is a 1 x 1000 double; [1 1 ... 100 100], each number present 10 times
            % expanding MontaguPop from 1 year age steps to 0.1 year age steps
            % each age group repeated 10 times therefore divide each entry by 10
            assert(isequal(size(MontaguPopExpand),[num_age_steps num_sexes]))
            ScalerMat = MontaguPopExpand ./ ModelPop;
            ScalerMat(isnan(ScalerMat)) = 0;
            ScalerMat(isinf(ScalerMat)) = 0;
            pop_scaler = repmat(reshape(ScalerMat, [1 num_age_steps num_sexes]), [num_disease_states 1 1 num_treat_blocks]);
            % MontaguPopExpand is sizes of the current year's population over 0.1 year age steps; a 1000 x 2 matrix of ages versus gender
            % add an extra dimension and duplicate it for each disease state and treatment method
            X = X .* pop_scaler;
            if(time==2020 ||time==2025)
            %disp([min(ScalerMat),max(ScalerMat)])
                disp("Uncomment the line below to show scalarmat")
                %%disp("Scalarmat here:")
                %%disp(ScalerMat)
            end
            % scale all parts of X, including dead people i.e. State i_HBVdeath=11
        end


        


        for k = 1:num_sexes % genders

            for ag = 1:max(agegroups_1yr) % 1:100
 
                if OutputEventNum > 1
                
                    %% MP: Magic numbers: sum over second (age) and 4th (treatment states):
                    state_prev_vec = squeeze(sum(sum(X(:, agegroups_1yr == ag, k, :), 2), 4)); % k is gender
                    assert(isequal(size(state_prev_vec),[num_disease_states 1]))

                    Tot_Pop_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec(i_alive));
                
                    Prev_Immune_Reactive_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_ImmReact);
                
                    Prev_Chronic_Hep_B_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_Chronic);

                    Prev_Comp_Cirr_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_CompCirr);

                    Prev_Decomp_Cirr_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_DecompCirr);

                    Prev_Liver_Cancer_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_HCC);
                
                    Prev_TDF_treat_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_TDFtreat);

                    %% MP TODO: Maybe remove this as the HBV deaths compartment has the same "edge of cliff" thing
                    %% where people aged 99 drop off the model once they turn 100 (so anyone who died age 99 is no longer counted).
                    %% Also HBV deaths is rescaled by ScalerMat, which is a bit funky in the older age groups.
                    %% It might be possible to patch it a bit (say truncate at age 95 to reduce the ScalerMat issue, then cumulatively count deaths
                    %% adding the new 95 yo dead people to an existing cumulative counter of people who are dead who would be 95+ now).
                    Prev_Deaths_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_HBVdeath);

                    NumSAg_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec([i_sAgpos]));
                
                    NumSAg_chronic_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec([i_sAgpos_chronic]));

                    yld_1yr(k, ag, OutputEventNum-1) = sum( state_prev_vec .* params.dwvec' );

                    %% MP: Magic numbers: sum over second (age) and 4th (treatment states):
                    Incid_chronic_all_1yr_approx(k,ag,OutputEventNum-1) = sum(sum(NewChronicCarriage(1, agegroups_1yr == ag, k, :), 2), 4);

                    Incid_Deaths_1yr_approx(k, ag, OutputEventNum-1) = sum(state_prev_vec .* Prog(:, i_HBVdeath));


                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% TAM: PAP mini-chunk 1
                    Prev_HCC_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_HCC);
                    NumEAg_chronic_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec(i_eAgpos_chronic));
                    NumEAg_chronic_acute_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec(i_eAgpos));
                    %% TAM: End of PAP mini-chunk 1
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                
            end % end agegroups_1yr for loop
				
            if (OutputEventNum > 1)

                %% MP: Magic number: the second index "1" represents the age group (0-year-olds)
                assert(Incid_chronic_all_1yr_approx(k, 1, OutputEventNum-1)==0) % 0-year-olds cannot get horizontal chronic infection (see FOI)							
                Incid_chronic_all_1yr_approx(i_female, 1, OutputEventNum-1) = female_multiplier * babies_ChronicCarriage;
                Incid_chronic_all_1yr_approx(i_male, 1, OutputEventNum-1) = male_multiplier * babies_ChronicCarriage;

            end
            

            for ag = 1:max(agegroups_5yr) % 1:20
														
                if OutputEventNum > 1
                    % model results assigned to a particular year at the beginning of that year, after which they are zeroed
                    NumSAg_5yr(k, ag, OutputEventNum-1) = sum(sum(sum(X([i_sAgpos], agegroups_5yr == ag, k, :))));
                
                    if NumSAg_5yr(k, ag, OutputEventNum-1)>0
                        PrevEAg_of_SAg_5yr(k,ag,OutputEventNum-1) = sum(sum(sum(X(2:3, agegroups_5yr == ag, k, :)))) / NumSAg_5yr(k, ag, OutputEventNum-1);
                        % Note that this is prevalence of e+ among s+
                    else
                        assert(NumSAg_5yr(k, ag, OutputEventNum-1)==0)
                        assert(sum(sum(sum(X(i_eAgpos_chronic, agegroups_5yr == ag, k, :))))==0)
                        PrevEAg_of_SAg_5yr(k,ag,OutputEventNum-1) = 0;
                    end

                    Incid_chronic_all_5yr_approx_no_VertTrans(k, ag, OutputEventNum-1) = sum(sum(NewChronicCarriage(1, agegroups_5yr == ag, k, :), 2), 4);

                end
            end

            %%if OutputEventNum > 1
                %% The second index in Incid_chronic_all_5yr_approx() is the age group (age 0-4).
                %% BUG: need to include incidence through child-child transmission:
                %% ORIGINAL CODE:
                %%Incid_chronic_all_5yr_approx(i_female, 1, OutputEventNum-1) = female_multiplier * babies_ChronicCarriage;
                %%Incid_chronic_all_5yr_approx(i_male, 1, OutputEventNum-1) = male_multiplier * babies_ChronicCarriage;
                
                %% FIXED CODE:
                %%Incid_chronic_all_5yr_approx(i_female, 1, OutputEventNum-1) = Incid_chronic_all_5yr_approx(i_female, 1, OutputEventNum-1) + female_multiplier * babies_ChronicCarriage;
                %%Incid_chronic_all_5yr_approx(i_male, 1, OutputEventNum-1) = Incid_chronic_all_5yr_approx(i_male, 1, OutputEventNum-1) + male_multiplier * babies_ChronicCarriage;
            %end
             
        end % end genders for loop


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TAM: PAP mini-chunk 2
        if OutputEventNum > 1
            %% Sum over all dimensions (2,3,4) except natural history:
            state_prev_vec = squeeze(sum(sum(sum(X,2),3),4)); % 15 x 1

            NumDecompCirr(OutputEventNum-1) = state_prev_vec(i_DecompCirr);
            NumLiverCancer(OutputEventNum-1) = state_prev_vec(i_HCC);
            num_births_toHbEAgWomenHVL_1yr_approx(OutputEventNum-1) = births_toHbEAgWomenHighVL;
            num_births_toHbEAgWomenLVL_1yr_approx(OutputEventNum-1) = births_toHbEAgWomenLowVL;
            num_births_toHbSAgWomenHVL_1yr_approx(OutputEventNum-1) = births_toHbSAgWomenHighVL;
            num_births_toHbSAgWomenLVL_1yr_approx(OutputEventNum-1) = births_toHbSAgWomenLowVL;
            num_births_1yr_approx(OutputEventNum-1) = births_Total;
            num_births_chronic_HbEAgWomenHVL_1yr_approx(OutputEventNum-1) = babiesChronic_from_HbEAgWomenHighVL;
            num_births_chronic_HbEAgWomenLVL_1yr_approx(OutputEventNum-1) = babiesChronic_from_HbEAgWomenLowVL;
            num_births_chronic_HbSAgWomenHVL_1yr_approx(OutputEventNum-1) = babiesChronic_from_HbSAgWomenHighVL;
            num_births_chronic_HbSAgWomenLVL_1yr_approx(OutputEventNum-1) = babiesChronic_from_HbSAgWomenLowVL;
            Incid_babies_chronic_1yr_approx(OutputEventNum-1) = babies_ChronicCarriage;
            RateBirthDoseVacc(OutputEventNum-1) = ratebirthdoses;
            RateInfantVacc(OutputEventNum-1) = squeeze(sum(sum(transfer_to_vacc,3),4)) * num_year_divisions;
            PregnantWomenNeedToScreen(OutputEventNum-1) = pregnantWomenNeedToScreen;
            PeripartumTreatment_HbEAg_HighVL_approx(OutputEventNum-1) = num_mothers_PAP_HbEAg_HighVL;
            PeripartumTreatment_HbEAg_LowVL_approx(OutputEventNum-1) = num_mothers_PAP_HbEAg_LowVL;
            PeripartumTreatment_HbSAg_HighVL_approx(OutputEventNum-1) = num_mothers_PAP_HbSAg_HighVL;
            PeripartumTreatment_HbSAg_LowVL_approx(OutputEventNum-1) = num_mothers_PAP_HbSAg_LowVL;
            RatePeripartumTreatment(OutputEventNum-1) = RateOfPAPInitiation;
            HBVPregnantWomenNeedToEvaluate(OutputEventNum-1) = HBVPositivePregnantWomenAtANC;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % MP: Use this to create a text file of all states (but summed into
        % 5 yr age groups to make it tractable).
        % Create a 4D array of size num_disease_states *  20 (5yr age gps
        % 0-99) * num_sexes * num_treat_blocks. 
        if(store_results_as_text==1)
            if OutputEventNum > 1
                % This is the (consolidated into 5 yr age gp) state variable at time t as
                % a 2D array. We want to store this as a row in 
                X_to_print_unshaped = DUMMY_VALUE * ones(max(agegroups_5yr), ncol_X_to_print_byage);
                deaths_to_print = DUMMY_VALUE * ones(1, max(agegroups_5yr));

                for ag = 1:max(agegroups_5yr) % 1:20
                    
                    temp_store = reshape(squeeze(sum(X(:, agegroups_5yr == ag, :, :),2)), [1,ncol_X_to_print_byage]); % k is gender
                    assert(isequal(size(temp_store),[1 60]))
                    X_to_print_unshaped(ag,:) = temp_store;

                    %% Now store incident deaths:
                    
                    death_age_group_indices = (5*(ag-1)+1):(5*(ag-1)+4);
                    temp_store_deaths = squeeze(sum(sum(Incid_Deaths_1yr_approx(:, death_age_group_indices, OutputEventNum-1),2),1));
                    
                    assert(isscalar(temp_store_deaths))
                    deaths_to_print(ag) = temp_store_deaths;
                end
                %% This is the state matrix (grouped into 5 yr age groups, and reshaped into a row vector):
                X_to_print = reshape(X_to_print_unshaped,[1, max(agegroups_5yr)*ncol_X_to_print_byage]);
                
                %% This gives the incidence in neonates (ie via vertical transmission), then in 5 yr age groups (summing over M+F)
                incidence_horizontal_transmission = squeeze(sum(Incid_chronic_all_5yr_approx_no_VertTrans(:, :, OutputEventNum-1),1));
                incidence_to_print = [babies_ChronicCarriage, incidence_horizontal_transmission];
                
                %% BD (standard/MAP/CPAD), infant vacc;
                %% PAP (by EAg+/- and VL), pregnantWomenNeedToScreen=number of women needed to screen to put women on PAP
                %% Treatment is done by state variable
                %% To do: testing to get people on treatment (both standard + birth cohort)
                resources_to_print = [ratebirthdoses, ratebirthdoses_MAP, ratebirthdoses_CPAD, RateInfantVacc(OutputEventNum-1),...
                    num_mothers_PAP_HbEAg_HighVL, num_mothers_PAP_HbEAg_LowVL, num_mothers_PAP_HbSAg_HighVL, num_mothers_PAP_HbSAg_LowVL, pregnantWomenNeedToScreen,...
                    number_starting_treatment_to_print];



                results_to_print(:,OutputEventNum-1) = [X_to_print, incidence_to_print, deaths_to_print, resources_to_print]';
                

            end

        end

        %if((OutputEventNum>115 && OutputEventNum<120) || OutputEventNum==160) 
            %%fprintf("AYear %5d Time %5d X%10.8f D2 %10.8f age2 %10.8f age3 %10.8f M%10.8f Treat%10.8f Pop %12.6f\n",1889+OutputEventNum, time, sum(X(1, agegroups_5yr == 1, 1, 1)), sum(X(2, agegroups_5yr == 1, 1, 1)), sum(X(1, agegroups_5yr == 2, 1, 1)), sum(X(1, agegroups_5yr == 3, 1, 1)), sum(X(1, agegroups_5yr == 1, 2, 1)), sum(X(1, agegroups_5yr == 1, 1, 2)), sum(sum(sum(sum(X(i_alive, :, :, :), 2), 4),1),3));
        %    if(OutputEventNum==116)
        %        writematrix(squeeze(X(i_alive, 1, 1, 1)),fullfile(basedir,'temp.csv'));
        %    end
            %disp([1889+OutputEventNum, sum(sum(sum(sum(X(i_alive, :, :, :), 2), 4),1),3)])
        %end

        %%disp([1889+OutputEventNum, sum(sum(sum(sum(X(i_HBVdeath, :, :, :), 2), 4),1),3)])

        if OutputEventNum > 1
            
            assert(isscalar(num_babies))
            num_births_1yr(OutputEventNum-1) = num_babies;
            
        end
        
        num_babies = 0;
        
        % Update counter
        OutputEventNum = OutputEventNum + 1;
        % increases every year
    end % end "rem(time, 1) == 0" if statement
    
    
    
    % Horizontal Transmission (infection and Chronic Carriage)
    
    FOI = zeros(1, num_age_steps, num_sexes, num_treat_blocks);
    % Calculate number of people in different age groups for use in FOI
    % denominator:
    n_child_1y_5y  = sum(sum(sum(sum(X(i_alive, i1y:(i5y - 1), :, :)))));
    n_child_1y_15y = sum(sum(sum(sum(X(i_alive, i1y:(i15y - 1), :, :)))));
    n_pop_5y_andabove = sum(sum(sum(sum(X(i_alive, i5y:end, :, :)))));
    % i: Transmission Between 1y-5y olds
    FOI(1, i1y:(i5y - 1), :, :) = ...
        beta_U5_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_notEagpos_notreat, i1y:(i5y - 1), :, :))))) / n_child_1y_5y ...
        + beta_U5_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos, i1y:(i5y - 1), :, :))))) / n_child_1y_5y;
    
    % ii: Transmission between 1-15 year olds
    FOI(1, i1y:(i15y - 1), :, :) = FOI(1, i1y:(i15y - 1), :, :) + ...
        beta_1to15_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_notEagpos_notreat, i1y:(i15y - 1), :, :))))) / n_child_1y_15y ...
        + beta_1to15_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos, i1y:(i15y - 1), :, :))))) / n_child_1y_15y;
    
    % iii: Transmission Between 5+ and Adults (Assuming equal risks for all persons 5y-100y)
    FOI(1, i5y:end, :, :) = FOI(1, i5y:end, :, :) + ...
        beta_5plus_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_notEagpos_notreat, i5y:end, :, :))))) / n_pop_5y_andabove ...
        + beta_5plus_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos, i5y:end, :, :))))) / n_pop_5y_andabove;
    
    
    % Disease Progression
    next_X = X;
    for tr = 1:length(Transactions.From)
        %%transaction_vals = Transactions.Values{tr};
        %%transaction_vals = transaction_vals(:);
        moving_btw_states(1, :, :, :) = X(Transactions.From(tr), :, :, :) .* Transactions.Values{tr};
        next_X(Transactions.From(tr), :, :, :) = next_X(Transactions.From(tr), :, :, :) + dt * ( -moving_btw_states ); % move people out of "from" state
        next_X(Transactions.To(tr), :, :, :) = next_X(Transactions.To(tr), :, :, :)  + dt * ( +moving_btw_states ); % move people into "to" state
    end % end Disease Progression for loop
    % multiply by dt because quantity is calculated every 0.1 years therefore needs to be divided by 10


        

    % (Time-dependent) Baseline Transition to TDF-Treatment
    assert(squeeze(sum(sum(sum(sum(X([i_3TCtreat i_3TCfailed], :, :, 1),1),2),3),4))==0)
    % (Time-dependent) Baseline Transition to TDF-Treatment
    
    %% SERNIK
    birth_cohort_testing_start = 2026;
    birth_cohort_testing_end = 2029;
    if (time >= birth_cohort_testing_start && time <= birth_cohort_testing_end)
        if(scenario_CohortTesting==I_COHORT_TEST)
            %%case I_NO_COHORT_TEST
            disp("Running birth cohort testing and treatment")
            disp(time)
            %%case I_COHORT_TEST
            if(time == birth_cohort_testing_start)
                i_cohortage_min = round((birth_cohort_testing_start-1992)/dt); % Individuals born before 1992
                i_cohortage_max = round(60/dt); % Individuals born before 1992
                % Check I haven't accidentally made the min age>max age:
                assert(i_cohortage_max>i_cohortage_min)

                birth_cohort_coverage = 0.8;
                %% Note we should use next_X rather than X here:
                moving_to_treatment_by_birthcohort_testing_per_timestep(i_treateligible, i_cohortage_min:i_cohortage_max, :, :) = dt * birth_cohort_coverage * next_X(i_treateligible, i_cohortage_min:i_cohortage_max, :, :)/(birth_cohort_testing_end-birth_cohort_testing_start); 
                disp("Eligible for cohort treatment")
                disp(sum(sum(sum(sum(moving_to_treatment_by_birthcohort_testing_per_timestep,1),2),3),4))
            end
            disp("At time")
            disp(time)
            
            assert(time<=birth_cohort_testing_end);
            i_birth_cohort_offset = round((time-birth_cohort_testing_start)/0.1);
            moving_to_treatment_by_birthcohort_testing_this_timestep(i_treateligible, (i_cohortage_min+i_birth_cohort_offset):i_cohortage_max, :, :) ...
                = moving_to_treatment_by_birthcohort_testing_per_timestep(i_treateligible, i_cohortage_min:(i_cohortage_max-i_birth_cohort_offset), :, :);
            %% Set the earlier age group elements to zero if needed:
            if(i_birth_cohort_offset>0)
                moving_to_treatment_by_birthcohort_testing_this_timestep(i_treateligible, i_cohortage_min:(i_cohortage_min+i_birth_cohort_offset-1), :, :) ...
                    = zeros(length(i_treateligible), i_birth_cohort_offset, num_sexes, num_treat_blocks);
            end
            %% Ensure we never go below 0:
            moving_to_treatment_by_birthcohort_testing_this_timestep(moving_to_treatment_by_birthcohort_testing_this_timestep>next_X) = next_X(moving_to_treatment_by_birthcohort_testing_this_timestep>next_X);
            next_X(i_treateligible, :, :, :) = next_X(i_treateligible, :, :, :) - moving_to_treatment_by_birthcohort_testing_this_timestep(i_treateligible, :, :, :);
            next_X(i_TDFtreat, :, :, :) = next_X(i_TDFtreat, :, :, :) + sum(moving_to_treatment_by_birthcohort_testing_this_timestep, 1);
            %%otherwise
            %%    disp("Error: Unknown value for scenario_CohortTesting. Exiting")
            %%    return
        end 
    end
    %%if (time >= birth_cohort_testing_start && time <= birth_cohort_testing_end)

    if (time >= treat_start_year && HAS_TREATMENT~=0)
    % 2016 must be the first year with nonzero treatment 
    % therefore start treating from 2015.9 onwards since prevalence is recorded at the top of the loop

        if ~initiated_treatment
            num_in_treatment = sum(sum(sum(sum(X(i_TDFtreat, :, :, :),1),2),3),4);
            assert(num_in_treatment==0) % no one is in treatment
            prev_pop = sum(sum(sum(sum(X([i_sAgpos], :, :, :),1),2),3),4); 

            total_num_to_move_to_treat = treat_coverage_in_2016 * prev_pop;
            eligible_pop = sum(sum(sum(sum(X(i_treateligible, :, :, :),1),2),3),4); 
            assert(total_num_to_move_to_treat<eligible_pop)
            scaling_num = total_num_to_move_to_treat / eligible_pop;
            next_X(i_treateligible,:,:,:)=next_X(i_treateligible,:,:,:) - X(i_treateligible,:,:,:) * scaling_num;

            
            % Every compartment in the eligible-for-treatment states in next_X must have a number subtracted from it 
            % such that the total number subtracted from the eligible-for-treatment states is in_treatment_2016
            % i.e. in_treatment_2016 = sum(sum(sum(sum(X(i_treateligible, :, :, :),1),2),3),4) * scaling_num = sum(sum(sum(sum(X(i_treateligible, :, :, :) * scaling_num,1),2),3),4)
            % Hence, scaling_num scales each compartment in X(i_treateligible, :, :, :) such that X(i_treateligible,:,:,:) * scaling_num subtracts the same proportion of people from each compartment in each of the eligible-for-treatment states in order to subtract a total of in_treatment_2016 from the eligible-for-treatment states.
            next_X(i_TDFtreat,:,:,:) = next_X(i_TDFtreat,:,:,:) + sum(X(i_treateligible,:,:,:) * scaling_num,1);

            num_in_treatment = sum(sum(sum(sum(next_X(i_TDFtreat, :, :, :),1),2),3),4);

            %% Since we just checked that the number in natural history state i_TDFtreat was zero before treatemnt started in the smulation, this represents the number of people starting treatment at this timestep.
            number_starting_treatment_to_print = num_in_treatment;
            
            %% eligible_pop includes those on treatment here:
            eligible_pop = sum(sum(sum(sum(X([i_treateligible, i_TDFtreat], :, :, :),1),2),3),4); 
            assert(num_in_treatment/eligible_pop >= treat_coverage_in_2016)
            % treatment coverage amongst treatment-eligible people will be greater than treatment coverage amongst HBsAg+ people, except if treatment coverage is 0
            treat_coverage_2016 = num_in_treatment / eligible_pop; %% MP: CHECK WITH SHEVANTHI - THIS IS CURRENTLY DEAD CODE.

            initiated_treatment = true;
        else
            assert(initiated_treatment) % ensure that, each time this code is encountered, treatment has already been initiated
            
            if (time<=treatment_rate_params.t_treatment_scaleup_start)
                treatment_rate = treatment_rate_params.Treatmentrate_2016;
            elseif (time>=treatment_rate_params.t_treatment_scaleup_end)
                treatment_rate = treatment_rate_params.Treatmentrate_final;
            else
                treatment_rate = treatment_rate_params.Treatmentrate_2016 + (treatment_rate_params.Treatmentrate_final - treatment_rate_params.Treatmentrate_2016) * (time-treatment_rate_params.t_treatment_scaleup_start)/(treatment_rate_params.t_treatment_scaleup_end - treatment_rate_params.t_treatment_scaleup_start);
            end
            assert(treatment_rate>=0)
            
            moving_to_treatment(i_treateligible, :, :, :) = X(i_treateligible, :, :, :) .* treatment_rate;
            next_X(i_treateligible, :, :, :) = next_X(i_treateligible, :, :, :) + dt * ( -moving_to_treatment(i_treateligible, :, :, :) );
            next_X(i_TDFtreat, :, :, :) = next_X(i_TDFtreat, :, :, :) + dt * ( +sum(moving_to_treatment, 1) );
            assert(max(moving_to_treatment(:))>=0)

            number_starting_treatment_to_print = squeeze(sum(sum(sum(sum(moving_to_treatment, 1), 2), 3), 4));
            assert(isscalar(number_starting_treatment_to_print))

        end
    end % end treatment if statement
    

    
    
    % Infection process
    NewInfections = X(i_Susc, :, :, :) .* FOI;
    % number of susceptibles times FOI i.e. number of new infections within population, excluding babies (since FOI is 0 for babies)
    % 1 x num_age_steps x 2 x 2 double i.e. 1 x 1000 x 2 x 2 double
    NewChronicCarriage = NewInfections .* p_ChronicCarriage;
    SevereAcute = NewInfections * theta;
    NonsevereAcute = NewInfections - SevereAcute;
    
    % Transitions dependent on a state that does not have a number and is therefore not in Prog or Transactions
    % multiply by dt, since FOI is an annual rate
    next_X(i_Susc, :, :, :) = next_X(i_Susc, :, :, :) + dt * ( -NewInfections );
    next_X(i_NonSevAcute, :, :, :) = next_X(i_NonSevAcute, :, :, :) + dt * ( +NonsevereAcute );
    next_X(i_SevereAcute, :, :, :) = next_X(i_SevereAcute, :, :, :) + dt * ( +SevereAcute );
    
    
    % Infant vaccination HepB3:
    % Do not multiply by dt, since one is vaccinating scenario_HepB3coverage(i_dt)% of people in next_X(1, i6mo, :, :), 
    % after which this cohort ages and moves to the next age bin
    % if divides all babies born in a year into 10 groups and vaccinatates scenario_HepB3coverage(i_dt)% of each group, 
    % then one will have vaccinated scenario_HepB3coverage(i_dt)% of all babies born in that year
    transfer_to_vacc = scenario_HepB3coverage(i_dt) * next_X(i_Susc, i6mo, :, :) * params.Efficacy_InfantVacc; % the 0.95 represent a take-type vaccine efficacy of 95%.
    next_X(i_Susc, i6mo, :, :) = next_X(i_Susc, i6mo, :, :) - transfer_to_vacc;
    next_X(i_Immune, i6mo, :, :) = next_X(i_Immune, i6mo, :, :) + transfer_to_vacc;
    
        
    
    % Natural Mortality
    % Do not apply background mortality to the HBV deaths state, since we want people in all countries to be treated as if they would have lived until 84 if they had not died of HBV. This is done outside of the model in the main script. 
    mu(i_HBVdeath, :, :, :, :)=0.0;
    next_X = next_X + dt * ( -next_X .* mu );
    % in the "for time=TimeSteps", which runs 10 times per year, therefore divide effect of mu by 10
    net_migration(i_HBVdeath, :, :, :, :)=0.0;
    next_X = next_X + dt * ( +next_X .* net_migration );
    % "+" because net_migration = (number of immigrants - number of emigrants) / population size
    
    % Update Stocks
    X = next_X;
    
    % Now age everyone (second index is the age index with 1=newborn in this
    % timestep).
    X(:, 2:num_age_steps, :, :) = X(:, 1:(num_age_steps - 1), :, :);
    X(:, 1, :, :) = 0; % set number of new babies (the age index "1") to 0 (babies will be born next)
    
    

    %% KAWA/TAM 2:
    % fill-out with new births in this time-step:
    %% MP: Magic numbers 1 and 4 mean sum over the listed natural history states and treatment states
    births_toNonInfectiousWomen = sum( fert' .* sum(sum(X([i_Susc i_Immune], :, i_female, :), 1), 4) ); % Susecptible, Immune
    %%births_toHbEAgWomen = sum(fert' .* sum(sum(X(i_eAgpos, :, i_female, :), 1), 4)); % Immune Tolerant, Immune Reactive
    %%births_toHbSAgWomen = sum(fert' .* sum(sum(X(i_sAgpos_notEagpos_notreat, :, i_female, :), 1), 4)); % All other stages (other infected women)
    
    %% In the PAP model these are incorporated in births_toNonInfectiousWomen. 
    %% As treatment will reduce VL we don't bother stratifying by high/low VL here:
    births_toTrWomen = sum(fert' .* sum(sum(X([i_TDFtreat i_3TCtreat], :, i_female, :), 1), 4)); % Women on Treatment
    
    births_toHbEAgWomenHighVL = PAP_VL_params.FracEPosHighVL    * sum(fert' .* sum(sum(X(i_eAgpos, :, i_female, :), 1), 4)); %Immune Tolerant, Immune Reactive, Acute
    births_toHbEAgWomenLowVL = (1-PAP_VL_params.FracEPosHighVL) * sum(fert' .* sum(sum(X(i_eAgpos, :, i_female, :), 1), 4)); %Immune Tolerant, Immune Reactive, Acute
    births_toHbSAgWomenHighVL = PAP_VL_params.FracSPosHighVL    * sum(fert' .* sum(sum(X(i_sAgpos_notEagpos_notreat, :, i_female, :), 1), 4)); %All other stages (other infected women not on treatment)
    births_toHbSAgWomenLowVL = (1-PAP_VL_params.FracSPosHighVL) * sum(fert' .* sum(sum(X(i_sAgpos_notEagpos_notreat, :, i_female, :), 1), 4)); %All other stages (other infected women not on treatment)
    births_Total = births_toNonInfectiousWomen + births_toTrWomen + births_toHbEAgWomenHighVL + ...
        births_toHbEAgWomenLowVL + births_toHbSAgWomenHighVL + births_toHbSAgWomenLowVL;
    assert(isscalar(births_Total))
    
    %%births_Total = births_toNonInfectiousWomen + births_toHbEAgWomen + births_toHbSAgWomen + births_toTrWomen;
    
    num_babies = num_babies + dt * births_Total;
    



    % babies_ChronicCarriage = p_ChronicCarriage(1, 1, 1, 1) * ( ... % a 1 x 1 double
    %     ...
    %     births_toHbSAgWomen * (1 - scenario_BDcoverage(i_dt) - scenario_BDcoverage_fromMAP_CPAD(i_dt))...
    %     * p_VerticalTransmission_HbSAg_NoBD ...
    %     + births_toHbSAgWomen * scenario_BDcoverage(i_dt) * p_VerticalTransmission_HbSAg_BirthDoseVacc ...
    %     + births_toHbSAgWomen * scenario_BDcoverage_fromMAP_CPAD(i_dt) * p_VerticalTransmission_HbSAg_BirthDose_MAP_CPAD ...
    %     ...
    %     + births_toHbEAgWomen * (1 - scenario_BDcoverage(i_dt)) * p_VerticalTransmission_HbEAg_NoBD ...
    %     + births_toHbEAgWomen * scenario_BDcoverage(i_dt) * p_VerticalTransmission_HbEAg_BirthDoseVacc ...
    %     + births_toHbEAgWomen * scenario_BDcoverage_fromMAP_CPAD(i_dt) * p_VerticalTransmission_HbEAg_BirthDose_MAP_CPAD ...
    %     ...
    %     + births_toTrWomen * (1 - scenario_BDcoverage(i_dt)) * p_VerticalTransmission_Tr_NoBD ...
    %     + births_toTrWomen * scenario_BDcoverage(i_dt) * p_VerticalTransmission_Tr_BirthDoseVacc ...
    %     + births_toTrWomen * scenario_BDcoverage_fromMAP_CPAD(i_dt) * p_VerticalTransmission_Tr_BirthDose_MAP_CPAD ...
    %     );


    %% KAWA:
    %% Interventions:
    %% BD (standard BD, MAP or CPAD)
    %% PAP or treatment (treatment is a separate compartment so dealt with separately).
    
    prop_no_BD_this_timestep = (1 - scenario_BDcoverage(i_dt) - scenario_BDcoverage_fromMAP(i_dt) - scenario_BDcoverage_fromCPAD(i_dt));
    
    %% Number of babies born with chronic Hep B from women EAg+ with high VL (and not on treatment):
    babiesChronic_from_HbEAgWomenHighVL = p_ChronicCarriage(1, 1, 1, 1) * births_toHbEAgWomenHighVL * ...
        ( ...
            prop_no_BD_this_timestep * (1-PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgHighVL(i_dt)) * p_VertTrans_HbEAgHighVL_NoIntv ...
            + prop_no_BD_this_timestep * PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgHighVL(i_dt) * p_VertTrans_HbEAgHighVL_PAP  ...
            + scenario_BDcoverage(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL(i_dt)) * p_VertTrans_HbEAgHighVL_BD  ...
            + scenario_BDcoverage(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL(i_dt) * p_VertTrans_HbEAgHighVL_BD_PAP  ...
            + scenario_BDcoverage_fromMAP(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL(i_dt)) * p_VertTrans_HbEAgHighVL_MAP  ...
            + scenario_BDcoverage_fromMAP(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL(i_dt) * p_VertTrans_HbEAgHighVL_MAP_PAP  ...
            + scenario_BDcoverage_fromCPAD(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL(i_dt)) * p_VertTrans_HbEAgHighVL_CPAD  ...
            + scenario_BDcoverage_fromCPAD(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL(i_dt) * p_VertTrans_HbEAgHighVL_CPAD_PAP  ...
        );
    % HVL HbEAg+ pregnant women who get PAP = (HVL HbEAg+ pregnant women who get PAP and whose babies do not receive BD) + (HVL HbEAg+ pregnant women who get PAP and whose babies receive BD),
    % where number of births is used to approximate number of mothers (a woman can have twins, which makes number of births not equal to number of mothers)
    num_mothers_PAP_HbEAg_HighVL = births_toHbEAgWomenHighVL * (prop_no_BD_this_timestep*PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgHighVL(i_dt) ...
        + (1-prop_no_BD_this_timestep) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL(i_dt) );


    %% Number of babies born with chronic Hep B from women EAg+ with low VL (and not on treatment):
    babiesChronic_from_HbEAgWomenLowVL = p_ChronicCarriage(1, 1, 1, 1) * births_toHbEAgWomenLowVL * ...
        ( ...
            prop_no_BD_this_timestep * (1-PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgLowVL(i_dt)) * p_VertTrans_HbEAgLowVL_NoIntv ...
            + prop_no_BD_this_timestep * PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgLowVL(i_dt) * p_VertTrans_HbEAgLowVL_PAP  ...
            + scenario_BDcoverage(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL(i_dt)) * p_VertTrans_HbEAgLowVL_BD  ...
            + scenario_BDcoverage(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL(i_dt) * p_VertTrans_HbEAgLowVL_BD_PAP  ...
            + scenario_BDcoverage_fromMAP(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL(i_dt)) * p_VertTrans_HbEAgLowVL_MAP  ...
            + scenario_BDcoverage_fromMAP(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL(i_dt) * p_VertTrans_HbEAgLowVL_MAP_PAP  ...
            + scenario_BDcoverage_fromCPAD(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL(i_dt)) * p_VertTrans_HbEAgLowVL_CPAD  ...
            + scenario_BDcoverage_fromCPAD(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL(i_dt) * p_VertTrans_HbEAgLowVL_CPAD_PAP  ...
        );
    % LVL HbEAg+ pregnant women who get PAP = (LVL HbEAg+ pregnant women who get PAP and whose babies do not receive BD) + (LVL HbEAg+ pregnant women who get PAP and whose babies receive BD),
    % where number of births is used to approximate number of mothers (a woman can have twins, which makes number of births not equal to number of mothers)
    num_mothers_PAP_HbEAg_LowVL = births_toHbEAgWomenLowVL * (prop_no_BD_this_timestep*PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgLowVL(i_dt) ...
        + (1-prop_no_BD_this_timestep) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL(i_dt) );


    %%%% CIASTECKO 2
    %% Number of babies born with chronic Hep B from women SAg+ (EAg-) with high VL (and not on treatment):
    babiesChronic_from_HbSAgWomenHighVL = p_ChronicCarriage(1, 1, 1, 1) * births_toHbSAgWomenHighVL * ...
        ( ...
            prop_no_BD_this_timestep * (1-PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgHighVL(i_dt)) * p_VertTrans_HbSAgHighVL_NoIntv ...
            + prop_no_BD_this_timestep * PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgHighVL(i_dt) * p_VertTrans_HbSAgHighVL_PAP  ...
            + scenario_BDcoverage(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL(i_dt)) * p_VertTrans_HbSAgHighVL_BD  ...
            + scenario_BDcoverage(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL(i_dt) * p_VertTrans_HbSAgHighVL_BD_PAP  ...
            + scenario_BDcoverage_fromMAP(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL(i_dt)) * p_VertTrans_HbSAgHighVL_MAP  ...
            + scenario_BDcoverage_fromMAP(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL(i_dt) * p_VertTrans_HbSAgHighVL_MAP_PAP  ...
            + scenario_BDcoverage_fromCPAD(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL(i_dt)) * p_VertTrans_HbSAgHighVL_CPAD  ...
            + scenario_BDcoverage_fromCPAD(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL(i_dt) * p_VertTrans_HbSAgHighVL_CPAD_PAP  ...
        );

    % HVL HbSAg+ pregnant women who get PAP = (HVL HbSAg+ pregnant women who get PAP and whose babies do not receive BD) + (HVL HbSAg+ pregnant women who get PAP and whose babies receive BD),
    % where number of births is used to approximate number of mothers (a woman can have twins, which makes number of births not equal to number of mothers)
    num_mothers_PAP_HbSAg_HighVL = births_toHbSAgWomenHighVL * (prop_no_BD_this_timestep*PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgHighVL(i_dt) ...
        + (1-prop_no_BD_this_timestep) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL(i_dt) );

    %%%% CIASTECKO 3
    %% Number of babies born with chronic Hep B from women SAg+ (EAg-) with low VL (and not on treatment):
    babiesChronic_from_HbSAgWomenLowVL = p_ChronicCarriage(1, 1, 1, 1) * births_toHbSAgWomenLowVL * ...
        ( ...
            prop_no_BD_this_timestep * (1-PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgLowVL(i_dt)) * p_VertTrans_HbSAgLowVL_NoIntv ...
            + prop_no_BD_this_timestep * PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgLowVL(i_dt) * p_VertTrans_HbSAgLowVL_PAP  ...
            + scenario_BDcoverage(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL(i_dt)) * p_VertTrans_HbSAgLowVL_BD  ...
            + scenario_BDcoverage(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL(i_dt) * p_VertTrans_HbSAgLowVL_BD_PAP  ...
            + scenario_BDcoverage_fromMAP(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL(i_dt)) * p_VertTrans_HbSAgLowVL_MAP  ...
            + scenario_BDcoverage_fromMAP(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL(i_dt) * p_VertTrans_HbSAgLowVL_MAP_PAP  ...
            + scenario_BDcoverage_fromCPAD(i_dt) * (1-PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL(i_dt)) * p_VertTrans_HbSAgLowVL_CPAD  ...
            + scenario_BDcoverage_fromCPAD(i_dt) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL(i_dt) * p_VertTrans_HbSAgLowVL_CPAD_PAP  ...
        );

    % LVL HbSAg+ pregnant women who get PAP = (LVL HbSAg+ pregnant women who get PAP and whose babies do not receive BD) + (LVL HbSAg+ pregnant women who get PAP and whose babies receive BD),
    % where number of births is used to approximate number of mothers (a woman can have twins, which makes number of births not equal to number of mothers)
    num_mothers_PAP_HbSAg_LowVL = births_toHbSAgWomenLowVL * (prop_no_BD_this_timestep*PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgLowVL(i_dt) ...
        + (1-prop_no_BD_this_timestep) * PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL(i_dt) );

    % Number of pregnant women who get PAP in this timestep:
    RateOfPAPInitiation = num_mothers_PAP_HbEAg_HighVL + num_mothers_PAP_HbEAg_LowVL + num_mothers_PAP_HbSAg_HighVL + num_mothers_PAP_HbSAg_LowVL;
    

    

    %%p_VertTrans_HbEAgHighVL_Treat      = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_Treat * p_VertTrans_HbEAgHighVL_NoIntv;
    %%p_VertTrans_HbEAgHighVL_BD_Treat   = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_BD_Treat * p_VertTrans_HbEAgHighVL_NoIntv;
    %%p_VertTrans_HbEAgHighVL_MAP_Treat  = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_MAP_Treat * p_VertTrans_HbEAgHighVL_NoIntv;
    %%p_VertTrans_HbEAgHighVL_CPAD_Treat = PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_CPAD_Treat * p_VertTrans_HbEAgHighVL_NoIntv;


    %% Now women on treatment:
    babiesChronic_from_HbEAgTrWomen = p_ChronicCarriage(1, 1, 1, 1) * births_toTrWomen * ...
        (...
            prop_no_BD_this_timestep * p_VertTrans_HbEAg_Treat ...
            + scenario_BDcoverage(i_dt) * p_VertTrans_HbEAg_Treat_BD ...
             + scenario_BDcoverage_fromMAP(i_dt) * p_VertTrans_HbEAg_Treat_MAP ...
             + scenario_BDcoverage_fromCPAD(i_dt) * p_VertTrans_HbEAg_Treat_CPAD ...
             );

        % (1 - scenario_BDcoverage(i_dt) - scenario_BDcoverage_fromMAP_CPAD(i_dt)) * p_VerticalTransmission_HbSAg_NoBD ...
        % + births_toHbSAgWomen * scenario_BDcoverage(i_dt) * p_VerticalTransmission_HbSAg_BD ...
        % + births_toHbSAgWomen * scenario_BDcoverage_fromMAP_CPAD(i_dt) * p_VerticalTransmission_HbSAg_BirthDose_MAP_CPAD ...

    ratebirthdoses = births_Total * scenario_BDcoverage(i_dt);
    ratebirthdoses_MAP = births_Total * scenario_BDcoverage_fromMAP(i_dt);
    ratebirthdoses_CPAD = births_Total * scenario_BDcoverage_fromCPAD(i_dt);
    

    babies_ChronicCarriage = babiesChronic_from_HbEAgWomenHighVL + babiesChronic_from_HbEAgWomenLowVL + ...
        babiesChronic_from_HbSAgWomenHighVL + babiesChronic_from_HbSAgWomenLowVL + babiesChronic_from_HbEAgTrWomen;

    babies_NotChronicCarriage = births_Total - babies_ChronicCarriage;


    assert(isscalar(babies_ChronicCarriage))
    assert(isscalar(babies_NotChronicCarriage))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TAM: PAP chunk 2 - should replace above.


    % number of chronic babies born to HVL HBeAg+ pregnant women =
    %     (probability of infection becoming chronic in babies) * (number of births to HVL HBeAg+ pregnant women) *
    %     (
    %       probability of a baby that does not receive BD born to a HVL HBeAg+ pregnant woman who does not receive PAP being infected
    %     + probability of a baby that does not receive BD born to a HVL HBeAg+ pregnant woman who receives PAP being infected
    %     + probability of a baby that receives BD born to a HVL HBeAg+ pregnant woman who does not receive PAP being infected
    %     + probability of a baby that receives BD born to a HVL HBeAg+ pregnant woman who receives PAP being infected
    %     )

                                                                     
    



    %% In the PAP model, but dead code.
    % tmp_MTCTRate_SPosPregWomen = (babiesChronic_from_HbSAgWomenHighVL + babiesChronic_from_HbSAgWomenLowVL) / (births_toHbSAgWomenLowVL + births_toHbSAgWomenHighVL);
    % 
    % tmp_MTCTRate_EPosPregWomen = (babiesChronic_from_HbEAgWomenHighVL + babiesChronic_from_HbEAgWomenLowVL) / (births_toHbEAgWomenLowVL + births_toHbEAgWomenHighVL);
    % 
    % tmp_MTCTRate_AllPosPregWomen = ...
    %                 (babiesChronic_from_HbEAgWomenHighVL + babiesChronic_from_HbEAgWomenLowVL + babiesChronic_from_HbSAgWomenHighVL + babiesChronic_from_HbSAgWomenLowVL) / ...
    %                 (births_toHbSAgWomenLowVL + births_toHbSAgWomenHighVL + births_toHbEAgWomenLowVL + births_toHbEAgWomenHighVL);

    pregnantWomenNeedToScreen = births_Total * max([PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL(i_dt),...
                                                    PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL(i_dt),...
                                                    PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL(i_dt),...
                                                    PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL(i_dt), ...
                                                    PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgHighVL(i_dt),...
                                                    PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgLowVL(i_dt),...
                                                    PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgHighVL(i_dt),...
                                                    PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgLowVL(i_dt)]);


    %% TO CHECK - this excludes women on treatment right now.
    HBVPositivePregnantWomenAtANC = (births_toHbSAgWomenHighVL + births_toHbSAgWomenLowVL + births_toHbEAgWomenHighVL + births_toHbEAgWomenLowVL) ...        %added 13/8/19
                                    * max([PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL(i_dt),...
                                           PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL(i_dt),...
                                           PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL(i_dt),...
                                           PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL(i_dt),...
                                           PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgHighVL(i_dt),...
                                           PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgLowVL(i_dt),...
                                           PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgHighVL(i_dt),...
                                           PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgLowVL(i_dt)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TAM: End of PAP Chunk 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    female_multiplier = 1 / (1 + sex_ratio);
    male_multiplier = sex_ratio / (1 + sex_ratio);
    % sex_ratio is number of male births per one female birth
    % 1 / (1 + sex_ratio) + sex_ratio / (1 + sex_ratio) = 1, hence total number of babies not changed
    % male births -> 0 => sex_ratio -> 0 => male_multiplier -> 0/1 = 0
    % male births -> infinity => sex_ratio -> infinity => male_multiplier -> 1
    % female births -> 0 => sex_ratio -> infinity => female_multiplier -> 0
    % female births -> infinity => sex_ratio -> 0 => female_multiplier -> 1
    
    X(i_Susc, 1, i_female, i_notreat) = female_multiplier * dt * babies_NotChronicCarriage;  % Suscpetible babies
    X(i_ImmTol, 1, i_female, i_notreat) = female_multiplier * dt * babies_ChronicCarriage;     % Babies with chronic carriage
    
    X(i_Susc, 1, i_male, i_notreat) = male_multiplier * dt * babies_NotChronicCarriage;  % Suscpetible babies
    X(i_ImmTol, 1, i_male, i_notreat) = male_multiplier * dt * babies_ChronicCarriage;     % Babies with chronic carriage


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TAM: PAP Chunk 3:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% In the PAP model, but dead code.
    % Prevalence HBsAg among pregnant women
    % tmp_PrevalenceAmongPregnantWomen = ...
    %     (births_toHbSAgWomenLowVL + births_toHbSAgWomenHighVL + births_toHbEAgWomenLowVL + births_toHbEAgWomenHighVL) / births_Total ; 
    % 
    % tmp_EPrevalenceAmongPregnantWomen = ...
    %     (births_toHbEAgWomenLowVL + births_toHbEAgWomenHighVL) / (births_toHbSAgWomenLowVL + births_toHbSAgWomenHighVL + births_toHbEAgWomenLowVL + births_toHbEAgWomenHighVL) ;     
    % 
    % % Mean year of birth of pregnant women
    % tmp_MeanYearOfBirthOfPregnantWomen = time - ... 
    %     (sum(ages .* fert' .* squeeze(sum(sum(X([1:10 12 13:15],:,1,:),1),4))) / ...
    %         sum(fert' .* squeeze(sum(sum(X([1:10 12 13:15],:,1,:),1),4))));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TAM: End of Chunk 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % increment the timestep index
    i_dt = i_dt + 1;
    % increases every 0.1 years

    %% PARANOID ANDROID
    % if((OutputEventNum>115 && OutputEventNum<120) || OutputEventNum==160) 
    %     fprintf("BYear %5d Time %5d X1=%10.8f x2=%10.8f x%10.8f Pop %12.6f\n",1889+OutputEventNum, time, sum(X(1, agegroups_5yr == 1, 1, 1)), sum(X(2, agegroups_5yr == 1, 1, 1),2), sum(X(4, agegroups_5yr == 3, 2, 1)), sum(sum(sum(sum(X(i_alive, :, :, :), 2), 4),1),3));
    % end
    
end % end "time = TimeSteps" for loop

output.Time = Time; % 1 x (num_years_simul + 1)
output.Tot_Pop_1yr = Tot_Pop_1yr; % 2 x 100 x (num_years_simul + 1)
output.num_births_1yr = num_births_1yr; % 1 x (num_years_simul + 1)
output.Incid_chronic_all_1yr_approx = Incid_chronic_all_1yr_approx; % 2 x 100 x (num_years_simul + 1)
output.Prev_Immune_Reactive_1yr = Prev_Immune_Reactive_1yr; % 2 x 100 x (num_years_simul + 1)
output.Prev_Chronic_Hep_B_1yr = Prev_Chronic_Hep_B_1yr; % 2 x 100 x (num_years_simul + 1)
output.Prev_Comp_Cirr_1yr = Prev_Comp_Cirr_1yr; % 2 x 100 x (num_years_simul + 1)
output.Prev_Decomp_Cirr_1yr = Prev_Decomp_Cirr_1yr; % 2 x 100 x (num_years_simul + 1)
output.Prev_TDF_treat_1yr = Prev_TDF_treat_1yr; % 2 x 100 x (num_years_simul + 1)
output.NumSAg_1yr = NumSAg_1yr; % 2 x 100 x (num_years_simul + 1)
output.NumSAg_chronic_1yr = NumSAg_chronic_1yr; % 2 x 100 x (num_years_simul + 1)
output.yld_1yr = yld_1yr; % 2 x 100 x (num_years_simul + 1)
output.Incid_Deaths_1yr_approx = Incid_Deaths_1yr_approx; % 2 x 100 x (num_years_simul + 1)
output.Prev_Deaths_1yr = Prev_Deaths_1yr; % 2 x 100 x (num_years_simul + 1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TAM: PAP Chunk 4:
t_PAPoutputs_start = start_year;
t_PAPoutputs_end = end_year;

i_PAPoutputs_start = find(Time >= t_PAPoutputs_start, 1);
i_PAPoutputs_end   = find(Time >= t_PAPoutputs_end, 1);


output.PrevEAg = PrevEAg_of_SAg_5yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x 20 x num_years_output
%% Note that this excludes Vertical transmission
output.NewChronicInfectionRate = Incid_chronic_all_5yr_approx_no_VertTrans(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x 20 x num_years_output
%%output.Tot_Pop_1yr = Tot_Pop_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x 100 x num_years_output
output.NewChronicInfectionRate_NeonatesOnly = Incid_babies_chronic_1yr_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.NumDecompCirr = NumDecompCirr(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output

%% HACK - USING (i_PAPoutputs_end-1) INSTEAD OF i_PAPoutputs_end
i_PAPoutputs_validation_end = i_PAPoutputs_end-1;
assert(max(abs(...
    squeeze(sum(sum(Prev_Decomp_Cirr_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_validation_end),1),2)) - ...
    NumDecompCirr(i_PAPoutputs_start:i_PAPoutputs_validation_end)'...
    )) < 1e-9); % squeeze(sum(sum(Prev_Decomp_Cirr_1yr,1),2)) is a num_years_output x 1 matrix
output.NumLiverCancer = NumLiverCancer(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
assert(max(abs(...
    squeeze(sum(sum(Prev_Liver_Cancer_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_validation_end),1),2)) - ...
    NumLiverCancer(i_PAPoutputs_start:i_PAPoutputs_validation_end)'...
    )) < 1e-8); % squeeze(sum(sum(Prev_Liver_Cancer_1yr,1),2)) is a num_years_output x 1 matrix
%%output.NumSAg_1yr = NumSAg_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x 100 x num_years_output
output.NumEAg_chronic_1yr = NumEAg_chronic_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x 100 x num_years_output
output.NumEAg_chronic_acute_1yr = NumEAg_chronic_acute_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x 100 x num_years_output
%%output.yld_1yr = yld_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x 100 x num_years_output
%%output.Incid_Deaths_1yr_approx = Incid_Deaths_1yr_approx(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x 100 x num_years_output
%%output.Prev_Deaths_1yr = Prev_Deaths_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x 100 x num_years_output
output.num_births_toHbEAgWomenHVL_1yr_approx = num_births_toHbEAgWomenHVL_1yr_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.num_births_toHbEAgWomenLVL_1yr_approx = num_births_toHbEAgWomenLVL_1yr_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.num_births_toHbSAgWomenHVL_1yr_approx = num_births_toHbSAgWomenHVL_1yr_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.num_births_toHbSAgWomenLVL_1yr_approx = num_births_toHbSAgWomenLVL_1yr_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.num_births_1yr_approx = num_births_1yr_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.num_births_chronic_HbEAgWomenHVL_1yr_approx = num_births_chronic_HbEAgWomenHVL_1yr_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.num_births_chronic_HbEAgWomenLVL_1yr_approx = num_births_chronic_HbEAgWomenLVL_1yr_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.num_births_chronic_HbSAgWomenHVL_1yr_approx = num_births_chronic_HbSAgWomenHVL_1yr_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.num_births_chronic_HbSAgWomenLVL_1yr_approx = num_births_chronic_HbSAgWomenLVL_1yr_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.RateBirthDoseVacc = RateBirthDoseVacc(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.RateInfantVacc = RateInfantVacc(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.PeripartumTreatment_HbEAg_HighVL_approx = PeripartumTreatment_HbEAg_HighVL_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.PeripartumTreatment_HbEAg_LowVL_approx = PeripartumTreatment_HbEAg_LowVL_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.PeripartumTreatment_HbSAg_HighVL_approx = PeripartumTreatment_HbSAg_HighVL_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.PeripartumTreatment_HbSAg_LowVL_approx = PeripartumTreatment_HbSAg_LowVL_approx(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.RatePeripartumTreatment = RatePeripartumTreatment(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output
output.PregnantWomenNeedToScreen = PregnantWomenNeedToScreen(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output; added 13.9.15
output.HBVPregnantWomenNeedToEvaluate = HBVPregnantWomenNeedToEvaluate(i_PAPoutputs_start:i_PAPoutputs_end); % 1 x num_years_output

output.beta_U5 = beta_U5;
output.p_HbSAg_av = p_HbSAg_av;

output.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL = p_VertTrans_HbSAgHighVL_NoIntv / p_VertTrans_HbSAgLowVL_NoIntv;
output.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL = p_VertTrans_HbEAgHighVL_NoIntv / p_VertTrans_HbEAgLowVL_NoIntv;

output.p_VerticalTransmission_HbSAgLowVL_NoIntv = p_VertTrans_HbSAgLowVL_NoIntv;
output.p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc = p_VertTrans_HbSAgLowVL_BD;
output.p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP = p_VertTrans_HbSAgLowVL_BD_PAP;
output.p_VerticalTransmission_HbSAgLowVL_PAP = p_VertTrans_HbSAgLowVL_PAP;
output.p_VerticalTransmission_HbSAgHighVL_NoIntv = p_VertTrans_HbSAgHighVL_NoIntv;
output.p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc = p_VertTrans_HbSAgHighVL_BD;
output.p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP = p_VertTrans_HbSAgHighVL_BD_PAP;
output.p_VerticalTransmission_HbSAgHighVL_PAP = p_VertTrans_HbSAgHighVL_PAP;

output.p_VerticalTransmission_HbEAgLowVL_NoIntv = p_VertTrans_HbEAgLowVL_NoIntv;
output.p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc = p_VertTrans_HbEAgLowVL_BD;
output.p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP = p_VertTrans_HbEAgLowVL_BD_PAP;
output.p_VerticalTransmission_HbEAgLowVL_PAP = p_VertTrans_HbEAgLowVL_PAP;
output.p_VerticalTransmission_HbEAgHighVL_NoIntv = p_VertTrans_HbEAgHighVL_NoIntv;
output.p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc = p_VertTrans_HbEAgHighVL_BD;
output.p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP = p_VertTrans_HbEAgHighVL_BD_PAP;
output.p_VerticalTransmission_HbEAgHighVL_PAP = p_VertTrans_HbEAgHighVL_PAP;


outputs_nums_cell_array = {...
    'beta_U5',...
    'p_HbSAg_av',...
    'p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL',...
    'p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL',...
    'p_VerticalTransmission_HbSAgLowVL_NoIntv',...
    'p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc',...
    'p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP',...
    'p_VerticalTransmission_HbSAgLowVL_PAP',...
    'p_VerticalTransmission_HbSAgHighVL_NoIntv',...
    'p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc',...
    'p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP',...
    'p_VerticalTransmission_HbSAgHighVL_PAP',...
    'p_VerticalTransmission_HbEAgLowVL_NoIntv',...
    'p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc',...
    'p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP',...
    'p_VerticalTransmission_HbEAgLowVL_PAP',...
    'p_VerticalTransmission_HbEAgHighVL_NoIntv',...
    'p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc',...
    'p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP',...
    'p_VerticalTransmission_HbEAgHighVL_PAP'...
    };
num_outputs_nums = length(outputs_nums_cell_array);
outputs_vectors_cell_array = {...
    'Time',...
    'NewChronicInfectionRate_NeonatesOnly',...
    'NumDecompCirr',...
    'NumLiverCancer',...
    'num_births_toHbEAgWomenHVL_1yr_approx',...
    'num_births_toHbEAgWomenLVL_1yr_approx',...
    'num_births_toHbSAgWomenHVL_1yr_approx',...
    'num_births_toHbSAgWomenLVL_1yr_approx',...
    'num_births_1yr_approx',...
    'num_births_1yr',...
    'num_births_chronic_HbEAgWomenHVL_1yr_approx',...
    'num_births_chronic_HbEAgWomenLVL_1yr_approx',...
    'num_births_chronic_HbSAgWomenHVL_1yr_approx',...
    'num_births_chronic_HbSAgWomenLVL_1yr_approx',...
    'RateBirthDoseVacc',...
    'RateInfantVacc',...
    'PeripartumTreatment_HbEAg_HighVL_approx',...
    'PeripartumTreatment_HbEAg_LowVL_approx',...
    'PeripartumTreatment_HbSAg_HighVL_approx',...
    'PeripartumTreatment_HbSAg_LowVL_approx',...
    'RatePeripartumTreatment',...
    'PregnantWomenNeedToScreen',...
    'HBVPregnantWomenNeedToEvaluate'...
    };
num_outputs_vectors = length(outputs_vectors_cell_array);
outputs_3D_cell_array = {...
    'PrevEAg',...
    'NewChronicInfectionRate',...
    'Tot_Pop_1yr',...
    'NumSAg_1yr',...
    'NumSAg_chronic_1yr',...
    'NumEAg_chronic_1yr',...
    'NumEAg_chronic_acute_1yr',...
    'yld_1yr',...
    'Incid_Deaths_1yr_approx',...
    'Incid_chronic_all_1yr_approx',...
    'Prev_Deaths_1yr',...
    'Prev_Immune_Reactive_1yr',...
    'Prev_Chronic_Hep_B_1yr',...
    'Prev_Comp_Cirr_1yr',...
    'Prev_Decomp_Cirr_1yr',...
    'Prev_TDF_treat_1yr',...
    };
num_outputs_3D = length(outputs_3D_cell_array);

assert(all(ismember(fields(output),[outputs_nums_cell_array,outputs_vectors_cell_array,outputs_3D_cell_array])))
assert(all(ismember([outputs_nums_cell_array,outputs_vectors_cell_array,outputs_3D_cell_array],fields(output)))) % cell arrays contain the same elements
assert(isequal(sort(fields(output)),sort([outputs_nums_cell_array,outputs_vectors_cell_array,outputs_3D_cell_array]')))

for ii=1:num_outputs_nums
    fieldname = outputs_nums_cell_array{ii};
    field_output = output.(fieldname);
    assert(isscalar(field_output))
end


n_years_PAPoutputs = t_PAPoutputs_end - t_PAPoutputs_start + 1;

for ii=1:num_outputs_vectors
    fieldname = outputs_vectors_cell_array{ii};
    field_output = output.(fieldname);
    assert(length(field_output)==n_years_PAPoutputs)
end

for ii=1:num_outputs_3D
    fieldname = outputs_3D_cell_array{ii};
    field_output = output.(fieldname);
    assert(size(field_output,3)==n_years_PAPoutputs)
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MP: comment this out for now as I don't know what it's doing.
% Diagnositic work-up required for this simulation:
% output.Dx_At_ANC_HBsAG = 0;
% output.Dx_At_ANC_HBeAG = 0;
% output.Dx_At_ANC_VL = 0; 
% 
% % Consider PAP for those who get BD:
% if  (0==cov_BirthDoseAndTDF_EAgHighVL) && (0==cov_BirthDoseAndTDF_SAgHighVL) && ...
%     (0==cov_BirthDoseAndTDF_EAgLowVL) && (0==cov_BirthDoseAndTDF_SAgLowVL) 
% 
%         % Zero coverage of all types of PAP 
%         output.Dx_At_ANC_HBsAG = 0;
%         output.Dx_At_ANC_HBeAG = 0;
%         output.Dx_At_ANC_VL = 0; 
% 
% elseif (   (cov_BirthDoseAndTDF_EAgHighVL==cov_BirthDoseAndTDF_SAgHighVL) && (cov_BirthDoseAndTDF_SAgHighVL==cov_BirthDoseAndTDF_EAgLowVL) ...
%                 && (cov_BirthDoseAndTDF_EAgLowVL==cov_BirthDoseAndTDF_SAgLowVL) && (cov_BirthDoseAndTDF_SAgLowVL==cov_BirthDoseAndTDF_SAgHighVL)  )
%         % PAP not differentiated by E or VL, so just use HBSAG
%         output.Dx_At_ANC_HBsAG = 1;
%         output.Dx_At_ANC_HBeAG = 0;
%         output.Dx_At_ANC_VL = 0; 
% 
% elseif (   (cov_BirthDoseAndTDF_EAgHighVL==cov_BirthDoseAndTDF_EAgLowVL) && (cov_BirthDoseAndTDF_SAgHighVL==cov_BirthDoseAndTDF_SAgLowVL) ...
%                 && (cov_BirthDoseAndTDF_SAgLowVL~=cov_BirthDoseAndTDF_EAgLowVL) && (cov_BirthDoseAndTDF_SAgHighVL~=cov_BirthDoseAndTDF_EAgHighVL)    )    
%         % PAP is differentiated only by E/S
%         output.Dx_At_ANC_HBsAG = 1;
%         output.Dx_At_ANC_HBeAG = 1;
%         output.Dx_At_ANC_VL = 0;  
% 
% 
% elseif (   (cov_BirthDoseAndTDF_SAgHighVL==cov_BirthDoseAndTDF_EAgHighVL) && (cov_BirthDoseAndTDF_SAgLowVL==cov_BirthDoseAndTDF_EAgLowVL) ...
%                 && (cov_BirthDoseAndTDF_SAgLowVL~=cov_BirthDoseAndTDF_SAgHighVL) && (cov_BirthDoseAndTDF_EAgLowVL~=cov_BirthDoseAndTDF_EAgHighVL) )
%          % PAP is differentiated only by VL   
%         output.Dx_At_ANC_HBsAG = 1;
%         output.Dx_At_ANC_HBeAG = 0;
%         output.Dx_At_ANC_VL = 1;         
% 
% elseif (   (cov_BirthDoseAndTDF_EAgHighVL~=cov_BirthDoseAndTDF_SAgHighVL) && (cov_BirthDoseAndTDF_SAgHighVL~=cov_BirthDoseAndTDF_EAgLowVL) ...
%                 && (cov_BirthDoseAndTDF_EAgLowVL~=cov_BirthDoseAndTDF_SAgLowVL) && (cov_BirthDoseAndTDF_SAgLowVL~=cov_BirthDoseAndTDF_SAgHighVL)  )
%         % PAP is differentiated by both E and VL
%         output.Dx_At_ANC_HBsAG = 1;
%         output.Dx_At_ANC_HBeAG = 1;
%         output.Dx_At_ANC_VL = 1;    
% 
% else
%         % Fail
%         assert(false)
% 
% end
% 
% 
% % Consider those who do not get BD
% if (0==cov_TDFOnly_EAgHighVL) && (0==cov_TDFOnly_SAgHighVL) && (0==cov_TDFOnly_EAgLowVL) && (0==cov_TDFOnly_SAgLowVL) 
%         % Zero coverage so do nothing
% 
% elseif (   (cov_TDFOnly_EAgHighVL==cov_TDFOnly_SAgHighVL) && (cov_TDFOnly_SAgHighVL==cov_TDFOnly_EAgLowVL) ...
%                 && (cov_TDFOnly_EAgLowVL==cov_TDFOnly_SAgLowVL) && (cov_TDFOnly_SAgLowVL==cov_TDFOnly_EAgHighVL) )
%         % PAP not differentiated by E or VL, so just use HBSAG
%         output.Dx_At_ANC_HBsAG = min(1,output.Dx_At_ANC_HBsAG+1);
% 
% elseif (   (cov_TDFOnly_EAgHighVL==cov_TDFOnly_EAgLowVL) && (cov_TDFOnly_SAgHighVL==cov_TDFOnly_SAgLowVL) ...
%                 && (cov_TDFOnly_SAgLowVL~=cov_TDFOnly_EAgLowVL) && (cov_TDFOnly_SAgHighVL~=cov_TDFOnly_EAgHighVL)     )
%         % PAP is differentiated only by E/S
%         output.Dx_At_ANC_HBsAG = min(1,output.Dx_At_ANC_HBsAG+1);
%         output.Dx_At_ANC_HBeAG = min(1,output.Dx_At_ANC_HBeAG+1);
% 
% elseif (   (cov_TDFOnly_SAgHighVL==cov_TDFOnly_SAgHighVL) && (cov_TDFOnly_SAgLowVL==cov_TDFOnly_SAgLowVL) ...
%                 && (cov_TDFOnly_SAgLowVL~=cov_TDFOnly_SAgHighVL) && (cov_TDFOnly_EAgLowVL~=cov_TDFOnly_EAgHighVL)  )
%         % PAP is differentiated only by VL   
%         output.Dx_At_ANC_HBsAG = min(1,output.Dx_At_ANC_HBsAG+1);
%         output.Dx_At_ANC_VL = min(1,output.Dx_At_ANC_VL+1);        
% 
% elseif (   (cov_TDFOnly_EAgHighVL==cov_TDFOnly_SAgHighVL) && (cov_TDFOnly_SAgHighVL==cov_TDFOnly_EAgLowVL) ...
%                 && (cov_TDFOnly_EAgLowVL==cov_TDFOnly_SAgLowVL) && (cov_TDFOnly_SAgLowVL==cov_TDFOnly_SAgHighVL) )
%         % PAP is differentiated by both E and VL
%         output.Dx_At_ANC_HBsAG = min(1,output.Dx_At_ANC_HBsAG+1);
%         output.Dx_At_ANC_HBeAG = min(1,output.Dx_At_ANC_HBeAG+1);
%         output.Dx_At_ANC_VL = min(1,output.Dx_At_ANC_VL+1);    
% 
% end
% 
% assert(isscalar(output.Dx_At_ANC_HBsAG))
% assert(isscalar(output.Dx_At_ANC_HBeAG))
% assert(isscalar(output.Dx_At_ANC_VL))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TAM: End of PAP Chunk 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DALYs = make_daly_mat(output,num_years_simul,num_year_1980_2100,life_expectancy);

size(sum(DALYs,1))
i_1950 = find(Time >= 1950, 1);
size(results_to_print(:,i_1950:end))

DALYs_summed = sum(DALYs,1);
disp("Sizes")
size(DALYs_summed)

i_1950 = find(Time >= 1950, 1);
size(results_to_print(:,i_1950:end))
disp("End")
if(store_results_as_text==1)
    if(stochas_run_str=="1")
        %disp("Making header.txt file")
        output_header = construct_header(agegroups_5yr, num_disease_states, num_sexes, num_treat_blocks);
        writelines(output_header, fullfile(basedir,'outputs',"header.txt"))
    end
    i_1950 = find(Time >= 1950, 1);
    filename_results_csv = strcat('results_',ISO,'_scenario',string(scenario_num),'_',sensitivity_analysis,'_run_', stochas_run_str, '.csv');
    disp(fullfile(basedir,'outputs',filename_results_csv))
    %%%size(results_to_print(:,i_1950:end))
    %%%size(Time(i_1950:end))
    writematrix([Time(i_1950:end);results_to_print(:,i_1950:end);DALYs_summed(i_1950:end)]',fullfile(basedir,'outputs',filename_results_csv));
end


end % end function HBVmodel_PPT

function output_labels=construct_header(agegroups, num_disease_states, num_sexes, num_treat_blocks)
    assert(num_sexes==2)
    assert(num_treat_blocks==2)

    % Create labels for disease stage:
    D_labels = strings(1, num_disease_states); for i = 1:num_disease_states; D_labels(i) = "D" + string(i); end
    
    n_age_groups = max(agegroups);
    age_width = 100/max(agegroups);
    age_labels = strings(1, n_age_groups); 
    for i = 1:n_age_groups
        age_min = string((i-1)*age_width);
        age_max = string(i*age_width-1);
        age_labels(i) = "Age" + age_min + "_"+age_max;
    end
    sex_labels = ["F","M"]; % F first in this model
    treat_labels = ["NoTreat","Treat"];

    %%output_labels = strings(num_disease_states,n_age_groups,num_sexes,num_treat_blocks);
    output_labels = "Year,";
    for t=1:num_treat_blocks
        for k=1:num_sexes
            for d=1:num_disease_states
                for a=1:n_age_groups
                    %%output_labels(d,a,k,t) = age_labels(a) + sex_labels(k) + "_" + D_labels(d) + treat_labels(t); 
                    temp_label = age_labels(a) + sex_labels(k) + "_" + D_labels(d) + treat_labels(t) +","; 
                    output_labels = output_labels+temp_label;
                end
            end
        end
    end
    %%output_labels = reshape(output_labels, [1,num_disease_states*n_age_groups*num_sexes*num_treat_blocks]);

    %% Incidence outputs:
    output_labels = output_labels + "Incidence_neonatal,";
    for a=1:n_age_groups
        temp_label = "Incid"+age_labels(a) +","; 
        output_labels = output_labels+temp_label;
    end

    %% Death outputs:
    for a=1:n_age_groups
        temp_label = "Death"+age_labels(a) +","; 
        output_labels = output_labels+temp_label;
    end

    %% Resources (for costing):
    output_labels = output_labels + "NBirthDose,NBD_MAP,NBD_CPAD,N_InfantVacc,N_PAP_EAgHVL,N_PAP_EAgLVL,N_PAP_SAgHVL,N_PAP_SAgLVL,N_screen_PAP,N_starting_treatment,DALYs";

               

end % End function output_labels
