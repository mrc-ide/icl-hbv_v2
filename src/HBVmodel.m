function output = HBVmodel(source_HBsAg,...
    num_year_divisions,dt,ages,num_age_steps,i_natural_hist,i_sexes,i_care, ...
    start_year,num_years_simul,...
    theta,ECofactor, ...
    treatment_rate_params, treat_start_year, ...
    params, PAP_VL_params, PAP_cov_params, ...
    Global_intervention_params, Intervention_data_thiscountry, ...
    p_ChronicCarriage,Prog,Transitions, ...
    scenario_BDcoverage, scenario_BDcoverage_fromMAP, ...
    scenario_BDcoverage_fromCPAD, scenario_HepB3coverage, ...
    scenario_Treatment, I_TREAT, scenario_treat_elig, scenario_data_ANCHBVtestingbyage_thiscountry, ...
    ISO, scenario_num, scenario_AddScreenIntervention, ...
    num_year_1980_2100, life_expectancy, ...
    stochas_run_str, sensitivity_analysis, basedir, store_results_as_text)


DUMMY_VALUE = -99;  % Used in initialising arrays to a dummy value (-99 should be easy to spot).

%% Establish basic simulation parameters
agegroups_5yr = 1 + floor(ages / 5); % categorises the ages into age-groups of 5 year width; 1 x 1000 double; [1 1 ... 20 20], each number present 50 times
agegroups_1yr = 1 + floor(ages); % categorises the ages into age-groups of 1 year width; 1 x 1000 double; [1 1 ... 100 100], each number present 10 times

%% This is the number of 1 year age groups (100):
num_1yr_age_gps = max(agegroups_1yr);

%% markers for key age boundaries
i6mo = find(ages >= 0.5, 1); 
i1y = find(ages >= 1, 1);
i5y = find(ages >= 5, 1);
i15y = find(ages >= 15, 1);
i30y = find(ages >= 30, 1); % age boundary for different treatment eligibility

end_year = start_year + num_years_simul; % 2101
TimeSteps = start_year:dt:end_year; % 1 x 2101 double; [1890 1890.1 1890.2 ... 2099.8 2099.9 2100 2100.1 ... 2100.8 2100.9 2101]


% X-stocks are (infection_state, age, sex(1=women, 2=men), accessible*)   {*accessible
% specifies whether this person can be reached by treatment progs, 1=no, 2=yes}
% So X() is 4D array of size num_disease_states *  num_age_steps * num_sexes * num_treat_blocks

num_sexes = i_sexes.n_sexes; % F=1, M=2.
i_female = i_sexes.F;
i_male = i_sexes.M;

% values for indices for 4th index (stratification as to whether would get tested/treated if needed)
% Code relating to this (providing a cap for treatment) is marked with LECZENIE.
%% ALPHA-1
num_treat_blocks = i_care.n_care_blocks; 
%%i_notseektreat = 1;
%%i_seektreat = 2;     % would get tested/treated
i_undiagnosed = i_care.undiagnosed;  % Undiagnosed (or never infected)
i_appropriate_management = i_care.appropriate_management;     % following diagnosis, gets appropriate management (monitoring if ineligible, on treatment and adherent if eligible)
i_incare_nonadherent = i_care.incare_nonadherent;         % following diagnosis, has sub-optimal management or goes on treatment but is non-adherent.
i_outofcare = i_care.outofcare;                  % following diagnosis, drops out of care. Would need new technology/guidelines to improve outcome to another compartment (e.g. treat-all

num_disease_states = i_natural_hist.n_nathist_states;
i_Susc = i_natural_hist.Susc;               % 'Susceptible', 
i_ImmTol = i_natural_hist.ImmTol;           % 'HBV: Immune Tolerant' : HBeAg+ with very high HBV DNA (>1e6IU/ml), normal ALT
i_ImmReact = i_natural_hist.ImmReact;       % 'HBV: Immune Reactive' :  HBeAg+ with high HBV DNA (>20000IU/ml), elevated ALT
i_AsymptCarr = i_natural_hist.AsymptCarr;   % 'HBV: Asymptomatic Carrier': HBeAg- low/undetectable HBV DNA, normal ALT
i_Chronic = i_natural_hist.Chronic;         % 'HBV: Chronic Hep B': HBeAg-, moderate to high HBV DNA levels, fluctuating/persistently elevated ALT 
i_CompCirr = i_natural_hist.CompCirr;       % 'HBV: Comp Cirrhosis',
i_DecompCirr = i_natural_hist.DecompCirr;   % 'HBV: Decomp Cirrhosis',
i_HCC = i_natural_hist.HCC;                 % 'HBV: Liver Cancer',
i_Immune = i_natural_hist.Immune;           % 'HBV: Immune (Rec. or vacc.)',
%% ALPHA - remove i_TDFtreat
i_TDFtreat = i_natural_hist.TDFtreat_LEGACY;    % 'HBV: TDF-Treatment',
i_HBVdeath = i_natural_hist.HBVdeath;       % 'Prematurely dead due to HBV', ... % 11
i_3TCtreat = i_natural_hist.i3TCtreat_LEGACY;   % '3TC-Treatment', ... % 12
i_3TCfailed = i_natural_hist.i3TCfailed_LEGACY; % 'Failed 3TC-Treatment', ...  % 13
i_NonSevAcute = i_natural_hist.NonSevAcute; % 'Non-severe acute', ...  % 14
i_SevereAcute = i_natural_hist.SevereAcute; % 'Severe acute' ...  % 15
%% ALPHA - update these counters to make sure they don't include i_TDFtreat = 10.
i_alive = setdiff(1:num_disease_states, i_HBVdeath);
i_acute = [i_NonSevAcute, i_SevereAcute];
i_eAgpos_chronic = [i_ImmTol, i_ImmReact];  %% Immune Tolerant, Immune Reactive, Non-severe + severe acute.

i_eAgpos = sort([i_eAgpos_chronic, i_acute]);  %% Immune Tolerant, Immune Reactive, Non-severe + severe acute.

%% Asymptomatic carrier, Chronic, Comp+Decom Cirr, HCC, failed 3TC. Represents the infectious (but less infectious than eAg+) stages
i_sAgpos_notEagpos_notreat = [i_AsymptCarr i_Chronic i_CompCirr i_DecompCirr i_HCC i_3TCfailed];

%% ALPHA - code no longer used
%%i_treateligible = [i_ImmReact i_Chronic i_CompCirr i_DecompCirr]; %% Immune Reactive, Chronic, Comp+Decomp Cirr
i_treat_legacy = [i_TDFtreat i_3TCtreat];  %% These are legacy states (likely to be repurposed for e.g. cure)

%% Includes 10 (TDFtreat) and 12 (3TCtreat) states
i_sAgpos_chronic = sort([i_sAgpos_notEagpos_notreat, i_eAgpos_chronic, i_treat_legacy]);
i_sAgpos = sort([i_sAgpos_chronic, i_acute]);

%i_sAgpos = [2:8 10 12:15];   %% Includes 10 (TDFtreat) and 12 (3TCtreat) states
%i_sAgpos_chronic = [2:8 10 12:13]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here we sort out the natural history states that are eligble for treatment in the current scenario:
% i_eAgpos_treatelig_under30 = i_ImmReact;  %% Under 30, immune reactive is eligible but not immune tolerant.
% i_eAgpos_treat_inelig_under30 = [i_ImmTol, i_acute];
% i_eAgpos_treatelig_30plus  = [i_ImmTol i_ImmReact];  %% Only acute is ineligible when aged 30+.
% i_eAgpos_treat_inelig_30plus = i_acute; 
%% sAg+ (not eAg+) states: 

%% Make sure these line up with the eligibility in get_treatment_eligible_ageindices.m:
if(strcmp(scenario_treat_elig,"Current treatment"))
    i_sAgpos_not_eAgpos_treatelig = [i_Chronic i_CompCirr i_DecompCirr];
    i_sAgpos_not_eAgpos_treat_inelig = [i_AsymptCarr i_HCC]; % Just asymptomatic (4) and HCC (8). Exclude the TDF treatment, 3TC treatment/failed treatment compartments.
    i_eAgpos_treatelig_under30 = i_ImmReact;  %% Under 30, immune reactive is eligible but not immune tolerant.
    i_eAgpos_treat_inelig_under30 = [i_ImmTol, i_acute];
    i_eAgpos_treatelig_30plus  = [i_ImmTol i_ImmReact];  %% Only acute is ineligible when aged 30+.
    i_eAgpos_treat_inelig_30plus = i_acute; 
elseif(strcmp(scenario_treat_elig,"Universal treatment"))
    i_sAgpos_not_eAgpos_treatelig = [i_Chronic i_CompCirr i_DecompCirr i_AsymptCarr];
    i_sAgpos_not_eAgpos_treat_inelig = i_HCC; % Just  HCC (8). 
    i_eAgpos_treatelig_under30 = [i_ImmTol, i_ImmReact];  %% Under 30, immune reactive is eligible but not immune tolerant.
    i_eAgpos_treat_inelig_under30 = i_acute;
    i_eAgpos_treatelig_30plus  = [i_ImmTol i_ImmReact];  %% Only acute is ineligible when aged 30+.
    i_eAgpos_treat_inelig_30plus = i_acute; 
else
    disp("Error - unknown value for scenario_treat_elig in HBVmodel.m. Exiting")
    return
end
i_treatelig_under30 = sort([i_sAgpos_not_eAgpos_treatelig i_eAgpos_treatelig_under30]);
i_treat_inelig_under30 = sort([i_sAgpos_not_eAgpos_treat_inelig i_eAgpos_treat_inelig_under30]);
i_treatelig_30plus = sort([i_sAgpos_not_eAgpos_treatelig i_eAgpos_treatelig_30plus]);
i_treat_inelig_30plus = sort([i_sAgpos_not_eAgpos_treat_inelig i_eAgpos_treat_inelig_30plus]);
 
    
    
    
    
    


%% Now check the above are consistent with the eligibility criteria in get_treatment_eligible_ageindices() 
%% - that is used to modify the natural history progression when on treatment

temp_ImmTol = get_treatment_eligible_ageindices(scenario_treat_elig, i_ImmTol, i_natural_hist, ages);
temp_ImmReact = get_treatment_eligible_ageindices(scenario_treat_elig, i_ImmReact, i_natural_hist, ages);
temp_AsymptCarr = get_treatment_eligible_ageindices(scenario_treat_elig, i_AsymptCarr, i_natural_hist, ages);
temp_Chronic = get_treatment_eligible_ageindices(scenario_treat_elig, i_Chronic, i_natural_hist, ages);
temp_CompCirr = get_treatment_eligible_ageindices(scenario_treat_elig, i_CompCirr, i_natural_hist, ages);
temp_DecompCirr = get_treatment_eligible_ageindices(scenario_treat_elig, i_DecompCirr, i_natural_hist, ages);
temp_HCC = get_treatment_eligible_ageindices(scenario_treat_elig, i_HCC, i_natural_hist, ages);
if(strcmp(scenario_treat_elig,"Current treatment"))
    if (~isequal(temp_ImmTol,i30y:num_age_steps) || ~isequal(temp_ImmReact,1:num_age_steps) || ...
        ~isempty(temp_AsymptCarr) || ...
        ~isequal(temp_Chronic,1:num_age_steps) || ~isequal(temp_CompCirr,1:num_age_steps) || ...
        ~isequal(temp_DecompCirr,1:num_age_steps) || ~isempty(temp_HCC))
        disp("Error - one or more of the treatment eligibility criteria (Current treatment) don't match. Exiting")
        return
    end
elseif(strcmp(scenario_treat_elig,"Universal treatment"))
    if (~isequal(temp_ImmTol,1:num_age_steps) || ~isequal(temp_ImmReact,1:num_age_steps) || ...
        ~isequal(temp_AsymptCarr,1:num_age_steps) || ...
        ~isequal(temp_Chronic,1:num_age_steps) || ~isequal(temp_CompCirr,1:num_age_steps) || ...
        ~isequal(temp_DecompCirr,1:num_age_steps) || ~isempty(temp_HCC))
        disp("Error - one or more of the treatment eligibility criteria (Universal treatment) don't match. Exiting")
        return
    end
else
    disp("Error - unknown value for scenario_treat_elig in HBVmodel.m when cross-checking eligibility criteria. Exiting")
    return
end




% else %% Current treatment:
%     RRtrans_effective_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRtrans_effective_TDFtreatment'),:).Value;
%     RRtrans_nonadherent_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRtrans_nonadherent_TDFtreatment'),:).Value;
% end


%% Note that params.dwvec is an external vector that (*SHOULD*) mimic the natural history states of the model. So we need to check that here:
%% There are three non-zero DALY weights (for *alive* states - the DALYs from death are dealt with separately in make_daly_mat.m).
%% **NEVER CHANGE THE 7,8,11,15 BELOW TO VARIABLES - THEY HAVE TO BE (MAGIC) NUMBERS***
assert(i_DecompCirr==7 && i_HCC==8 && i_SevereAcute==15)
assert(i_HBVdeath==11)
if ~isequal(find(params.dwvec>0),[7,8,15])
    disp("Error - params.dwvec has been modified so that non-zero indices no longer correspond to 7,8,15. Exiting")
    return
end



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
%%%p_VerticalTransmission_HbEAg_BD = p_VerticalTransmission_HbEAg_NoBD * (1 - params.Efficacy_BirthDoseVacc_HbEAg); 
%%%p_VerticalTransmission_HbEAg_BirthDose_MAP_CPAD = p_VerticalTransmission_HbEAg_NoBD * (1 - efficacy_MAP_CPAD_HbEAg); 
%%%p_VerticalTransmission_Tr_NoBD = p_VerticalTransmission_HbEAg_NoBD * (1 - Efficacy_Treatment_MTCT); % probability of transmission from an HBeAg+ mother on treatment to her baby without intervention


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PAP chunk 1:
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
    %% Now replaced with e.g. PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL
    % % Coverage of PAP among those with BD
    % [cov_BirthDoseAndTDF_EAgHighVL_itt, ...
    %     cov_BirthDoseAndTDF_EAgLowVL_itt, ...
    %     cov_BirthDoseAndTDF_SAgHighVL_itt, ...
    %     cov_BirthDoseAndTDF_SAgLowVL_itt] = ...
    %     deal(zeros(size(TimeSteps)));

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
    RateOfPAPInitiation, HBVPositivePregnantWomenAtANC,...
    num_starting_treatment_as_eligible_this_year] = deal(0);

  %% For now let's just count the number of people starting treatment (unstratified by age/sex):
  number_starting_treatment_to_print = 0;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of PAP chunk 1
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
%% However it is currently set (by the parameters) to be unused legacy code.
CalendarTimeReducInTransmission = params.ReducInTransmission;              % Fractional reduction in transmission. Currently set to 0
YearCalendarReducInTransmission = params.YearReducInTransmission;      % Turning point year for reduction. Currently set to 2100
DurCalendarTimeReducInTransmission = 15;                                  % Time taken to complete change - MP: note that this is not actually that. 
%% MP: beta_scaler seems to be legacy code. Currently CalendarTimeReducInTransmission is 0. Otherwise (even with YearCalendarReducInTransmission=2100)
%% we still get some reduction in beta, and quite a large reduction after 2080 (reaching 50% reduction in 2100).

%% Check the number of timesteps per year (num_year_divisions) is an integer (and not infinity/non-integer/NaN):
assert(mod(num_year_divisions, 1)==0)
beta_scaler = CalendarTimeReducInTransmission ./ (1 + exp( (TimeSteps - (YearCalendarReducInTransmission)) ./ (DurCalendarTimeReducInTransmission / num_year_divisions) ));

beta_U5_SAg = beta_U5 * (1 - CalendarTimeReducInTransmission) + beta_U5 * zeros(size(beta_scaler));  % NOT TIME DEPENDENT as no beta_scaler.
beta_U5_EAg = min(1.0, beta_U5_SAg * ECofactor); % probabilities have to be capped at 1.

beta_1to15_SAg = beta_1to15 * (1 - CalendarTimeReducInTransmission) + beta_1to15 * beta_scaler;  % TIME DEPENDENT as beta_scaler is.
beta_1to15_EAg = min(1.0, beta_1to15_SAg * ECofactor); % probabilities have to be capped at 1.

beta_5plus_SAg = beta_5plus * (1 - CalendarTimeReducInTransmission) + beta_5plus * beta_scaler;  % TIME DEPENDENT as beta_scaler is.
beta_5plus_EAg = min(1.0, beta_5plus_SAg * ECofactor); % probabilities have to be capped at 1.



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
    % params.HBsAg_prevs_middle_year_1 is a 18 x 2 double of age group by gender
    % params.HBsAg_prevs_middle_year_1 age groups: 0--4 5--9 10--14 15--19 20--24 25--29 30--34 35--39 40--44 45--49 50--54 55--59 60--64 65--69 70--74 75--79 80--84 85+
    % So the "3" is because we replicate the 85+ prevalence (so in the
    % model it is the prevalence in 85-89, 90-94, 95-99)
    StartPrev_byAgeGroups = [StartPrev_byAgeGroups(1:end-1,:); repmat(StartPrev_byAgeGroups(end,:),3,1)];
    %% Magic number 20 - probably length(unique(agegroups_5yr))
    assert(isequal(size(StartPrev_byAgeGroups),[20 num_sexes]))

    NumSAg = StartPrev_byAgeGroups(agegroups_5yr, :) .* StartPop;
    % agegroups_5yr is a 1 x 1000 double; [1 1 ... 20 20], each number present 50 times
    % expanding StartPrev_byAgeGroups from 5 year age steps to 0.1 year age steps
    NumNotSAg = (1 - StartPrev_byAgeGroups(agegroups_5yr, :)) .* StartPop;
elseif strcmp(source_HBsAg,'CDA')
    %% MP: Magic numbers 99.9, 6, 1, 2
    StartPrev_byAgeGroups = [repmat(params.country_HBsAg_prevalences_by_ages_mid_1_young_old(1),num_year_divisions*(5.9-0.0)+1,2); ...
        repmat(params.country_HBsAg_prevalences_by_ages_mid_1_young_old(2),num_year_divisions*(99.9-6.0)+1, num_sexes)];
    % apply prevalence in 5-year-olds to 0 to 6 year olds; apply prevalence in all ages to 6 to 99 year olds
    assert(isequal(size(StartPrev_byAgeGroups),[num_age_steps num_sexes]))

    NumSAg = StartPrev_byAgeGroups .* StartPop;
    NumNotSAg = (1 - StartPrev_byAgeGroups) .* StartPop;
elseif strcmp(source_HBsAg,'WHO')
    under_5_pos_vec_len = length(find(ages<=5.0));
    over_5_pos_vec_len = length(find(ages>5.0));
    assert(under_5_pos_vec_len+over_5_pos_vec_len==num_age_steps)
    %% MP: Magic numbers 1, 2 (these correspond to the U5 prevalence and the >5 prevalences in the data).
    StartPrev_byAgeGroups = [ ...
        repmat(params.country_HBsAg_prevalences_by_ages_prevacc_young_old(1),under_5_pos_vec_len,num_sexes); ...
        repmat(params.country_HBsAg_prevalences_by_ages_prevacc_young_old(2),over_5_pos_vec_len,num_sexes) ...
        ];
    assert(isequal(size(StartPrev_byAgeGroups),[num_age_steps num_sexes]))

    NumSAg = StartPrev_byAgeGroups .* StartPop;
    NumNotSAg = (1 - StartPrev_byAgeGroups) .* StartPop;
end    

%% LECZENIE:

%% This is the initial distribution of people:
X(i_Susc, :, :, i_undiagnosed) = NumNotSAg;
X(i_ImmReact, :, :, i_undiagnosed) = 0.5 * NumSAg;
X(i_AsymptCarr, :, :, i_undiagnosed) = 0.5 * NumSAg;

%% Previously incorporated the stratification i_notseektreat/i_seektreat which determines if someone will not/will seek treatment (currently TDF) for chronic HBV in future.
%X(i_Susc, :, :, i_undiagnosed) = NumNotSAg * (1-(treatment_rate_params.prop_diagnosed_t0*treatment_rate_params.prop_treatifdiag_t0));
%X(i_Susc, :, :, i_seektreat) = NumNotSAg * (treatment_rate_params.prop_diagnosed_t0*treatment_rate_params.prop_treatifdiag_t0);

%% MP: Magic number 0.5 (and 1-0.5) - put in main_script.m
%X(i_ImmReact, :, :, i_undiagnosed) = 0.5 * NumSAg * (1-(treatment_rate_params.prop_diagnosed_t0*treatment_rate_params.prop_treatifdiag_t0));
%X(i_ImmReact, :, :, i_seektreat) = 0.5 * NumSAg * (treatment_rate_params.prop_diagnosed_t0*treatment_rate_params.prop_treatifdiag_t0);
%X(i_AsymptCarr, :, :, i_undiagnosed) = 0.5 * NumSAg * (1-(treatment_rate_params.prop_diagnosed_t0*treatment_rate_params.prop_treatifdiag_t0));
%X(i_AsymptCarr, :, :, i_seektreat) = 0.5 * NumSAg * (treatment_rate_params.prop_diagnosed_t0*treatment_rate_params.prop_treatifdiag_t0);


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

[NumSAg_5yr, PrevEAg_of_SAg_5yr] = deal(-99 * ones(2, max(agegroups_5yr), (num_years_simul+1))); 

%% Note that this excludes Vertical transmission
Incid_chronic_all_5yr_approx_no_VertTrans = zeros(2, max(agegroups_5yr), num_years_simul+1);
%% end of mini-chunk

[...
    Tot_Pop_1yr, Prev_treatment_eligible_1yr, ...
    Prev_Liver_Cancer_1yr, Prev_Decomp_Cirr_1yr, Prev_TDF_treat_1yr, NumSAg_1yr, NumSAg_chronic_1yr, yld_1yr, Prev_Deaths_1yr...
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
    
    nstates_deaths = length(unique(agegroups_5yr)); %% - deaths per year - 5 yr age groups = 20
    nstates_newcases_chroniccarriage = length(unique(agegroups_5yr)) + 1; %% - new cases of chronic carriage/yr (neonates, plus 5 yr age gps) = 21
    ncol_X_to_print_byage = num_disease_states*  num_sexes* num_treat_blocks;
    %% The final 10 are "resources_to_print" outputs (note that the below has 11 outputs - but DALYs are dealt with outside of "results_to_print" and are appended when writing the output via writematrix()):
    %% NBirthDose,NBD_MAP,NBD_CPAD,N_InfantVacc,N_PAP_EAgHVL,N_PAP_EAgLVL,N_PAP_SAgHVL,N_PAP_SAgLVL,N_screen_PAP,N_starting_treatment,DALYs
    n_resource_cols = 10;
    ncol_results_to_print = max(agegroups_5yr) * ncol_X_to_print_byage + nstates_newcases_chroniccarriage + nstates_deaths + n_resource_cols;
    results_to_print = DUMMY_VALUE * ones(ncol_results_to_print, num_years_simul + 1);
    
    %%X_to_print = DUMMY_VALUE * ones(max(agegroups_5yr)*ncol_X_to_print, num_years_simul + 1);
end

%% Get HepB3:
%% Infection stage susceptible (so x1), Age gp - age 6m (so x1), by sex (so x2), by treatment stratum:
%% ALPHA (DONE) - changed for new treatment structure. Previously was transfer_to_HepB3vacc = zeros(1, 1, num_sexes, 2);
%% transfer_to_HepB3vacc is the proportion of infants getting HepB3 who become immune.

transfer_to_HepB3vacc = zeros(1, 1, num_sexes, num_treat_blocks);


%% MP: Magic number 1 is because arrays NewChronicCarriage, moving_btw_states should have same number of dimensions as X()
%% but the first index is single (because not indexing over natural history states).
%% Note - both of these will be across all treatment blocks (though in practice we're assuming infections won't occur if adherent).
[...
    NewChronicCarriage, moving_btw_states, ...
    ] = deal(zeros(1, num_age_steps, num_sexes, num_treat_blocks));



% single output per year
%% PAP - extra PAP-model-specific outputs included here:
[Time, RateInfantVacc, RateBirthDoseVacc, RatePeripartumTreatment, num_starting_treatment_as_eligible,... 
    num_births_1yr, NumDecompCirr, NumLiverCancer, ...
    PregnantWomenNeedToScreen, HBVPregnantWomenNeedToEvaluate] = deal(DUMMY_VALUE * ones(1, num_years_simul+1));
 
[num_births_toHbEAgWomenHVL_1yr_approx, num_births_toHbEAgWomenLVL_1yr_approx, num_births_toHbSAgWomenHVL_1yr_approx, num_births_toHbSAgWomenLVL_1yr_approx, ... 
    num_births_1yr_approx, ...
    num_births_chronic_HbEAgWomenHVL_1yr_approx, num_births_chronic_HbEAgWomenLVL_1yr_approx, num_births_chronic_HbSAgWomenHVL_1yr_approx, num_births_chronic_HbSAgWomenLVL_1yr_approx, ... 
    Incid_babies_chronic_1yr_approx, ...
    PeripartumTreatment_HbEAg_HighVL_approx, PeripartumTreatment_HbEAg_LowVL_approx, PeripartumTreatment_HbSAg_HighVL_approx, PeripartumTreatment_HbSAg_LowVL_approx...
    ] = deal(DUMMY_VALUE * ones(1, num_years_simul+1));


% ----- Simulation -----

i_dt = 1; % i_dt increase every time i.e. every 0.1 years; goes from 1 to 2101 (length of TimeSteps)
OutputEventNum = 1; % OutputEventNum increase every year; goes from 1 to 212
%% ALPHA-2 - DONE. Renamed from "moving_to_treatment" to "moving_to_diagnosed".
moving_to_diagnosed = zeros(size(X));
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


%%moving_to_diagnosed_by_birthcohort_testing = zeros(size(X));
%% ALPHA-3 - DONE. Renamed "treatment" "diagnosed".
%% The below move people to "diagnosed" - some of those will be out of care, some will start treatment.
moving_to_diagnosed_by_birthcohort_testing_per_timestep = zeros(size(X));
moving_to_diagnosed_by_birthcohort_testing_this_timestep = zeros(size(X));
moving_to_diagnosed_by_community_screening_per_timestep = zeros(size(X));
moving_to_diagnosed_by_community_screening_this_timestep = zeros(size(X));

%% Note that moving_to_diagnosed_by_ANC_testing_this_timestep is defined separately 
%% (as it's an ongoing intervention rather than a fixed-period one).

for time = TimeSteps 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Modify treatment-related parameters if needed (natural history, transmission):
    %% Set effectiveness of treatemnt in reducing transmission (different for long-acting treatment):
    if(scenario_Treatment==I_TREAT.LA)
        
        
        if(time<treatment_rate_params.t_treatment_scaleup_start)
            RRtrans_effective_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRtrans_effective_TDFtreatment'),:).Value;
            RRtrans_nonadherent_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRtrans_nonadherent_TDFtreatment'),:).Value;
        else
            RRtrans_effective_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRtrans_effective_LAtreatment'),:).Value;
            RRtrans_nonadherent_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRtrans_nonadherent_LAtreatment'),:).Value;
        end
    else %% Current treatment:
        RRtrans_effective_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRtrans_effective_TDFtreatment'),:).Value;
        RRtrans_nonadherent_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRtrans_nonadherent_TDFtreatment'),:).Value;
    end
    
    %if(scenario_Treatment==I_TREAT_LA)
    if(scenario_Treatment==I_TREAT.LA)
        if(time<treatment_rate_params.t_treatment_scaleup_start)
            prop_adhere_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'prop_adhere_TDF'),:).Value;
            %%prop_nonadhere_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'prop_nonadhere_TDF'),:).Value;
        else %% LA treatment should change the proportions adhering:
            prop_adhere_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'prop_adhere_LAtreat'),:).Value;
            %%prop_nonadhere_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'prop_nonadhere_LAtreat'),:).Value;
        end
    else %% Current treatment:
        prop_adhere_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'prop_adhere_TDF'),:).Value;
        %%prop_nonadhere_treatment = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'prop_nonadhere_TDF'),:).Value;
    end

    %% Intervention-specific parameters - these are the proportion of people having a test who don't enter care:
    birthcohort_prop_Dx_outofcare = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'BirthCohort_prop_Dx_outofcare'),:).Value;
    communityscreening_prop_Dx_outofcare = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'CommunityScreening_prop_Dx_outofcare'),:).Value;
    ANCtesting_prop_Dx_outofcare = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'ANCtesting_prop_Dx_outofcare'),:).Value;
    
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
    %% MP: Magic number 1 - in the previous step we converted fert to a num_1yr_age_gps (100) vector (the 1 is just because it is treating it as a 100x1 matrix):
    assert(isequal(size(fert), [num_1yr_age_gps 1]))
    %% MP: Magic number 1
    fert = repmat(fert',num_year_divisions,1);
    fert = fert(:);  % Reshape fert into a 1D vector from a matrix
    %% MP: Magic number 1
    assert(isequal(size(fert),[num_age_steps 1]))
    assert(length(params.net_migration)==(num_years_simul+1))
    net_migration = params.net_migration(OutputEventNum); %% Single number per year.
    sex_ratio = params.sex_ratios(OutputEventNum); %% Single number per year.
    assert(isscalar(net_migration))
    net_migration = repmat(net_migration, [num_disease_states num_age_steps num_sexes num_treat_blocks]);

    
    % Compute Outputs once per year
    if rem(time, 1) == 0 % only saves variables in this "for" loop every 10 time steps (or once a year, since dt=0.1)
        Time(OutputEventNum) = time;
 
        
        % Rescale population sizes of each age group and gender
        if (time >= 1950)
            base_year_montagu = 1949; %% So 1950 corresponds to index 1.
            n_years_montagu_rescaling = end_year - base_year_montagu; %% 152 years.
            %% Sum over first (natural history) and 4th (treatment) strata - so this is the current population divided by age and sex:
            ModelPop = squeeze(sum(sum(X(i_alive,:,:,:), 1), 4));
            assert(isequal(size(ModelPop),[num_age_steps num_sexes]))
            % sum over disease state of alive people and treatment; ModelPop is 1000 x 2 i.e. age groups versus gender
            assert(isequal(size(params.total_pop_female),[101 n_years_montagu_rescaling]))
            % params.total_pop_female is a 101 x 152 matrix of n_years_montagu_rescaling years (1950 to 2101 inclusive) for 101 age groups (0--0, 1--1, 2--2,..., 98--98, 99--99, 100--100)
            assert(isequal(size(params.total_pop_male),[101 n_years_montagu_rescaling]))
            col_index = time - base_year_montagu;
            %% MP: magic numbers 1:num_1yr_age_gps represent indexes in params.total_pop for ages 0-99
            MontaguPopFemale = params.total_pop_female(1:num_1yr_age_gps,col_index);
            % only want ages 0--99
            MontaguPopMale = params.total_pop_male(1:num_1yr_age_gps,col_index);
            MontaguPop = [MontaguPopFemale MontaguPopMale];
            assert(isequal(size(MontaguPop),[num_1yr_age_gps num_sexes]))    %% MP: num_1yr_age_gps is number of 1-year age gps 0-99.
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
            % add an extra dimension and duplicate it for each disease
            % state and treatment stratum
            X = X .* pop_scaler;
            if(time==2025)
            %disp([min(ScalerMat),max(ScalerMat)])
                disp("Uncomment the line below to show scalarmat")
                %%disp("Scalarmat here:")
                %%disp(ScalerMat)
            end
            % scale all parts of X, including dead people i.e. State i_HBVdeath=11
        end


        


        for k = 1:num_sexes % genders

            for ag = 1:num_1yr_age_gps % 1:100
 
                if OutputEventNum > 1
                
                    %% MP: Magic numbers: sum over second (age) and 4th (treatment states):
                    state_prev_vec = squeeze(sum(sum(X(:, agegroups_1yr == ag, k, :), 2), 4)); % k is gender
                    assert(isequal(size(state_prev_vec),[num_disease_states 1]))

                    Tot_Pop_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec(i_alive));                
                    
                    Prev_Liver_Cancer_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_HCC);
                    Prev_Decomp_Cirr_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_DecompCirr);
                    
                    %% ALPHA - Prev_TDF_treat_1yr is now the sum over groups in care (note that this will depend on eligibility):
                    if(ag<30)
                        i_treatelig_thisage = i_treatelig_under30;
                    else
                        i_treatelig_thisage = i_treatelig_30plus;
                    end
                    %%i_treatelig_thisage = get_treatment_eligible_nathistindices(scenario_treat_elig, ag, i_natural_hist, ages);
                    Prev_TDF_treat_1yr(k, ag, OutputEventNum-1) = squeeze(sum(sum(sum(X(i_treatelig_thisage, agegroups_1yr == ag, k, [i_appropriate_management,i_incare_nonadherent]), 1), 2), 4));
                    Prev_treatment_eligible_1yr(k, ag, OutputEventNum-1) = squeeze(sum(sum(sum(X(i_treatelig_thisage, agegroups_1yr == ag, k, :), 1), 2), 4));
                    
                    %%Prev_TDF_treat_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_TDFtreat);

                    %% MP TODO: Maybe remove this as the HBV deaths compartment has the same "edge of cliff" thing
                    %% where people aged 99 drop off the model once they turn 100 (so anyone who died age 99 is no longer counted).
                    %% Also HBV deaths is rescaled by ScalerMat, which is a bit funky in the older age groups.
                    %% It might be possible to patch it a bit (say truncate at age 95 to reduce the ScalerMat issue, then cumulatively count deaths
                    %% adding the new 95 yo dead people to an existing cumulative counter of people who are dead who would be 95+ now).
                    Prev_Deaths_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_HBVdeath);
                    NumSAg_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec(i_sAgpos));
                    NumSAg_chronic_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec(i_sAgpos_chronic));

                    %% params.dwvec is a length-15 vector (so contains states). For now this is OK (I've added an assert statement which will hopefully stop me making unwanted changes to the state variable indices). 
                    %% "yld" = "(disability-adjusted life-)years living with disease. Deaths (yll_spread) is calculated in make_daly_mat.m using Prev_Deaths_1yr.
                    yld_1yr(k, ag, OutputEventNum-1) = sum( state_prev_vec .* params.dwvec' );

                    %% MP: Magic numbers: sum over second (age) and 4th (treatment states):
                    %% The first index in NewChronicCarriage() has to be 1 (it does not represent susceptibles!) - NewChronicCarriage ...= deal(zeros(1, num_age_steps, num_sexes, num_treat_blocks));
                    Incid_chronic_all_1yr_approx(k,ag,OutputEventNum-1) = sum(sum(NewChronicCarriage(1, agegroups_1yr == ag, k, :), 2), 4);
                    Incid_Deaths_1yr_approx(k, ag, OutputEventNum-1) = sum(state_prev_vec .* Prog(:, i_HBVdeath));

                    Prev_HCC_1yr(k, ag, OutputEventNum-1) = state_prev_vec(i_HCC);
                    NumEAg_chronic_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec(i_eAgpos_chronic));
                    NumEAg_chronic_acute_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec(i_eAgpos));
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
                    % model results assigned to a particular year at the beginning of that year, after which they are zero'd
                    NumSAg_5yr(k, ag, OutputEventNum-1) = sum(sum(sum(X(i_sAgpos, agegroups_5yr == ag, k, :))));
                
                    if NumSAg_5yr(k, ag, OutputEventNum-1)>0
                        PrevEAg_of_SAg_5yr(k,ag,OutputEventNum-1) = sum(sum(sum(X(i_eAgpos_chronic, agegroups_5yr == ag, k, :)))) / NumSAg_5yr(k, ag, OutputEventNum-1);
                        % Note that this is prevalence of e+ among s+
                    else
                        disp(NumSAg_5yr(k, ag, OutputEventNum-1))
                        assert(NumSAg_5yr(k, ag, OutputEventNum-1)==0)
                        assert(sum(sum(sum(X(i_eAgpos_chronic, agegroups_5yr == ag, k, :))))==0)
                        PrevEAg_of_SAg_5yr(k,ag,OutputEventNum-1) = 0;
                    end
                    %% The first index in NewChronicCarriage() has to be 1 (it does not represent susceptibles!) - NewChronicCarriage ...= deal(zeros(1, num_age_steps, num_sexes, num_treat_blocks));
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
        %% PAP mini-chunk 2
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
            %% transfer_to_HepB3vacc is a 4d array, but the first two dimensions are trivial (as only susceptible infants 6m get HepB3 vaccination). So we sum over sex and treatment cascade status (though in the current model infants are all in the "undiag" stratum).
            RateInfantVacc(OutputEventNum-1) = squeeze(sum(sum(transfer_to_HepB3vacc,3),4)) * num_year_divisions;
            PregnantWomenNeedToScreen(OutputEventNum-1) = pregnantWomenNeedToScreen;
            PeripartumTreatment_HbEAg_HighVL_approx(OutputEventNum-1) = num_mothers_PAP_HbEAg_HighVL;
            PeripartumTreatment_HbEAg_LowVL_approx(OutputEventNum-1) = num_mothers_PAP_HbEAg_LowVL;
            PeripartumTreatment_HbSAg_HighVL_approx(OutputEventNum-1) = num_mothers_PAP_HbSAg_HighVL;
            PeripartumTreatment_HbSAg_LowVL_approx(OutputEventNum-1) = num_mothers_PAP_HbSAg_LowVL;
            RatePeripartumTreatment(OutputEventNum-1) = RateOfPAPInitiation;
            HBVPregnantWomenNeedToEvaluate(OutputEventNum-1) = HBVPositivePregnantWomenAtANC;
            num_starting_treatment_as_eligible(OutputEventNum-1) = num_starting_treatment_as_eligible_this_year; 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % MP: Use this to create a text file of all states (but summed into 5 yr age groups to make it tractable).
        % Create a 4D array of size num_disease_states *  20 (5yr age gps 0-99) * num_sexes * num_treat_blocks. 
        if(store_results_as_text==1)
            if OutputEventNum > 1
                % This is the (consolidated into 5 yr age gp) state variable at time t as
                % a 2D array. We want to store this as a row in 
                X_to_print_unshaped = DUMMY_VALUE * ones(max(agegroups_5yr), ncol_X_to_print_byage);
                deaths_to_print = DUMMY_VALUE * ones(1, max(agegroups_5yr));

                for ag = 1:max(agegroups_5yr) % 1:20
                    
                    temp_store = reshape(squeeze(sum(X(:, agegroups_5yr == ag, :, :),2)), [1,ncol_X_to_print_byage]); % k is gender
                    X_to_print_unshaped(ag,:) = temp_store;

                    %% Now store incident deaths:
                    %% The "(5*(ag-1)+1):(5*(ag-1)+5)" represents the indices converting from 1 year age gps to 5 yr age groups.
                    %% BUG FIX - was previously "(5*(ag-1)+1):(5*(ag-1)+4)"
                    death_age_group_indices = (5*(ag-1)+1):(5*(ag-1)+5);
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
    
    %% Children 1-5:
    FOI(1, i1y:(i5y - 1), :, :) = ...
        ... %% In care (either optimally or non-optimally) but ineligible for treatment:
        + beta_U5_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_not_eAgpos_treat_inelig, i1y:(i5y - 1), :, [i_appropriate_management i_incare_nonadherent]))))) / n_child_1y_5y ...
        + beta_U5_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treat_inelig_under30, i1y:(i5y - 1), :, [i_appropriate_management i_incare_nonadherent]))))) / n_child_1y_5y ...
        ... %% Undiagnosed or out of care
        + beta_U5_SAg(i_dt) * sum(sum(sum(sum(X([i_sAgpos_not_eAgpos_treatelig i_sAgpos_not_eAgpos_treat_inelig], i1y:(i5y - 1), :, [i_undiagnosed, i_outofcare]))))) / n_child_1y_5y ...
        + beta_U5_EAg(i_dt) * sum(sum(sum(sum(X([i_eAgpos_treatelig_under30 i_eAgpos_treat_inelig_under30], i1y:(i5y - 1), :, [i_undiagnosed, i_outofcare]))))) / n_child_1y_5y ...
        ... %% On effective treatment:
        + RRtrans_effective_treatment * beta_U5_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_not_eAgpos_treatelig, i1y:(i5y - 1), :, i_appropriate_management))))) / n_child_1y_5y ...
        + RRtrans_effective_treatment * beta_U5_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treatelig_under30, i1y:(i5y - 1), :, i_appropriate_management))))) / n_child_1y_5y ...
        ... %% On treatment but imperfectly:
        + RRtrans_nonadherent_treatment * beta_U5_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_not_eAgpos_treatelig, i1y:(i5y - 1), :, i_incare_nonadherent))))) / n_child_1y_5y ...
        + RRtrans_nonadherent_treatment * beta_U5_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treatelig_under30, i1y:(i5y - 1), :, i_incare_nonadherent))))) / n_child_1y_5y;
    
    % ii: Transmission between 1-15 year olds
    FOI(1, i1y:(i15y - 1), :, :) = FOI(1, i1y:(i15y - 1), :, :) + ...
        ... %%  In care (either optimally or non-optimally) but ineligible for treatment:
        beta_1to15_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_not_eAgpos_treat_inelig, i1y:(i15y - 1), :, [i_appropriate_management i_incare_nonadherent]))))) / n_child_1y_15y ...
        + beta_1to15_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treat_inelig_under30, i1y:(i15y - 1), :, [i_appropriate_management i_incare_nonadherent]))))) / n_child_1y_15y ...
        ... %% Undiagnosed or out of care
        + beta_1to15_SAg(i_dt) * sum(sum(sum(sum(X([i_sAgpos_not_eAgpos_treatelig i_sAgpos_not_eAgpos_treat_inelig], i1y:(i15y - 1), :, [i_undiagnosed i_outofcare]))))) / n_child_1y_15y ...
        + beta_1to15_EAg(i_dt) * sum(sum(sum(sum(X([i_eAgpos_treatelig_under30 i_eAgpos_treat_inelig_under30], i1y:(i15y - 1), :, [i_undiagnosed i_outofcare]))))) / n_child_1y_15y ...
        ... %% On effective treatment:
        + RRtrans_effective_treatment * beta_1to15_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_not_eAgpos_treatelig, i1y:(i15y - 1), :, i_appropriate_management))))) / n_child_1y_15y ...
        + RRtrans_effective_treatment * beta_1to15_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treatelig_under30, i1y:(i15y - 1), :, i_appropriate_management))))) / n_child_1y_15y ...
        ... %% On treatment but imperfectly:
        + RRtrans_nonadherent_treatment * beta_1to15_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_not_eAgpos_treatelig, i1y:(i15y - 1), :, i_incare_nonadherent))))) / n_child_1y_15y ...
        + RRtrans_nonadherent_treatment * beta_1to15_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treatelig_under30, i1y:(i15y - 1), :, i_incare_nonadherent))))) / n_child_1y_15y;
    
    % iii: Transmission Between 5+ and Adults (Assuming equal risks for all persons 5y-100y)
    %% ALPHA-4
    %% Note that we need to split the eAg+ transmission into <30 and >=30 as treatment eligibility differs:
    FOI(1, i5y:end, :, :) = FOI(1, i5y:end, :, :) + ...
        ... %%  In care (either optimally or non-optimally) but ineligible for treatment:
        beta_5plus_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_not_eAgpos_treat_inelig, i5y:end, :, [i_appropriate_management, i_incare_nonadherent]))))) / n_pop_5y_andabove ...
        ... %% eAg positive eligibility differs by age group (>=30 immune tolerant are eligible):
        + beta_5plus_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treat_inelig_under30, i5y:(i30y-1), :, [i_appropriate_management, i_incare_nonadherent]))))) / n_pop_5y_andabove ...
        + beta_5plus_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treat_inelig_30plus, i30y:num_age_steps, :, [i_appropriate_management, i_incare_nonadherent]))))) / n_pop_5y_andabove ...
        ... %% Undiagnosed or out of care. Note that here we use i_eAgpos to mean all eAg+ (chronic+acute) because we don't need to split by age:
        + beta_5plus_SAg(i_dt) * sum(sum(sum(sum(X([i_sAgpos_not_eAgpos_treatelig i_sAgpos_not_eAgpos_treat_inelig], i5y:end, :, [i_undiagnosed, i_outofcare]))))) / n_pop_5y_andabove ...
        + beta_5plus_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos, i5y:end, :, [i_undiagnosed, i_outofcare]))))) / n_pop_5y_andabove ...   
        ... %% On effective treatment:
        + RRtrans_effective_treatment * beta_5plus_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_not_eAgpos_treatelig, i5y:end, :, i_appropriate_management))))) / n_pop_5y_andabove ...
        ... %% eAg positive eligibility differs by age group (>=30 immune tolerant are eligible):
        + RRtrans_effective_treatment * beta_5plus_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treatelig_under30, i5y:(i30y-1), :, i_appropriate_management))))) / n_pop_5y_andabove ...   
        + RRtrans_effective_treatment * beta_5plus_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treatelig_30plus, i30y:num_age_steps, :, i_appropriate_management))))) / n_pop_5y_andabove ...   
        ... %% On treatment but imperfectly:
        + RRtrans_nonadherent_treatment * beta_5plus_SAg(i_dt) * sum(sum(sum(sum(X(i_sAgpos_not_eAgpos_treatelig, i5y:end, :, i_incare_nonadherent))))) / n_pop_5y_andabove ...
        ... %% eAg positive eligibility differs by age group (>=30 immune tolerant are eligible):
        + RRtrans_nonadherent_treatment * beta_5plus_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treatelig_under30, i5y:(i30y-1), :, i_incare_nonadherent))))) / n_pop_5y_andabove ...   
        + RRtrans_nonadherent_treatment * beta_5plus_EAg(i_dt) * sum(sum(sum(sum(X(i_eAgpos_treatelig_30plus, i30y:num_age_steps, :, i_incare_nonadherent))))) / n_pop_5y_andabove;

    

    % Disease Progression
    %% ALPHA-5. 
    %% There are now different progression matrices Transitions.Values_baseline/future (*TREATMENT STRATUM AND AGE-SPECIFIC* - as treatment eligibility is currently age-specific)
    next_X = X;
    starting_treatment_as_eligible = 0;
    for tr = 1:length(Transitions.From)
        %%transaction_vals = Transitions.Values{tr};
        %%transaction_vals = transaction_vals(:);
        if(time<treatment_rate_params.t_treatment_scaleup_start)
            Current_transition_matrix = Transitions.Values_baseline{tr};
        elseif(time>=treatment_rate_params.t_treatment_scaleup_end)
            Current_transition_matrix = Transitions.Values_future{tr};
        else
            %% Interpolate - mimics effects of scale-up of treatment
            f = (time-treatment_rate_params.t_treatment_scaleup_start) / (treatment_rate_params.t_treatment_scaleup_end - treatment_rate_params.t_treatment_scaleup_start);
            assert(f>=0 & f<=1)
            Current_transition_matrix = (1-f)*Transitions.Values_baseline{tr} + f*Transitions.Values_future{tr};
        end

        
        %%assert(all(all(all(Current_transition_matrix>=Transitions.Values_baseline{tr}))))

        % multiply by dt as Transitions.Values are annual rates.
        moving_btw_states(1, :, :, :) = X(Transitions.From(tr), :, :, :) .* Current_transition_matrix;
        next_X(Transitions.From(tr), :, :, :) = next_X(Transitions.From(tr), :, :, :) - dt * moving_btw_states; % move people out of "from" state
        next_X(Transitions.To(tr), :, :, :)   = next_X(Transitions.To(tr), :, :, :)   + dt * moving_btw_states; % move people into "to" state

        %% Check if this is a treatment initiation (only after 2016):
        if(time>2016)
            if(ismember(Transitions.From(tr),i_treat_inelig_under30) && ismember(Transitions.To(tr), i_treatelig_under30))
                starting_treatment_as_eligible = starting_treatment_as_eligible + ...
                    dt * sum(sum(sum(sum(moving_btw_states(1,1:(i30y-1),:,[i_appropriate_management i_incare_nonadherent])))));
            end
            if(ismember(Transitions.From(tr),i_treat_inelig_30plus) && ismember(Transitions.To(tr), i_treatelig_30plus))
                starting_treatment_as_eligible = starting_treatment_as_eligible +...
                    dt * sum(sum(sum(sum(moving_btw_states(1,i30y:num_age_steps,:,[i_appropriate_management i_incare_nonadherent])))));
            end
        end


    end % end Disease Progression for loop

    %% Stored as a model output:
    num_starting_treatment_as_eligible_this_year = num_starting_treatment_as_eligible_this_year + starting_treatment_as_eligible;


    % Check - no lamivudine treatment any more (and TDF treatment is now dealt with separately) so right now these *must* be
    % zero.
    assert(squeeze(sum(sum(sum(sum(X([i_TDFtreat i_3TCtreat i_3TCfailed], :, :, 1),1),2),3),4))==0)

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ADDITIONAL TESTING (BIRTH COHORT, ANC, COMMUNITY SCREENING)
    %% ALPHA-6 - make this the testing section:
    birth_cohort_testing_start = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'Dx_T_birthcohort_start'),:).Value;
    birth_cohort_testing_end = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'Dx_T_birthcohort_end'),:).Value;
   

    %% Check if any additional screening needed:
    if(~strcmp(scenario_AddScreenIntervention,"No additional screening"))
        if strcmp(scenario_AddScreenIntervention,"Birth cohort screening")
            if (time >= birth_cohort_testing_start && time <= birth_cohort_testing_end)        
                %%case I_NO_COHORT_TEST
                disp("Running birth cohort testing and treatment")
                disp(time)
                BirthCohort_extrayears = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'BirthCohort_extrayears'),:).Value;

                BirthCohort_youngest_birth_year = Intervention_data_thiscountry.BirthCohortTest_year_first_BD + BirthCohort_extrayears;
                %%case I_BIRTHCOHORT_SCREENING
                %% Here we set up the number of people to test at each timestep for birth cohort testing:
                if(time == birth_cohort_testing_start)
                    % This is the index corresponding to the minumum age.
                    % Note that we add 1 as an offset (as age 0 <-> index=1).
                    i_cohortage_min = round((birth_cohort_testing_start-BirthCohort_youngest_birth_year)/dt) + 1; 
                    i_cohortage_max = num_age_steps; % All individuals born before 1992 (no upper age limit).
                    % Check I haven't accidentally made the min age>max age:
                    assert(i_cohortage_max>i_cohortage_min)
                    %% This is the diagnosis coverage we want to achieve in the birth cohort as a whole:
                    birth_cohort_coverage = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'Dx_birthcohort_coverage'),:).Value;
                    %% Duration of birth cohort testing intervention (in years - note we multiply ):
                    duration_birth_cohort_testing = birth_cohort_testing_end-birth_cohort_testing_start;
                    %% Check the duration of testing is OK:
                    assert(duration_birth_cohort_testing>0 && duration_birth_cohort_testing<num_years_simul)
                    
                    %% Note we should use next_X rather than X here:
                    
                    %% Note - the commented code below can be used if we want to set a treatment coverage target rather than a diagnosis target.
                    %% To use - replace N_currentDx_in_birth_cohort with N_current_treatment_in_birth_cohort; with N_treateligible_in_birth_cohort.
                    % if(i_cohortage_min<i30y)
                    %     N_current_treatment_in_birth_cohort = squeeze(sum(sum(sum(sum(next_X(i_treatelig_under30, i_cohortage_min:(i30y-1), :, [i_appropriate_management i_incare_nonadherent]) ...
                    %         + next_X(i_treatelig_30plus, i30y:i_cohortage_max, :, [i_appropriate_management i_incare_nonadherent]))))));
                    %     N_treateligible_in_birth_cohort = squeeze(sum(sum(sum(sum(next_X(i_treatelig_under30, i_cohortage_min:(i30y-1), :, :) ...
                    %         + next_X(i_treatelig_30plus, i30y:i_cohortage_max, :, :))))));
                    % else
                    %    %% Only age30+ in cohort, so eligibility is the same for all:
                    %     N_current_treatment_in_birth_cohort = squeeze(sum(sum(sum(sum(next_X(i_treatelig_30plus, i_cohortage_min:i_cohortage_max, :, [i_appropriate_management i_incare_nonadherent]))))));
                    %      N_treateligible_in_birth_cohort = squeeze(sum(sum(sum(sum(next_X(i_treatelig_30plus, i_cohortage_min:i_cohortage_max, :, :))))));
                    % end

                    %% Diagnosis if sAG positive chronic infection:
                    N_currentDx_in_birth_cohort = squeeze(sum(sum(sum(sum(next_X(i_sAgpos_chronic, i_cohortage_min:i_cohortage_max, :, [i_appropriate_management i_incare_nonadherent i_outofcare]))))));
                    N_sAgpos_in_birth_cohort = squeeze(sum(sum(sum(sum(next_X(i_sAgpos_chronic, i_cohortage_min:i_cohortage_max, :, :))))));

                    assert(N_currentDx_in_birth_cohort>=0)
                    assert(N_sAgpos_in_birth_cohort>0)
                    assert(N_sAgpos_in_birth_cohort>=N_currentDx_in_birth_cohort)
                    
                    %% This is the percentage point increase required (among those not currently on treatment) to reach coverage target.
                    proportion_notcurrentlyDx_toDx = ((N_sAgpos_in_birth_cohort*birth_cohort_coverage) - N_currentDx_in_birth_cohort) / ...
                        (N_sAgpos_in_birth_cohort - N_currentDx_in_birth_cohort);
                    if(proportion_notcurrentlyDx_toDx<0)
                        %% No birth cohort testing if Dx too high - e.g. in China
                        proportion_notcurrentlyDx_toDx = 0;
                    end
                    %% assert(proportion_notcurrentlyDx_toDx>0) %% This doesn't hold for China
                    
                    %% Note that "diagnosis" here means either diagnosing undiagnosed, or finding those out of care and potentially supporting them into care (the cascade is leaky so they may still not reenter care).
                    moving_to_diagnosed_by_birthcohort_testing_per_timestep(i_sAgpos_chronic, i_cohortage_min:i_cohortage_max, :, [i_undiagnosed i_outofcare]) ...
                        = dt * proportion_notcurrentlyDx_toDx * next_X(i_sAgpos_chronic, i_cohortage_min:i_cohortage_max, :, [i_undiagnosed i_outofcare])/duration_birth_cohort_testing;                         
                    disp("Eligible for cohort treatment")
                    disp(squeeze(sum(sum(sum(sum(moving_to_diagnosed_by_birthcohort_testing_per_timestep,1),2),3),4)))
                    disp("At time")
                    disp(time)
                end
                
                %% These are the steps carried out every timestep while the birth cohort testing is happening:

                %% This is the adjustment factor for the age index to account for the time that has passed since the birth cohort started.
                i_birth_cohort_offset = round((time-birth_cohort_testing_start)/dt);
                moving_to_diagnosed_by_birthcohort_testing_this_timestep(:, (i_cohortage_min+i_birth_cohort_offset):i_cohortage_max, :, :) ...
                    = moving_to_diagnosed_by_birthcohort_testing_per_timestep(:, i_cohortage_min:(i_cohortage_max-i_birth_cohort_offset), :, :);
                %% Set the earlier age group elements to zero if needed:
                if(i_birth_cohort_offset>0)
                    moving_to_diagnosed_by_birthcohort_testing_this_timestep(:, i_cohortage_min:(i_cohortage_min+i_birth_cohort_offset-1), :, :) ...
                        = zeros(num_disease_states, i_birth_cohort_offset, num_sexes, num_treat_blocks);
                end
                %% Cap the number of people to move from a given compartment to be at most the number of people in that compartment right now:
                moving_to_diagnosed_by_birthcohort_testing_this_timestep(moving_to_diagnosed_by_birthcohort_testing_this_timestep>next_X) = next_X(moving_to_diagnosed_by_birthcohort_testing_this_timestep>next_X);

                %% Now move people:
                moving_to_diagnosed_by_birthcohort_testing_this_timestep_by_nathist_age_sex = sum(moving_to_diagnosed_by_birthcohort_testing_this_timestep, 4);
                next_X(:, :, :, [i_undiagnosed i_outofcare]) = next_X(:, :, :, [i_undiagnosed i_outofcare]) - moving_to_diagnosed_by_birthcohort_testing_this_timestep(:, :, :, [i_undiagnosed i_outofcare]);
                next_X(:, :, :, i_appropriate_management) = next_X(:, :, :, i_appropriate_management) + prop_adhere_treatment*(1-birthcohort_prop_Dx_outofcare)*moving_to_diagnosed_by_birthcohort_testing_this_timestep_by_nathist_age_sex;
                next_X(:, :, :, i_incare_nonadherent) = next_X(:, :, :, i_incare_nonadherent) + (1-prop_adhere_treatment)*(1-birthcohort_prop_Dx_outofcare)*moving_to_diagnosed_by_birthcohort_testing_this_timestep_by_nathist_age_sex;
                next_X(:, :, :, i_outofcare) = next_X(:, :, :, i_outofcare) + birthcohort_prop_Dx_outofcare*moving_to_diagnosed_by_birthcohort_testing_this_timestep_by_nathist_age_sex;
                %%otherwise
                %%    disp("Error: Unknown value for scenario_AddScreenIntervention. Exiting")
                %%    return
                %% ALPHA-7 - make sure I count the number of tests done. And branch diagnoses into the 3 categories.
                ntests = sum(sum(sum(moving_to_diagnosed_by_birthcohort_testing_this_timestep_by_nathist_age_sex)));
            end 
        elseif strcmp(scenario_AddScreenIntervention,"ANC screening")
            t_ANC = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'Start_ANCtesting'),:).Value;
            if (time >= t_ANC)
                if(time==t_ANC)
                    disp("Running ANC screening")
                end
                %% This is the % of women undergoing ANC testing this timestep in each age group:
                %% ALPHA-8 - bring in the actual ANC testing rates, and make the time ANC testing is introduced (and acceptance) placeholders.
                

                %% scenario_data_ANCHBVtestingbyage_thiscountry has % of women (by 5 yr age group, 15-49) who would go to ANC and would get a test (using HIV testing as a proxy for HBV).
                %% The denominator is all women. This is outputted every 5 years (the original dataset is 1950-2095 but I cut it to 2025-2100 - alter ANC_testing_analysis.R if extra years needed).
                %% Firstly we need to work out what two timepoints (every 5 yrs) we need to interpolate between:
                
                %% For time after 2095 we use the 2095 value.
                if(time<=2095)
                    delta_t = mod(time,5); %% How much time has passed since last multiple of 5 yrs.
                    t0 = time - delta_t;
                    t1 = t0 + 5; %% Next multiple of 5.
                    %% columns 3-end have the ANC testing data:
                    temp_y0 = table2array(scenario_data_ANCHBVtestingbyage_thiscountry(scenario_data_ANCHBVtestingbyage_thiscountry.Year==t0,3:end));
                    temp_y1 = table2array(scenario_data_ANCHBVtestingbyage_thiscountry(scenario_data_ANCHBVtestingbyage_thiscountry.Year==t0,3:end));
                    annual_ANC_test_acceptHBVtest_15_49 = temp_y1*delta_t/5.0 + temp_y0*(5-delta_t)/5.0;
                else
                    annual_ANC_test_acceptHBVtest_15_49 = table2array(scenario_data_ANCHBVtestingbyage_thiscountry(scenario_data_ANCHBVtestingbyage_thiscountry.Year==2095,3:end));
                end

                assert(length(annual_ANC_test_acceptHBVtest_15_49)==7)
                %%annual_ANC_test_acceptHBVtest_15_49 = [0.02,0.04,0.05,0.05,0.02,0.01,0.005];

                %% Add in the 0-14 and 50+ age groups:
                annual_ANC_test_acceptHBVtest = [0,0,0,annual_ANC_test_acceptHBVtest_15_49,0,0,0,0,0,0,0,0,0,0];
                %% Check that this is 20 (0-99 in 5 year age groups)
                assert(length(annual_ANC_test_acceptHBVtest)==20)

                %% percentage of treatment-eligible women who undergo ANC testing and will accept treatment:
                %% PLACEHOLDER - 0.05
                %annual_ANC_test_accepttreat = annual_ANC_test * 0.05;
                

                ANC_testing_by_age = annual_ANC_test_acceptHBVtest(agegroups_5yr);
                %% Percentage of people (by age and sex) getting ANC tested - men + women not going to ANC are zero.
                ANC_testing = zeros(num_disease_states, num_age_steps, num_sexes, num_treat_blocks);

                %% Testing over chronic stages:
                ANC_testing(i_ImmTol, : ,i_female, i_undiagnosed) = ANC_testing_by_age;
                ANC_testing(i_ImmReact,:,i_female,i_undiagnosed) = ANC_testing_by_age;
                ANC_testing(i_AsymptCarr,:,i_female,i_undiagnosed) = ANC_testing_by_age;
                ANC_testing(i_Chronic,:,i_female,i_undiagnosed) = ANC_testing_by_age;
                ANC_testing(i_CompCirr,:,i_female,i_undiagnosed) = ANC_testing_by_age;
                ANC_testing(i_DecompCirr,:,i_female,i_undiagnosed) = ANC_testing_by_age;
                ANC_testing(i_HCC,:,i_female,i_undiagnosed) = ANC_testing_by_age;
                %% Assume ANC testing also reaches those previously diagnosed but out of care:
                ANC_testing(i_ImmTol, : ,i_female, i_outofcare) = ANC_testing_by_age;
                ANC_testing(i_ImmReact,:,i_female,i_outofcare) = ANC_testing_by_age;
                ANC_testing(i_AsymptCarr,:,i_female,i_outofcare) = ANC_testing_by_age;
                ANC_testing(i_Chronic,:,i_female,i_outofcare) = ANC_testing_by_age;
                ANC_testing(i_CompCirr,:,i_female,i_outofcare) = ANC_testing_by_age;
                ANC_testing(i_DecompCirr,:,i_female,i_outofcare) = ANC_testing_by_age;
                ANC_testing(i_HCC,:,i_female,i_outofcare) = ANC_testing_by_age;


                %% POWIETRZE
                %% Note we should use next_X rather than X here:
                moving_to_diagnosed_by_ANC_testing_this_timestep = dt * ANC_testing .* next_X; 
                %% This line caps the number of people moving at this timestep in a given compartment to be at most next_X in that compartment.
                moving_to_diagnosed_by_ANC_testing_this_timestep(moving_to_diagnosed_by_ANC_testing_this_timestep>next_X) = next_X(moving_to_diagnosed_by_ANC_testing_this_timestep>next_X);
                               
                % if(time<2028)
                %     fprintf("Eligible for ANC treatment: %6.4f at time %6.4f\n",sum(sum(sum(sum(moving_to_diagnosed_by_ANC_testing_this_timestep,1),2),3),4), time)
                % end
                %disp("Eligible for ANC treatment")
                %disp(sum(sum(sum(sum(moving_to_diagnosed_by_ANC_testing_per_timestep,1),2),3),4))
                %disp("At time")
                %disp(time)
                %fprintf("Eligible for ANC treatment %d at time %d",sum(sum(sum(sum(moving_to_diagnosed_by_ANC_testing_this_timestep,1),2),3),4),time)

                sum_moving_to_diagnosed_by_ANC_testing_this_timestep = sum(moving_to_diagnosed_by_ANC_testing_this_timestep, 4);
                next_X(:, :, :, [i_undiagnosed i_outofcare]) = next_X(:, :, :, [i_undiagnosed i_outofcare]) - moving_to_diagnosed_by_ANC_testing_this_timestep(:, :, :, [i_undiagnosed i_outofcare]);
                next_X(:, :, :, i_appropriate_management) = next_X(:, :, :, i_appropriate_management) + prop_adhere_treatment*(1-ANCtesting_prop_Dx_outofcare)*sum_moving_to_diagnosed_by_ANC_testing_this_timestep;
                next_X(:, :, :, i_incare_nonadherent) = next_X(:, :, :, i_incare_nonadherent) + (1-prop_adhere_treatment)*(1-ANCtesting_prop_Dx_outofcare)*sum_moving_to_diagnosed_by_ANC_testing_this_timestep;
                next_X(:, :, :, i_outofcare) = next_X(:, :, :, i_outofcare) + ANCtesting_prop_Dx_outofcare*sum_moving_to_diagnosed_by_ANC_testing_this_timestep;
            end
        elseif (strcmp(scenario_AddScreenIntervention,"Community screening") || strcmp(scenario_AddScreenIntervention,"Perfect community screening"))
            if strcmp(scenario_AddScreenIntervention,"Community screening")
                community_screening_coverage = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'DX_community_screening_coverage'),:).Value;
            elseif strcmp(scenario_AddScreenIntervention,"Perfect community screening")
                community_screening_coverage = 1.0;
            else
                disp("Error - Unknown screening scenario. Exiting")
                return
            end
           

            community_screening_start = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'Dx_T_community_screening_start'),:).Value;
            community_screening_end = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'Dx_T_community_screening_end'),:).Value;
            community_screening_min_age = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'Dx_community_screening_min_age'),:).Value;
            community_screening_max_age = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'Dx_community_screening_max_age'),:).Value;
            i_screeningage_min = round(community_screening_min_age/dt); % Individuals born before 1992
            i_screeningage_max = round(community_screening_max_age/dt); % All individuals born before 1992 (no upper age limit).
            % Check I haven't accidentally made the min age>max age:
            assert(i_screeningage_max>i_screeningage_min)
            assert(i_screeningage_max<=num_age_steps)
            
            if (time >= community_screening_start && time <= community_screening_end)        
                %% At the beginning work out how many people need to move each timestep:
                if(time == community_screening_start)
                    duration_community_screening = community_screening_end-community_screening_start;
                    %% Check the duration of testing is OK:
                    assert(duration_community_screening>0 && duration_community_screening<num_years_simul)
                    
                    %% Use this code if the target for community screening is % on treatment rather than % diagnosed.
                    % if(i_screeningage_min<i30y)
                    %     N_current_treatment_by_community_screening = squeeze(sum(sum(sum(sum(next_X(i_treatelig_under30, i_screeningage_min:(i30y-1), :, [i_appropriate_management i_incare_nonadherent]) ...
                    %         + next_X(i_treatelig_30plus, i30y:i_screeningage_max, :, [i_appropriate_management i_incare_nonadherent]))))));
                    %     N_treateligible_by_community_screening = squeeze(sum(sum(sum(sum(next_X(i_treatelig_under30, i_screeningage_min:(i30y-1), :, :) ...
                    %         + next_X(i_treatelig_30plus, i30y:i_screeningage_max, :, :))))));
                    % else
                    %     %% Only age30+ in cohort, so eligibility is the same for all:
                    %     N_current_treatment_by_community_screening = squeeze(sum(sum(sum(sum(next_X(i_treatelig_30plus, i_screeningage_min:i_screeningage_max, :, [i_appropriate_management i_incare_nonadherent]))))));
                    %     N_treateligible_by_community_screening = squeeze(sum(sum(sum(sum(next_X(i_treatelig_30plus, i_screeningage_min:i_screeningage_max, :, :))))));
                    % end

                    N_currentDx_in_community_screening = squeeze(sum(sum(sum(sum(next_X(i_sAgpos_chronic, i_screeningage_min:i_screeningage_max, :, [i_appropriate_management i_incare_nonadherent i_outofcare]))))));
                    N_sAgpos_in_community_screening = squeeze(sum(sum(sum(sum(next_X(i_sAgpos_chronic, i_screeningage_min:i_screeningage_max, :, :))))));
                    
                    assert(N_currentDx_in_community_screening>=0)
                    assert(N_sAgpos_in_community_screening>0)
                    assert(N_sAgpos_in_community_screening>=N_currentDx_in_community_screening)

                    %% This is the percentage point increase required (among those not currently diagnosed) to reach coverage target.
                    proportion_notcurrentlyDx_toDx = ((N_sAgpos_in_community_screening*community_screening_coverage) - N_currentDx_in_community_screening) / ...
                        (N_sAgpos_in_community_screening - N_currentDx_in_community_screening);
                    %% No community screening if Dx too high - e.g. in China
                    if(proportion_notcurrentlyDx_toDx<0)
                        proportion_notcurrentlyDx_toDx = 0;
                    end
                    %%assert(proportion_notcurrentlyDx_toDx>0)
                    
                    %% Note that "diagnosis" here means either diagnosing undiagnosed, or finding those out of care and potentially supporting them into care (the cascade is leaky so they may still not reenter care).
                    %% Note we should use next_X rather than X here:
                    moving_to_diagnosed_by_community_screening_per_timestep(i_sAgpos_chronic, i_screeningage_min:i_screeningage_max, :, [i_undiagnosed i_outofcare]) ...
                        = dt * proportion_notcurrentlyDx_toDx * next_X(i_sAgpos_chronic, i_screeningage_min:i_screeningage_max, :, [i_undiagnosed i_outofcare])/duration_community_screening;

                    disp("Community screening eligibles")
                    disp(squeeze(sum(sum(sum(sum(moving_to_diagnosed_by_community_screening_per_timestep,1),2),3),4)))
                    disp("At time")
                    disp(time)
                end
                
                assert(time<=community_screening_end);
                %% This is the adjustment factor for the age index to account for the time that has passed since the birth cohort started.
                i_community_screening = round((time-community_screening_start)/dt);
                i_max = min((i_screeningage_max+i_community_screening),num_age_steps);
                %% For interventions with no maximum, this ensures we don't exceed the highest index:
                i_trunc = max((i_screeningage_max+i_community_screening-i_max),0);
                
                %% Copy the number of people who need to be screened this timestep (taking into account the fact that they have aged):
                moving_to_diagnosed_by_community_screening_this_timestep(:, (i_screeningage_min+i_community_screening):i_max, :, :) ...
                    = moving_to_diagnosed_by_community_screening_per_timestep(:, i_screeningage_min:(i_screeningage_max-i_trunc), :, :);


                %% Set the earlier age group elements to zero if needed:
                if(i_community_screening>0)
                    moving_to_diagnosed_by_community_screening_this_timestep(:, i_screeningage_min:(i_screeningage_min+i_community_screening-1), :, :) ...
                        = zeros(num_disease_states, i_community_screening, num_sexes, num_treat_blocks);
                end
                %% Ensure we never go below 0:
                %% Firstly, for any elements of moving_to_diagnosed_by_community_screening_this_timestep which are > than the corresponding element in next_X, set that to be the value in next_X:
                moving_to_diagnosed_by_community_screening_this_timestep(moving_to_diagnosed_by_community_screening_this_timestep>next_X) = next_X(moving_to_diagnosed_by_community_screening_this_timestep>next_X);

                
                %% Now move people:
                moving_to_diagnosed_by_community_screening_this_timestep_by_nathist_age_sex = sum(moving_to_diagnosed_by_community_screening_this_timestep, 4);
                next_X(:, :, :, [i_undiagnosed i_outofcare]) = next_X(:, :, :, [i_undiagnosed i_outofcare]) - moving_to_diagnosed_by_community_screening_this_timestep(:, :, :, [i_undiagnosed i_outofcare]);
                next_X(:, :, :, i_appropriate_management) = next_X(:, :, :, i_appropriate_management) + prop_adhere_treatment*(1-communityscreening_prop_Dx_outofcare)*moving_to_diagnosed_by_community_screening_this_timestep_by_nathist_age_sex;
                next_X(:, :, :, i_incare_nonadherent) = next_X(:, :, :, i_incare_nonadherent) + (1-prop_adhere_treatment)*(1-communityscreening_prop_Dx_outofcare)*moving_to_diagnosed_by_community_screening_this_timestep_by_nathist_age_sex;
                next_X(:, :, :, i_outofcare) = next_X(:, :, :, i_outofcare) + communityscreening_prop_Dx_outofcare*moving_to_diagnosed_by_community_screening_this_timestep_by_nathist_age_sex;

                %%otherwise
                %%    disp("Error: Unknown value for scenario_AddScreenIntervention. Exiting")
                %%    return

            end %% End of community screening

        end %% End of if statement looping through different screening options.
    end   %% End of scenario_AddScreenIntervention!="No screening

    %% TREATMENT:
    if (time >= treat_start_year && scenario_Treatment>0)
    % 2016 must be the first year with nonzero treatment 
    % therefore start treating from 2015.9 onwards since prevalence is recorded at the top of the loop

        if ~initiated_treatment
            %% Count up the number of people on treatment in 2016 and check this is zero:
            num_in_treatment = sum(sum(sum(sum(X(i_treatelig_under30, 1:(i30y-1), :, [i_appropriate_management i_incare_nonadherent]),1),2),3),4) + ...
                                    sum(sum(sum(sum(X(i_treatelig_30plus, i30y:num_age_steps, :, [i_appropriate_management i_incare_nonadherent]),1),2),3),4);
            assert(num_in_treatment==0) % no one is in treatment

            %% If data says treatment coverage is >0 in 2016, work out how many people need to be on treatment (and transfer them).
            %% Note that the denom for treatment_rate_params.Tx_coverage_2016 is all sAgpos (including acute!):
            %% in the MJdV code the line is prev_pop = sum(sum(sum(sum(X([2:8 10 12:15], :, :, :),1),2),3),4); - so includes TDF treatment (10) and acute (14:15)
            if(treatment_rate_params.Tx_coverage_2016>0)
                prev_pop = sum(sum(sum(sum(X(i_sAgpos, :, :, :),1),2),3),4); %% Whole pop of sAg+ (including those not on treatment/never diagnosed)
    
                %% Note - prev_pop is whole pop, so that treatment_rate_params.Tx_coverage_2016 is coverage in the whole population of sAg+.
                total_num_to_move_to_treat = treatment_rate_params.Tx_coverage_2016 * prev_pop;

                %% Eligible pop is undiagnosed (in 2016 this should be everyone as we haven't modelled testing pre-2016):
                eligible_pop = squeeze(sum(sum(sum(X(i_treatelig_under30, 1:(i30y-1), :, i_undiagnosed),1),2),3)) + ...
                                    squeeze(sum(sum(sum(X(i_treatelig_30plus, i30y:num_age_steps, :, i_undiagnosed),1),2),3)); 
                
                %% For low coverages it is possible to have 0 eligible but >0 coverage (due to rounding) so only care if >1%:
                % if(treatment_rate_params.Tx_coverage_2016>=0.01)
                %     if((total_num_to_move_to_treat>eligible_pop) && (total_num_to_move_to_treat<(6*eligible_pop)))
                %         total_num_to_move_to_treat = eligible_pop;
                %     end

                assert(total_num_to_move_to_treat<=eligible_pop)
                % end

                %% Scaling_num is an adjustment for the fact that not everyone who is sAg+ (and chronic) is eligible:
                if(eligible_pop>0)
                    scaling_num = total_num_to_move_to_treat / eligible_pop;
                else
                    scaling_num = 0;
                end

                %% Number of people who get moved to treatment when treatment first starts (in 2016) to match data on coverage at that time:
                n_to_move_treatment_start_year = zeros(size(X));              

                % n_to_move_treatment_start_year = min(next_X(i_treateligible,:,:,i_seektreat), X(i_treateligible,:,:,i_seektreat) * scaling_num);
                %%next_X(i_treateligible,:,:,i_seektreat)=next_X(i_treateligible,:,:,i_seektreat) - X(i_treateligible,:,:,i_seektreat) * scaling_num;
                % next_X(i_treateligible,:,:,i_seektreat)=next_X(i_treateligible,:,:,i_seektreat) - n_to_move_treatment_start_year;

                %% Calculate the number of people who need to move, and then
                n_to_move_treatment_start_year(i_treatelig_under30, 1:(i30y-1), :, i_undiagnosed) = X(i_treatelig_under30, 1:(i30y-1), :, i_undiagnosed) * scaling_num;
                n_to_move_treatment_start_year(i_treatelig_30plus, i30y:num_age_steps, :, i_undiagnosed) = X(i_treatelig_30plus, i30y:num_age_steps, :, i_undiagnosed) * scaling_num;
                %% Previous version: n_to_move_treatment_start_year = X(i_treateligible,:,:,i_seektreat) * scaling_num;


                %% Remove people from undiagnosed to either appropriate_management (adheres to treatment) or incare_noadherent (doesn't adhere)
                next_X(:,:,:,i_undiagnosed) = next_X(:,:,:,i_undiagnosed) - n_to_move_treatment_start_year(:,:,:,i_undiagnosed);
                next_X(:,:,:,i_appropriate_management) = next_X(:,:,:,i_appropriate_management) + prop_adhere_treatment*n_to_move_treatment_start_year(:,:,:,i_undiagnosed);
                next_X(:,:,:,i_incare_nonadherent)     = next_X(:,:,:,i_incare_nonadherent) + (1-prop_adhere_treatment)*n_to_move_treatment_start_year(:,:,:,i_undiagnosed);

                
                % Every compartment in the eligible-for-treatment states in next_X must have a number subtracted from it 
                % such that the total number subtracted from the eligible-for-treatment states is in_treatment_2016
                % i.e. in_treatment_2016 = sum(sum(sum(sum(X(i_treateligible, :, :, :),1),2),3),4) * scaling_num = sum(sum(sum(sum(X(i_treateligible, :, :, :) * scaling_num,1),2),3),4)
                % Hence, scaling_num scales each compartment in X(i_treateligible, :, :, :) such that X(i_treateligible,:,:,:) * scaling_num subtracts the same proportion of people from each compartment in each of the eligible-for-treatment states in order to subtract a total of in_treatment_2016 from the eligible-for-treatment states.

                %next_X(i_TDFtreat,:,:,i_seektreat) = next_X(i_TDFtreat,:,:,i_seektreat) + sum(X(i_treateligible,:,:,i_seektreat) * scaling_num,1);

                num_in_treatment = sum(sum(sum(sum(next_X(i_treatelig_under30, 1:(i30y-1), : , [i_appropriate_management i_incare_nonadherent]),1),2),3),4) ...
                    + sum(sum(sum(sum(next_X(i_treatelig_30plus, i30y:num_age_steps, :, [i_appropriate_management i_incare_nonadherent]),1),2),3),4);

                %% This represents the number of people starting treatment at this timestep (as noone is on treatment in the model before 2016).
                number_starting_treatment_to_print = num_in_treatment;
                
                %% Now just double-check everything again:
                eligible_pop = squeeze(sum(sum(sum(sum(X(i_treatelig_under30, 1:(i30y-1), :, :),1),2),3),4)) + ...
                                    squeeze(sum(sum(sum(sum(X(i_treatelig_30plus, i30y:num_age_steps, :, :),1),2),3),4)); 
                assert(num_in_treatment/eligible_pop >= treatment_rate_params.Tx_coverage_2016)

                % treatment coverage amongst treatment-eligible people will be greater than treatment coverage amongst HBsAg+ people, except if treatment coverage is 0
                %%treat_coverage_2016 = num_in_treatment / eligible_pop; %% MP: CHECK WITH SHEVANTHI - THIS IS CURRENTLY DEAD CODE.
            end
            initiated_treatment = true;
        else
            assert(initiated_treatment) % ensure that, each time this code is encountered, treatment has already been initiated
            %% Treatment rate is the rate that people start treatment. *FIXME* - make this comment more informative.
            if (time<=treatment_rate_params.t_treatment_scaleup_start)
                annual_increase_Dx = treatment_rate_params.annual_increase_Dx_past;
                annual_increase_TxifDx = treatment_rate_params.annual_increase_TxifDx_past;
            elseif (time>=treatment_rate_params.t_treatment_scaleup_end)
                annual_increase_Dx = treatment_rate_params.annual_increase_Dx_future;
                annual_increase_TxifDx = treatment_rate_params.annual_increase_TxifDx_future;
            else
                temp_tscale = (time-treatment_rate_params.t_treatment_scaleup_start)/(treatment_rate_params.t_treatment_scaleup_end - treatment_rate_params.t_treatment_scaleup_start);
                annual_increase_Dx = treatment_rate_params.annual_increase_Dx_past + (treatment_rate_params.annual_increase_Dx_future - treatment_rate_params.annual_increase_Dx_past) * temp_tscale;
                annual_increase_TxifDx = treatment_rate_params.annual_increase_TxifDx_past + (treatment_rate_params.annual_increase_TxifDx_future - treatment_rate_params.annual_increase_TxifDx_past) * temp_tscale;

            end
            assert(isscalar(annual_increase_Dx))
            assert(annual_increase_Dx>=0)
            assert(isscalar(annual_increase_TxifDx))
            assert(annual_increase_TxifDx>=0)
            
            %% LECZENIE:
            %% Determine if we are now increasing the number of people who would seek treatment if necessary (by removing barriers to testing/treatment e.g. through decentralisation, integration):
            % if(time>=treatment_rate_params.t_remove_treatment_barriers)
            %     %% Check if the proportion currently seeking treatment is already above the threshold:
            %     prop_currently_seek_treat = sum(sum(sum(sum(next_X(:, :, :, i_seektreat),1),2),3),4)/sum(sum(sum(sum(next_X(:, :, :, :),1),2),3),4);
            %     disp(scenario_treat_elig)
            %     % if(prop_currently_seek_treat<max_treatment_coverage)
            %         prop_inc_seek_treatment = dt*(treatment_rate_params.annual_increase_diagnosis*treatment_rate_params.annual_increase_treatifdiag);
            % end
            
            %% Total chronic HBV (denominator for Dx)
            n_chronic = sum(sum(sum(sum(X(i_sAgpos_chronic, :, :, :),1),2),3),4);
            %%n_undiagnosed = squeeze(sum(sum(sum(X(i_sAgpos_chronic, :, :, i_undiagnosed),1),2),3));

            %% Total diagnosed (denominator for TxifDx):
            %% n_diagnosed = n_chronic - n_undiagnosed; %% Note - I tried comparing this with the calculation for n_diagnosed below, and they were different by 1e-11 - presumably numerical error.
            %% Stupid check:
            n_diagnosed = sum(sum(sum(sum(X(i_sAgpos_chronic, :, :, [i_appropriate_management, i_incare_nonadherent, i_outofcare])))));
            % sum(sum(sum(sum(n_diagnosed-temp))))
            % assert(n_diagnosed==temp)

            
            %% This is the number of people who aren't on treatment but are eligible (so the pool of people we could move onto treatment).
            %% Note - moving these people onto treatment is a mix of people getting diagnosed and treated immediately (as eligible) and people out of care getting linked back into care and starting Tx.
            n_eligible_notonTx = sum(sum(sum(sum(X(i_treatelig_under30, 1:(i30y-1), :, [i_undiagnosed, i_outofcare]),1),2),3),4) ...
                + sum(sum(sum(sum(X(i_treatelig_30plus, i30y:num_age_steps, :, [i_undiagnosed, i_outofcare]),1),2),3),4);

            %% We multiply by dt to calculate the number of people to move this timestep later on:
            n_to_diagnose_thisyear = n_chronic*annual_increase_Dx;
            %% This is the number of people who are diagnosed and (more or less) immediately start treatment.
            %% We have to remove those who were already in care and started treatment.
            %% Note that starting_treatment_as_eligible is in this timestep (hence divide by dt here), while annual_increase_TxifDx is an annual rate. We convert n_to_start_treatment_directly_thisyear to be per timestep later.
            n_to_start_treatment_directly_thisyear = n_diagnosed*annual_increase_TxifDx - (starting_treatment_as_eligible/dt);

            if(n_to_start_treatment_directly_thisyear>0)
                prop_eligbutnotTx_to_treat = n_to_start_treatment_directly_thisyear/n_eligible_notonTx;
                if(prop_eligbutnotTx_to_treat>1)
                    %% Cap diagnosis at 1;
                    prop_eligbutnotTx_to_treat = 1;
                end
            else
                %% This always occurs for SQ scenario, but print warning if it happens at other times:
                if(~scenario_Treatment==I_TREAT.SQ)
                    disp(time)
                    disp("Warning: Nobody left to treat!")
                end
                prop_eligbutnotTx_to_treat = 0;
            end
            
            n_undiagnosed = squeeze(sum(sum(sum(X(i_sAgpos_chronic, :, :, i_undiagnosed),1),2),3));
            if(n_undiagnosed>0)
                prop_undiagnosed_to_diagnose = n_to_diagnose_thisyear/n_undiagnosed;
                if(prop_undiagnosed_to_diagnose>1)
                    %% Cap diagnosis at 1;
                    prop_undiagnosed_to_diagnose = 1;
                end
            else
                %% Hopefully this will never happen - if the warning appears, check what is causing this.
                disp(time)
                disp("Warning: Nobody left to diagnose!")
                prop_undiagnosed_to_diagnose = 0;
            end

            moving_to_diagnosed(i_sAgpos_chronic, :, :, i_undiagnosed) = prop_undiagnosed_to_diagnose * X(i_sAgpos_chronic, :, :, i_undiagnosed);
            moving_to_treatment(i_treatelig_under30, 1:(i30y-1), :, [i_undiagnosed, i_outofcare]) = prop_eligbutnotTx_to_treat * X(i_treatelig_under30, 1:(i30y-1), :, [i_undiagnosed, i_outofcare]);
            moving_to_treatment(i_treatelig_30plus, i30y:num_age_steps, :, [i_undiagnosed, i_outofcare]) = prop_eligbutnotTx_to_treat * X(i_treatelig_30plus, i30y:num_age_steps, :, [i_undiagnosed, i_outofcare]);

            %% We need to have a number for the % of people who remain in care (versus leave care) - for those who aren't treatment-eligible.
            %% We approximate this as the fraction of those diagnosed who are treatment-eligible who start treatment (so remain in care)
            %% Firstly store the number of people who are treatment-eligible who get diagnosed at this timestep:
            n_treatelig_whoarediagnosed = sum(sum(sum(sum(moving_to_diagnosed(i_treatelig_under30, 1:(i30y-1), :, i_undiagnosed))))) + sum(sum(sum(moving_to_diagnosed(i_treatelig_30plus, i30y:num_age_steps, :, i_undiagnosed))));
            %% This is the number of people who get diagnosed this timestep who start treatment:
            n_starttreat_whoarediagnosed = sum(sum(sum(sum(moving_to_treatment(i_treatelig_under30, 1:(i30y-1), :, i_undiagnosed))))) + sum(sum(sum(sum(moving_to_treatment(i_treatelig_30plus, i30y:num_age_steps, :, i_undiagnosed)))));
            if(n_treatelig_whoarediagnosed>0)
                prop_remain_in_care = n_starttreat_whoarediagnosed/n_treatelig_whoarediagnosed;
            end

            %% Now remove anyone who goes from undiagnosed direct to treatment 
            moving_to_diagnosed(i_treatelig_under30, 1:(i30y-1), :, i_undiagnosed) = moving_to_diagnosed(i_treatelig_under30, 1:(i30y-1), :, i_undiagnosed) - moving_to_treatment(i_treatelig_under30, 1:(i30y-1), :, i_undiagnosed);
            moving_to_diagnosed(i_treatelig_30plus, i30y:num_age_steps, :, i_undiagnosed) = moving_to_diagnosed(i_treatelig_30plus, i30y:num_age_steps, :, i_undiagnosed) - moving_to_treatment(i_treatelig_30plus, i30y:num_age_steps, :, i_undiagnosed);


            %% Ensure the number moving is not negative:
            moving_to_diagnosed(moving_to_diagnosed<0) = 0;
            %% Stupid checks:
            % if(max(max(max(max(moving_to_diagnosed))))<=0)
            %     fprintf("At t=%6.4f there is %6.4f-%6.4f to treat",time,min(min(min(min(moving_to_diagnosed)))), max(max(max(max(moving_to_diagnosed)))))
            % end
            assert(max(max(max(max(moving_to_diagnosed))))>=0)
            assert(min(min(min(min(moving_to_diagnosed))))>=0)

            
            next_X(:, :, :, i_undiagnosed)         = next_X(:, :, :, i_undiagnosed) - dt * moving_to_diagnosed(:, :, :, i_undiagnosed) - dt * moving_to_treatment(:, :, :, i_undiagnosed);
            next_X(:,:,:,i_appropriate_management) = next_X(:,:,:,i_appropriate_management) + dt * prop_remain_in_care * prop_adhere_treatment * moving_to_diagnosed(:,:,:,i_undiagnosed) ...
                + dt * prop_adhere_treatment * moving_to_treatment(:, :, :, i_undiagnosed);
            next_X(:,:,:,i_incare_nonadherent)     = next_X(:,:,:,i_incare_nonadherent) + dt * prop_remain_in_care * (1-prop_adhere_treatment) * moving_to_diagnosed(:,:,:,i_undiagnosed) ...
                + dt * (1-prop_adhere_treatment) * moving_to_treatment(:, :, :, i_undiagnosed);
            next_X(:,:,:,i_outofcare)              = next_X(:,:,:,i_outofcare) + dt * (1-prop_remain_in_care) * moving_to_diagnosed(:,:,:,i_undiagnosed);
            
            number_starting_treatment_to_print = squeeze(sum(sum(sum(sum(moving_to_diagnosed, 1), 2), 3), 4));
            assert(isscalar(number_starting_treatment_to_print))

        end
    end % end treatment if statement

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Infection process
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NewInfections = X(i_Susc, :, :, :) .* FOI;
    % number of susceptibles times FOI i.e. number of new infections within population, excluding babies (since FOI is 0 for babies)
    % 1 x num_age_steps x 2 x 2 double i.e. 1 x 1000 x 2 x 2 double
    NewChronicCarriage = NewInfections .* p_ChronicCarriage;
    SevereAcute = NewInfections * theta;
    NonsevereAcute = NewInfections - SevereAcute;
    
    % Transition dependent on a state that does not have a number and is therefore not in Prog or Transitions
    % multiply by dt, since FOI is an annual rate
    next_X(i_Susc, :, :, :) = next_X(i_Susc, :, :, :) - dt * NewInfections;
    next_X(i_NonSevAcute, :, :, :) = next_X(i_NonSevAcute, :, :, :) + dt * NonsevereAcute;
    next_X(i_SevereAcute, :, :, :) = next_X(i_SevereAcute, :, :, :) + dt * SevereAcute;
    
    
    % Infant vaccination HepB3:
    % Do not multiply by dt, since one is vaccinating scenario_HepB3coverage(i_dt)% of people in next_X(1, i6mo, :, :), 
    % after which this cohort ages and moves to the next age bin
    % if divides all babies born in a year into 10 groups and vaccinatates scenario_HepB3coverage(i_dt)% of each group, 
    % then one will have vaccinated scenario_HepB3coverage(i_dt)% of all
    % babies born in that year.
    %% In the current model, babies are all in the i_undiagnosed treatment/care stratum, so we could simplify the expression below slightly.
    %% However, keep as-is in case we ever modify the treatment/care strata again.
    transfer_to_HepB3vacc = scenario_HepB3coverage(i_dt) * next_X(i_Susc, i6mo, :, :) * params.Efficacy_InfantVacc; % the 0.95 represent a take-type vaccine efficacy of 95%.
    next_X(i_Susc, i6mo, :, :) = next_X(i_Susc, i6mo, :, :) - transfer_to_HepB3vacc;
    next_X(i_Immune, i6mo, :, :) = next_X(i_Immune, i6mo, :, :) + transfer_to_HepB3vacc;
    
        
    
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
    
    % Now age everyone by one timestep (second index is the age index with 1=newborn in this timestep).
    X(:, 2:num_age_steps, :, :) = X(:, 1:(num_age_steps - 1), :, :);
    X(:, 1, :, :) = 0; % set number of new babies (the age index "1") to 0 (babies will be born next)
    
    

    % fill-out with new births in this time-step:
    %% MP: Magic numbers 1 and 4 mean sum over the listed natural history states and treatment states
    births_toNonInfectiousWomen = sum( fert' .* sum(sum(X([i_Susc i_Immune], :, i_female, :), 1), 4) ); % Susecptible, Immune
    %%births_toHbEAgWomen = sum(fert' .* sum(sum(X(i_eAgpos, :, i_female, :), 1), 4)); % Immune Tolerant, Immune Reactive
    %%births_toHbSAgWomen = sum(fert' .* sum(sum(X(i_sAgpos_notEagpos_notreat, :, i_female, :), 1), 4)); % All other stages (other infected women)
    
    %% In the PAP model these are incorporated in births_toNonInfectiousWomen. 
    %% As treatment will reduce VL we don't bother stratifying by high/low VL here:
    %% ALPHA-11 - change to strata of treatment (by eligibility!):
    
    % params.fert is a 1000 x (num_years_simul+1) matrix; ages in 0.1 year jumps versus 212 years
    %% fert = params.fert(1:10:end, OutputEventNum);

    %%births_toTrWomen = sum(fert' .* sum(sum(X([i_TDFtreat i_3TCtreat], :, i_female, :), 1), 4)); % Women on Treatment
    % Women on Treatment - split into under 30 and 30+:
    births_toTrWomen = sum(fert(1:(i30y-1))' .* sum(sum(X(i_treatelig_under30, 1:(i30y-1), i_female, [i_appropriate_management i_incare_nonadherent]), 1), 4)) ...
        + sum(fert(i30y:num_age_steps)' .* sum(sum(X(i_treatelig_30plus, i30y:num_age_steps, i_female, [i_appropriate_management i_incare_nonadherent]), 1), 4));
    
    n_births_toHbEAgWomen_not_on_treatment = sum(fert' .* sum(sum(X(i_eAgpos, :, i_female, [i_undiagnosed i_outofcare]), 1), 4)) ...
        + sum(fert(1:(i30y-1))' .* sum(sum(X(i_eAgpos_treat_inelig_under30, 1:(i30y-1), i_female, [i_appropriate_management i_incare_nonadherent]), 1), 4)) ...
        + sum(fert(i30y:num_age_steps)' .* sum(sum(X(i_eAgpos_treat_inelig_30plus, i30y:num_age_steps, i_female, [i_appropriate_management i_incare_nonadherent]), 1), 4));
    
    n_births_toHbsAg_not_eAg_women_not_on_treatment = sum(fert' .* sum(sum(X(i_sAgpos_not_eAgpos_treatelig, :, i_female, [i_undiagnosed i_outofcare]), 1), 4)) ...
        + sum(fert' .* sum(sum(X(i_sAgpos_not_eAgpos_treat_inelig, :, i_female, :), 1), 4));
    
    births_toHbEAgWomenHighVL = PAP_VL_params.FracEPosHighVL    * n_births_toHbEAgWomen_not_on_treatment;
    births_toHbEAgWomenLowVL = (1-PAP_VL_params.FracEPosHighVL) * n_births_toHbEAgWomen_not_on_treatment; 
    births_toHbSAgWomenHighVL = PAP_VL_params.FracSPosHighVL    * n_births_toHbsAg_not_eAg_women_not_on_treatment;
    births_toHbSAgWomenLowVL = (1-PAP_VL_params.FracSPosHighVL) * n_births_toHbsAg_not_eAg_women_not_on_treatment;
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

    %% ratebirthdoses is normal BD (i.e. not using MAP/CPAD):
    ratebirthdoses = births_Total * scenario_BDcoverage(i_dt);
    ratebirthdoses_MAP = births_Total * scenario_BDcoverage_fromMAP(i_dt);
    ratebirthdoses_CPAD = births_Total * scenario_BDcoverage_fromCPAD(i_dt);
    

    babies_ChronicCarriage = babiesChronic_from_HbEAgWomenHighVL + babiesChronic_from_HbEAgWomenLowVL + ...
        babiesChronic_from_HbSAgWomenHighVL + babiesChronic_from_HbSAgWomenLowVL + babiesChronic_from_HbEAgTrWomen;

    babies_NotChronicCarriage = births_Total - babies_ChronicCarriage;


    assert(isscalar(babies_ChronicCarriage))
    assert(isscalar(babies_NotChronicCarriage))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PAP chunk 2 - should replace above.


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
    %% End of PAP Chunk 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    female_multiplier = 1 / (1 + sex_ratio);
    male_multiplier = sex_ratio / (1 + sex_ratio);
    % sex_ratio is number of male births per one female birth
    % 1 / (1 + sex_ratio) + sex_ratio / (1 + sex_ratio) = 1, hence total number of babies not changed
    % male births -> 0 => sex_ratio -> 0 => male_multiplier -> 0/1 = 0
    % male births -> infinity => sex_ratio -> infinity => male_multiplier -> 1
    % female births -> 0 => sex_ratio -> infinity => female_multiplier -> 0
    % female births -> infinity => sex_ratio -> 0 => female_multiplier -> 1

    % Susceptible babies
    X(i_Susc, 1, i_female, i_undiagnosed) = female_multiplier * dt * babies_NotChronicCarriage;
    X(i_Susc, 1, i_male, i_undiagnosed)   = male_multiplier * dt * babies_NotChronicCarriage;

    % Babies with chronic carriage
    X(i_ImmTol, 1, i_female, i_undiagnosed) = female_multiplier * dt * babies_ChronicCarriage;
    X(i_ImmTol, 1, i_male, i_undiagnosed)   = male_multiplier * dt * babies_ChronicCarriage;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PAP Chunk 3:
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
    %% End of Chunk 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % increment the timestep index
    i_dt = i_dt + 1;
    % increases every dt (=0.1) years

    
end % end "time = TimeSteps" for loop

output.Time = Time; % 1 x (num_years_simul + 1)
output.Tot_Pop_1yr = Tot_Pop_1yr; % 2 x num_1yr_age_gps x (num_years_simul + 1)
output.num_births_1yr = num_births_1yr; % 1 x (num_years_simul + 1)
output.Incid_chronic_all_1yr_approx = Incid_chronic_all_1yr_approx; % 2 x num_1yr_age_gps x (num_years_simul + 1)
%%output.Prev_Immune_Reactive_1yr = Prev_Immune_Reactive_1yr; % 2 x num_1yr_age_gps x (num_years_simul + 1)
%%output.Prev_Chronic_Hep_B_1yr = Prev_Chronic_Hep_B_1yr; % 2 x num_1yr_age_gps x (num_years_simul + 1)
%%output.Prev_Comp_Cirr_1yr = Prev_Comp_Cirr_1yr; % 2 x num_1yr_age_gps x (num_years_simul + 1)
%%output.Prev_Decomp_Cirr_1yr = Prev_Decomp_Cirr_1yr; % 2 x num_1yr_age_gps x (num_years_simul + 1)
output.Prev_TDF_treat_1yr = Prev_TDF_treat_1yr; % 2 x num_1yr_age_gps x (num_years_simul + 1)
output.Prev_treatment_eligible_1yr = Prev_treatment_eligible_1yr; % 2 x num_1yr_age_gps x (num_years_simul + 1)
output.NumSAg_1yr = NumSAg_1yr; % 2 x num_1yr_age_gps x (num_years_simul + 1)
output.NumSAg_chronic_1yr = NumSAg_chronic_1yr; % 2 x num_1yr_age_gps x (num_years_simul + 1)
output.yld_1yr = yld_1yr; % 2 x num_1yr_age_gps x (num_years_simul + 1)
output.Incid_Deaths_1yr_approx = Incid_Deaths_1yr_approx; % 2 x num_1yr_age_gps x (num_years_simul + 1)
output.Prev_Deaths_1yr = Prev_Deaths_1yr; % 2 x num_1yr_age_gps x (num_years_simul + 1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PAP Chunk 4:
t_PAPoutputs_start = start_year;
t_PAPoutputs_end = end_year;

i_PAPoutputs_start = find(Time >= t_PAPoutputs_start, 1);
i_PAPoutputs_end   = find(Time >= t_PAPoutputs_end, 1);


output.PrevEAg = PrevEAg_of_SAg_5yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x 20 x num_years_output
%% Note that this excludes Vertical transmission
output.NewChronicInfectionRate = Incid_chronic_all_5yr_approx_no_VertTrans(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x 20 x num_years_output
%%output.Tot_Pop_1yr = Tot_Pop_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x num_1yr_age_gps x num_years_output
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
%%output.NumSAg_1yr = NumSAg_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x num_1yr_age_gps x num_years_output
output.NumEAg_chronic_1yr = NumEAg_chronic_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x num_1yr_age_gps x num_years_output
output.NumEAg_chronic_acute_1yr = NumEAg_chronic_acute_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x num_1yr_age_gps x num_years_output
%%output.yld_1yr = yld_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x num_1yr_age_gps x num_years_output
%%output.Incid_Deaths_1yr_approx = Incid_Deaths_1yr_approx(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x num_1yr_age_gps x num_years_output
%%output.Prev_Deaths_1yr = Prev_Deaths_1yr(:,:,i_PAPoutputs_start:i_PAPoutputs_end); % 2 x num_1yr_age_gps x num_years_output
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
output.num_starting_treatment_as_eligible = num_starting_treatment_as_eligible(i_PAPoutputs_start:i_PAPoutputs_end);

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
    'HBVPregnantWomenNeedToEvaluate',...
    'num_starting_treatment_as_eligible'...
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
    'Prev_TDF_treat_1yr',...
    'Prev_treatment_eligible_1yr'...
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
%% End of PAP Chunk 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DALYs = make_daly_mat(output,num_years_simul,num_year_1980_2100,life_expectancy);

DALYs_summed = sum(DALYs,1);

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
    %% The "4" below is used as follows - round( . , 4) rounds to 4 decimal places.
    writematrix(round([Time(i_1950:end);results_to_print(:,i_1950:end);DALYs_summed(i_1950:end)]',4),fullfile(basedir,'outputs',filename_results_csv));
end


end % end function HBVmodel_PPT

function output_labels=construct_header(agegroups, num_disease_states, num_sexes, num_treat_blocks)
    sex_labels = ["F","M"]; % F first in this model

    treat_labels = ["Undiagnosed","AppropriateManage","IncareNonadherent","OutOfCare"];
    assert(num_sexes==length(sex_labels))
    assert(num_treat_blocks==length(treat_labels))

    % Create labels for disease stage:
    D_labels = strings(1, num_disease_states); for i = 1:num_disease_states; D_labels(i) = "D" + string(i); end
    
    n_age_groups = max(agegroups);
    
    
    age_width = 1;  %% We are outputting in 1 year age groups.
    age_labels = strings(1, n_age_groups); 
    for i = 1:n_age_groups
        age_min = string((i-1)*age_width);
        age_max = string(i*age_width-1);
        age_labels(i) = "Age" + age_min + "_"+age_max;
    end

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
