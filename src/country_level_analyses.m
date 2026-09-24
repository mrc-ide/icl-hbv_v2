function country_level_analyses(sensitivity_analysis,...
    stochas_run_str,...
    ListOfISOs,countries_to_run,...
    BD_table,HepB3_table,...
    num_in_treatment_2016_map,pop_size_HBsAg_treatment_map,treatment_rates_map,...
    country_s_e_HCCdeaths_map,...
    params_map,stochas_params_mat,country_start_cols,...
    WUENIC2024BDdata, WUENIC2024HepB3data, ...
    Countrylevel_intervention_params, Global_intervention_params, ...
    GHO_infacilitybirthproportion_map, ANC_coverage_map, ...
    Polaris_diagnosis_coverage_map, Polaris_treat_coverage_map, ...
    basedir,i_natural_hist,i_sexes, i_care, ...
    num_year_divisions,dt,ages,num_age_steps,start_year,num_years_simul,end_year,...
    theta,CFR_Acute,rate_6months,ECofactor,p_ChronicCarriage,life_expectancy, Prog,...
    scenario_data_ANCHBVtestingbyage)

    if nargin < 1
       error('No input')
    end

    

    
    %% Intervention with HepB3 takes 3 years and can start now:
    T_INTERVENTION_START_HepB3 = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'Start_ContImp+HepB3'),:).Value;
    T_INTERVENTION_END_HepB3 = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'End_ContImp+HepB3'),:).Value;
    HepB3_WHO_target_coverage = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'HepB3_WHO_target_coverage'),:).Value;
    assert(T_INTERVENTION_START_HepB3>2024 && T_INTERVENTION_END_HepB3>2024)
    assert(HepB3_WHO_target_coverage>=0 && HepB3_WHO_target_coverage<=1)

    T_INTERVENTION_START = 2026.0;
    T_INTERVENTION_END = 2029.0;
    %% Birth dose takes 5 years:
    T_INTERVENTION_START_BD = 2026.0;
    T_INTERVENTION_END_BD  = 2031.0;

    RRprogress_effective_TDFtreatment_nonCC = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRprogress_effective_TDFtreatment_nonCC'),:).Value;
    RRprogress_effective_TDFtreatment_CC = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRprogress_effective_TDFtreatment_CC'),:).Value;
    RRprogress_nonadherent_TDFtreatment_nonCC = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRprogress_nonadherent_TDFtreatment_nonCC'),:).Value;
    RRprogress_nonadherent_TDFtreatment_CC = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRprogress_nonadherent_TDFtreatment_CC'),:).Value;
    
    RRprogress_effective_LAtreatment_nonCC = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRprogress_effective_LAtreatment_nonCC'),:).Value;
    RRprogress_effective_LAtreatment_CC = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRprogress_effective_LAtreatment_CC'),:).Value;
    RRprogress_nonadherent_LAtreatment_nonCC = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRprogress_nonadherent_LAtreatment_nonCC'),:).Value;
    RRprogress_nonadherent_LAtreatment_CC = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'RRprogress_nonadherent_LAtreatment_CC'),:).Value;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Modelled scenarios:
    %% Firstly we list the scenarios, then we list the possible options for each prevention/treatment.
    %% A given scenario consists of specifying the option chosen for each prevention/treatment.
    
    i_scenario_SQ = 1;  %% Post-2025 all coverage (BD,HepB3,PAP,diagnosis,treatment) kept fixed at 2025 levles.
    i_scenario_ContImp = 2;         %% Continued improvement - HepB3 and BD remain at current levels. PAP increases to 5% (2026-2030) of HVL. Treatment and diagnosis increase annually by region-specific rate
    i_scenario_ContImp_plusB3 = 3;  %% ContImp + HepB3 to 90% (or current coverage if higher) increasing T_INTERVENTION_START_HepB3-T_INTERVENTION_END_HepB3
    i_scenario_ContImp_plusBD_IF = 4;  %% ContImp + in-facility (IF) BD optimization: BD to IF ceiling in countries with BD at present
    i_scenario_ContImp_plusBD_OOF = 5;  %% ContImp + BD increases to also reach to 60% of non-facility births (BD introduced in non-BD countries).
    i_scenario_ContImp_plusBD_IF_OOF = 6;  %% ContImp + combined IF/OOF BD.
    i_scenario_ContImp_plusPAP_HVL_targeted = 7;  %% ContImp + X% PAP coverage of high VL (X%=access to VL monitoring). 
    i_scenario_ContImp_plusPAP_PoC = 8;  %% ContImp + PoC test to determine PAP eligibility. Assume 90% sensitivity, and 90% specificity.
    i_scenario_ContImp_plusPAP_all = 9;  %% ContImp + PAP ANC-1 coverage of all pregnant women (HVL+LVL)
    i_scenario_ContImp_plusDx_ANCscreening = 10;  %% ContImp with additional diagnosis among pregnant women during ANC (ANC-1 capped at HIV screening %)
    i_scenario_ContImp_plusDx_BirthCohort = 11;  %% ContImp with additional diagnosis through birth cohort screening (those born within 5 years of country introduction of HepB3)
    i_scenario_ContImp_plusDx_IFscreening = 12; %% ContImp with diagnosis through in-facility screening (e.g. acute care). Use WHO 
    i_scenario_ContImp_plusDx_CommunityScreening = 13; %% ContImp with additional diagnosis through (realistic) community screening similar to PROLIFICA.
    i_scenario_ContImp_plusDx_CommunityScreening_perfect = 14; %% ContImp with additional diagnosis through perfect community screening.
    i_scenario_ContImp_plusDx_IntegratedServices = 15;  %% CompInt with integrated services (TBD?)
    i_scenario_ContImp_plusTx_PoCeligibility = 16; %% CompInt with PoC treatment eligibility 
    i_scenario_ContImp_plusTx_treatall = 17;
    i_scenario_ContImp_plusTx_LA = 18;
    i_scenario_ContImp_plusDecentralisedDxTx = 19;
    i_scenario_ContImp_plusTx_cure_Bepi = 20;
    i_scenario_ContImp_plusTx_cure_improved = 21;

    %% Possible options for BD: different BD trends, changes in how BD is introduced etc.
    I_BD_WUENIC2025 = 100;  %% Follow WUENIC2025 and after 2024 coverage remains at last (2024) value
    I_BD_contimp = 101;  %% Possible increase beyond WUENIC2025 at a slow rate (currently 0). Introduction of BD to GAVI-approved countries
    I_BD_IFexpansion = 102; %% BD increases to current in-facilility % (or current coverage if higher) in all countries with BD/GAVI-approved.
    I_BD_OOFexpansion = 103; %% I_BD_IFexpansion with additional intervention to reach BD_oof_coverage=60% of OOF births (e.g. CHW, MAP). BD introduced in countries without BD.
    I_BD_IF_OOFexpansion = 104; %% Combining IF expansion + OOF 

    %% Possible options for HepB3 (scenario_HepB3)
    I_HEPB3_WUENIC2025 = 201;  %% Fixed at WUENIC 2025
    I_HEPB3_contimp = 202;      %% Possible increase beyond WUENIC2025 at a slow rate (currently 0).
    I_HEPB3_WHOtarget = 203;        %% Increase to WHO target.
    
    %% Possible options for PAP (scenario_PAP variable): peripartum antiviral prophylaxis (PAP) treatment for HBsAg+ mothers (treat all, treat high VL etc).
    I_PAP_SQ = 301;      %% No PAP unless already present.
    I_PAP_HVL_contimp = 302; %% Increase PAP to 20% of HVL unless already present.
    I_PAP_HVL_targeted = 303; %% Increase PAP to X% of HVL unless already present (X% placeholder but will be based on VL testing availability).
    I_PAP_PoC = 304;      %% Eligibility based on PoC test with given sensitivity, specificity and (globally constant) coverage.
    I_PAP_all = 305; %% All eligible (no VL criteria). Uptake to national ANC-1 value.
    
    %% scenario_Treatment: governs how treatment happens:
    % Modified by treat-all, introduction of PoC HBcrAg, PoC ALT tests, 
    I_TREAT = struct('SQ', 401,...  %% Current treatment (Capped at most recent treatemnt data - currently Polaris 2025).
        'continuedimprovement', 402,...     % 'HBV: Immune Tolerant' : HBeAg+ with very high HBV DNA (>1e6IU/ml), normal ALT
        'IFscreening', 403,...
        'IntegratedServices', 404,...
        'PoCeligibility', 405,...
        'universal', 406,...
        'LA', 407,...
        'decentralised', 408,...
        'cureBepi', 409,...
        'curev2', 410);

    
    % TUTAJ:
    num_scenarios = 21;
    %start_scenario = 17;
    start_scenario = 10;

    %%assert(ismember(sensitivity_analysis,{'default','infant_100','treat_medium','treat_high'}))

    stochas_run_num = str2double(stochas_run_str);

    %%filename_results = ['results_countries_', sensitivity_analysis, '_stochastic_run_', stochas_run_str, '.mat'];


    %%begin_time_run_num = datetime('now');
    %disp(append('Run number ',stochas_run_str,' (of ',num2str(num_stochas_runs),') started at ',string(begin_time_run_num)));

    outMap = containers.Map; 
    % for this particle (stochas_run_num), outMap contains the num_scenario scenarios, each of which contains countryMap, which contains the model results (lastrun) for 110 countries
    %%label_array = cell(1,num_scenarios);
    scenario_hours_vec = repmat(duration(0,0,0),1,num_scenarios);

    

    for scenario_num = start_scenario:num_scenarios
    %%for scenario_num = [3,9,10]
        disp("Running scenario")
        disp(scenario_num)
        % Make a copy of "Prog" for the given scenario - we can change
        % Prog_scenario in this loop if needed.
        Prog_scenario = Prog;

        begin_time_scenario = datetime('now');
        %%scenario = ListOfScenarios{scenario_num};
        %disp(scenario)
        disp(append('The "',num2str(scenario_num),'" scenario (run number ',stochas_run_str,') started at ',string(begin_time_scenario)));
        %%diary off
        %%diary(fullfile(basedir,'outputs',filename_diaries))

        countryMap = containers.Map; 
        % for this particle (stochas_run_num) and scenario (scenario), countryMap contains the model results (lastrun) of each of the 110 countries in this scenario

        for country_num = countries_to_run
            ISO = ListOfISOs{country_num};

            treatment_boundaries_vec = treatment_rates_map(ISO);
            assert(isequal(size(treatment_boundaries_vec),[1 6])) 
            % treatment_boundaries_vec contains the following (NOTE - this is *not* a complete list, but only treatment_boundaries_vec([1 2 3 5]) are used in the code:
	        % (1) treatment year, 
            % (2) rate to keep number of people in treatment constant,
	        % (3) rate to have 40% of eligible people in treatment by 2030 (if the constant rate is less than the 40% rate)
            % (5) rate to have 80% of eligible people in treatment by 2030 (if the constant rate is less than the 80% rate)
            treatment_start_year = treatment_boundaries_vec(1);
            assert(treatment_start_year==2016)
            in_treatment_2016_CDA = num_in_treatment_2016_map(ISO);
            assert(in_treatment_2016_CDA>=0)
            if in_treatment_2016_CDA>0
                pop_size_HBsAg_treatment_2016_vec = pop_size_HBsAg_treatment_map(ISO);
                %% This is the proportion of people who are on treatment in 2016 
                %% (it is different to the *rate* of treatment initiation 
                %% treatment_rate_params.annual_increase_TxifDx_past = stochas_params_mat(stochas_run_num,country_start_col+7))
                HBsAg_treat_cov_all_ages = pop_size_HBsAg_treatment_2016_vec(4);
                assert(HBsAg_treat_cov_all_ages>0 && HBsAg_treat_cov_all_ages<1, "HBsAg_treat_cov_all_ages must be between 0 and 1")
                assert(in_treatment_2016_CDA==pop_size_HBsAg_treatment_2016_vec(5))
            else
                HBsAg_treat_cov_all_ages = 0;
            end

            scenario_data_ANCHBVtestingbyage_thiscountry = scenario_data_ANCHBVtestingbyage(scenario_data_ANCHBVtestingbyage.ISO==ISO,:);

            %% Load country-specific parameters from calibration:
            params = params_map(ISO);
            %params = rmfield(params,'Efficacy_BirthDoseVacc_HbEAg');
            %params = rmfield(params,'Efficacy_InfantVacc');


            %% The previous code (copied from the PAP_scripts github repo branch) was the following:
            %% pr_VerticalTransmission_HbSAgHighVL_PAP and pr_VerticalTransmission_HbEAgHighVL_PAP
            %% pr_VerticalTransmission_HbSAgHighVL_PAP takes values 0.05, 0.06, or 0.1-1.0 in steps of 0.1
            %%parameter_to_vary_values_vec = [0.05 0.06 0.1:0.1:1];
            %%parameter_to_vary_num_vals = length(parameter_to_vary_values_vec);
            %%other_parameter_to_vary = strrep(parameter_to_vary,'HbSAg','HbEAg');
            %%effparams.(parameter_to_vary) = parameter_to_vary_values_vec(param_val_num);
            %%effparams.(other_parameter_to_vary) = parameter_to_vary_values_vec(param_val_num);
           
            

            Prog_scenario(i_natural_hist.HCC, i_natural_hist.HBVdeath) = params.CancerDeathRate;  % HCC to HBV death.


            % Parameters
            assert(params.Efficacy_BirthDoseVacc_HbSAg==0.95)
            assert(params.p_VerticalTransmission_HbEAg_NoIntv==0.9)
            %% MP: moved to main_script.m %% params.dwvec = dwvec;

            %% Load country-specific (non-stochastic - ie same across runs) parameters:
            country_start_col = country_start_cols(strcmp(ISO,ListOfISOs));
            params.beta_U5 = stochas_params_mat(stochas_run_num,country_start_col);

            %% MP: moved to main_script.m %% params.SpeedUpELoss_F = 9.5;

            %% Load country-specific stochastic parameters
            params.SpeedUpELoss_Beta = stochas_params_mat(stochas_run_num,country_start_col+1);
            
            %% ***WARNING*** This overwrites the value (0.0762) initially assigned when we load from params_map.mat in main_script.m
            params.p_VerticalTransmission_HbSAg_NoIntv = stochas_params_mat(stochas_run_num,country_start_col+2);
            params.cancer_rate_coeff = stochas_params_mat(stochas_run_num,country_start_col+3);
            params.cirrh_rate_coeff = stochas_params_mat(stochas_run_num,country_start_col+4);
            %% MP: moved to main_script.m %% params.CancerRate_WomenCoFactor = 1;
            %% MP: moved to main_script.m %% params.CirrhosisRate_WomenCoFactor = 1;
            params.CancerRate_MenCoFactor = stochas_params_mat(stochas_run_num,country_start_col+5);
            params.CirrhosisRate_MenCoFactor = stochas_params_mat(stochas_run_num,country_start_col+6);

            params.Efficacy_BirthDoseVacc_HbEAg = stochas_params_mat(stochas_run_num,end-1);
            params.Efficacy_InfantVacc = stochas_params_mat(stochas_run_num,end);


            %% PAP: assign PAP/VL parameters:
            %% ***WARNING***: because we overwrite params.p_VerticalTransmission_HbSAg_NoIntv using stochas_params_mat,
            %% we need to call assign_PAP_VL_params() here, otherwise it uses the wrong value of p_VerticalTransmission_HbSAg_NoIntv.
            PAP_VL_params = assign_PAP_VL_params(params);  %% PAP_VL_params was previously called effparams in the PAP branch.


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Load HBsAg prevalence data. This is used to initialise HBV prevalence in StartPrev_byAgeGroups
            country_ref_data_struct = country_s_e_HCCdeaths_map(ISO);
            source_HBsAg = country_ref_data_struct.source_HBsAg;
            %% HBsAg_prevs_year_1 = country_ref_data_struct.HBsAg_prevs_year_1; %% MP: this is not used.
            params.country_HBsAg_prevalences_by_ages_mid_1_young_old = country_ref_data_struct.country_HBsAg_prevalences_by_ages_mid_1_young_old;
            if strcmp(source_HBsAg,'Cui')
                params.HBsAg_prevs_middle_year_1 = country_ref_data_struct.HBsAg_prevs_middle_year_1;
            else
                assert(~ismember('HBsAg_prevs_middle_year_1',fields(country_ref_data_struct)))
            end
            if strcmp(source_HBsAg,'WHO')
                params.country_HBsAg_prevalences_by_ages_prevacc_young_old = country_ref_data_struct.country_HBsAg_prevalences_by_ages_prevacc_young_old;
            else
                assert(~ismember('country_HBsAg_prevalences_by_ages_prevacc_young_old',fields(country_ref_data_struct)))
            end




            % Summarise the Prog_scenario matrix as transactions lists.
            % Note - with the new treatment stratification, we have three
            % Transitions structures
            [transactions_from, transactions_to] = find(Prog_scenario > 0);

            % Load these into a data-structure (note that each transition element is stratified by same strata as state variable X):
            Transitions = [];
            %%Transitions = [];
            for tr = 1:length(transactions_from)

                %% This is the rate of transition per year (which we will then modify by treatment status):
                temparray = repmat(Prog_scenario(transactions_from(tr), transactions_to(tr)), [1, num_age_steps, i_sexes.n_sexes, i_care.n_care_blocks]);
                %%Transitions.Values(tr) = {temparray}; 
                
                % %% This is a vector containing the age indices people can be on treatment for the given starting natural history state:
                % treat_eligibility_ageindices = get_treatment_eligible_ageindices(ages, transactions_from(tr), i_natural_hist);
                % %% Different relative rate multiplier when on treatment (vs not) for CC/death versus other disease progressions:
                % if(transactions_to(tr)==i_natural_hist.HCC || transactions_to(tr)==i_natural_hist.HBVdeath)
                %     thisRR_effective_treatment = RRprogress_effective_treatment_CC;
                %     thisRR_nonadherent_treatment = RRprogress_nonadherent_treatment_CC;
                % elseif(transactions_from(tr)==i_natural_hist.ImmReact || transactions_from(tr)==i_natural_hist.Chronic ...
                %         || transactions_from(tr)==i_natural_hist.CompCirr || transactions_from(tr)==i_natural_hist.DecompCirr)
                %     %% These are the "non-CC" transitions:
                %     thisRR_effective_treatment = RRprogress_effective_treatment_nonCC;
                %     thisRR_nonadherent_treatment = RRprogress_nonadherent_treatment_nonCC;
                % else
                %     %% No effect (because not on treatment):
                %     thisRR_effective_treatment = 1;
                %     thisRR_nonadherent_treatment = 1;
                % end 
                % 
                % %% Now modify the rate of progression of temparray for any age groups that can be in treatment:
                % if(~isempty(treat_eligibility_ageindices))  %% Checks if any age groups can be in treatment for this natural history state
                %     temparray(:,treat_eligibility_ageindices,:,i_care.appropriate_management) = thisRR_effective_treatment*temparray(:,treat_eligibility_ageindices,:,i_care.appropriate_management);
                %     temparray(:,treat_eligibility_ageindices,:,i_care.incare_nonadherent) = thisRR_nonadherent_treatment*temparray(:,treat_eligibility_ageindices,:,i_care.incare_nonadherent);
                % end

                % For each to-from pair, form a 1 x num_age_steps x 2 x 2 double containing the progression parameter for that to-from transition
                % Arrange this sequence of matrices in a cell array called Transitions.Values, which is contained in Transitions
                Transitions.From(tr) = transactions_from(tr);
                Transitions.To(tr) = transactions_to(tr);
                Transitions.Values_withouttreat{tr} = temparray; 
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Now we add in age-varying natural history transitions:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Add age-specific progression acute severe/non-severe to immune tolerant or immune (14, 15)-->(2, 9)
            %% None of these vary by treatment (because not on treatment if acute)
            Transitions.From(length(Transitions.From) + 1) = i_natural_hist.NonSevAcute;
            Transitions.To(length(Transitions.To) + 1) = i_natural_hist.ImmTol;
            Transitions.Values_withouttreat{length(Transitions.Values_withouttreat) + 1} = p_ChronicCarriage * rate_6months;

            Transitions.From(length(Transitions.From) + 1) = i_natural_hist.NonSevAcute;
            Transitions.To(length(Transitions.To) + 1) = i_natural_hist.Immune;
            Transitions.Values_withouttreat{length(Transitions.Values_withouttreat) + 1} = (1 - p_ChronicCarriage) * rate_6months;

            Transitions.From(length(Transitions.From) + 1) = i_natural_hist.SevereAcute;
            Transitions.To(length(Transitions.To) + 1) = i_natural_hist.ImmTol;
            Transitions.Values_withouttreat{length(Transitions.Values_withouttreat) + 1} = p_ChronicCarriage * (1 - CFR_Acute) * rate_6months;
           
            Transitions.From(length(Transitions.From) + 1) = i_natural_hist.SevereAcute;
            Transitions.To(length(Transitions.To) + 1) = i_natural_hist.Immune;
            Transitions.Values_withouttreat{length(Transitions.Values_withouttreat) + 1} = (1 - p_ChronicCarriage) * (1 - CFR_Acute) * rate_6months;

            % Add age-specific progression immune tolerant -> immune reactive -> asymptomatic (2,3) and (3,4)
            AgeSpecELossFunction=params.SpeedUpELoss_F*exp(-params.SpeedUpELoss_Beta*ages);

            %% Immune tolerant -> Immune Reactive (age-specific)
            Transitions.From(length(Transitions.From)+1) = i_natural_hist.ImmTol;
            Transitions.To(length(Transitions.To)+1)     = i_natural_hist.ImmReact;
            %% Now make care-stratum specific (+age+sex) progression rate:
            temparray = repmat(Prog_scenario(i_natural_hist.ImmTol,i_natural_hist.ImmReact)*AgeSpecELossFunction,[1, 1, i_sexes.n_sexes, i_care.n_care_blocks]);
            Transitions.Values_withouttreat{length(Transitions.Values_withouttreat)+1} = temparray;

            %% Immune Reactive >- Asymptomatic (age-specific):
            Transitions.From(length(Transitions.From)+1) = i_natural_hist.ImmReact;
            Transitions.To(length(Transitions.To)+1) = i_natural_hist.AsymptCarr;
            %% Now make care-stratum specific (+age+sex) progression rate:
            temparray = repmat(Prog_scenario(i_natural_hist.ImmReact,i_natural_hist.AsymptCarr)*AgeSpecELossFunction,[1, 1, i_sexes.n_sexes, i_care.n_care_blocks]);
            Transitions.Values_withouttreat{length(Transitions.Values_withouttreat)+1} = temparray;

            % Modify immune reactive to chronic (3,5)  to an age-specific progression
            shortcut_rate_age = 20;
            i_shortcut_rate_age = find(ages >= shortcut_rate_age, 1); 

            %%indicator_vec = [zeros(1,shortcut_rate_age*10) ones(1,num_age_steps - shortcut_rate_age*10)];
            tmp_pos = find((Transitions.From==i_natural_hist.ImmReact) & (Transitions.To==i_natural_hist.Chronic));
            assert(isscalar(tmp_pos))
            tmp_mat = Transitions.Values_withouttreat{tmp_pos};
            %%tmp_mat = tmp_mat{1};
            assert(min(min(min(min(tmp_mat))))==max(max(max(max(tmp_mat)))))

            tmp_mat(:,1:i_shortcut_rate_age,:,:) = 0;
            %%Transitions.Values_withouttreat(tmp_pos) = {repmat(Prog_scenario(i_natural_hist.ImmReact, i_natural_hist.Chronic)*indicator_vec,...
            %%    [1, 1, i_sexes.n_sexes, i_care.n_care_blocks])}; % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
            Transitions.Values_withouttreat{tmp_pos} = tmp_mat;

            trans_rate_age = 25; % generic setting - at this age the rate of transition from chronic Hep B to comp cirrhosis is minimised (it's 0).
            trans_rate_by_age = (params.cirrh_rate_coeff*(ages - trans_rate_age)).^2; 
            trans_rate_by_age = trans_rate_by_age .* [zeros(1,trans_rate_age*10) ones(1,num_age_steps - trans_rate_age*10)];

            % Modify Chronic Hep B to Comp Cirrhosis (5,6) to an
            % age-specific progression:
            trans_rate_by_age_chronic_to_compcirr = Prog_scenario(i_natural_hist.Chronic,i_natural_hist.CompCirr)*trans_rate_by_age;
            tmp_pos = find((Transitions.From == i_natural_hist.Chronic) & (Transitions.To == i_natural_hist.CompCirr));
            assert(isscalar(tmp_pos))
            tmp_mat = Transitions.Values_withouttreat{tmp_pos};
            %%tmp_mat = tmp_mat{1};
            assert(min(min(min(min(tmp_mat))))==max(max(max(max(tmp_mat)))))
            tmp_mat = repmat(trans_rate_by_age_chronic_to_compcirr,[1, 1, i_sexes.n_sexes, i_care.n_care_blocks]); % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
            tmp_mat(:,:,i_sexes.F,:) = tmp_mat(:,:,i_sexes.F,:)*params.CirrhosisRate_WomenCoFactor;
            tmp_mat(:,:,i_sexes.M,:) = tmp_mat(:,:,i_sexes.M,:)*params.CirrhosisRate_MenCoFactor;
            %% Cap value at 5:
            tmp_mat = min(5, tmp_mat);
            Transitions.Values_withouttreat{tmp_pos} = tmp_mat;

            % Modify Immune Reactive to Comp Cirrhosis (3,6) to an age-specific progression
            trans_rate_by_age_immreact_to_compcirr = Prog_scenario(i_natural_hist.ImmReact, i_natural_hist.CompCirr)*trans_rate_by_age;
            tmp_pos = find((Transitions.From == i_natural_hist.ImmReact) & (Transitions.To == i_natural_hist.CompCirr));
            assert(isscalar(tmp_pos))
            tmp_mat = Transitions.Values_withouttreat{tmp_pos};
            %%tmp_mat = tmp_mat{1};
            assert(min(min(min(min(tmp_mat))))==max(max(max(max(tmp_mat)))))
            tmp_mat = repmat(trans_rate_by_age_immreact_to_compcirr,[1, 1, i_sexes.n_sexes, i_care.n_care_blocks]); % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
            tmp_mat(:,:,i_sexes.F,:) = tmp_mat(:,:,i_sexes.F,:)*params.CirrhosisRate_WomenCoFactor;
            tmp_mat(:,:,i_sexes.M,:) = tmp_mat(:,:,i_sexes.M,:)*params.CirrhosisRate_MenCoFactor;
            %% Again cap at max value 5:
            tmp_mat = min(5, tmp_mat);
            Transitions.Values_withouttreat{tmp_pos} = tmp_mat;


            % Add sex-specific co-factor to clearance (asympt -> immune) (4, 9)
            tmp_pos = find((Transitions.From == i_natural_hist.AsymptCarr) & (Transitions.To == i_natural_hist.Immune));
            assert(isscalar(tmp_pos))
            tmp = Transitions.Values_withouttreat{tmp_pos}; % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
            tmp(:, :, i_sexes.F, :) = tmp(:, :, i_sexes.F, :) * params.ClearanceRateWomenCoFactor;
            Transitions.Values_withouttreat{tmp_pos} = tmp;


            % Add age-specific progresion to HCC (2, 3, 4, 5, 6)-->8
            cancer_rate_age = 10;
            cancer_rate_by_age = (params.cancer_rate_coeff*(ages - cancer_rate_age)).^2; 
            %% Force cancer rate to be zero if age<(cancer_rate_age):
            cancer_rate_by_age = cancer_rate_by_age .* [zeros(1,cancer_rate_age*10) ones(1,num_age_steps - cancer_rate_age*10)];
            AgeSpecCancerFunction_Women = min(1, params.CancerRate_WomenCoFactor * cancer_rate_by_age);
            AgeSpecCancerFunction_Men = min(1, params.CancerRate_MenCoFactor * cancer_rate_by_age);

            AgeSpecificProgToCancer = zeros(1, num_age_steps, i_sexes.n_sexes, i_care.n_care_blocks);
            AgeSpecificProgToCancer(:, :, i_sexes.F, :) = repmat(AgeSpecCancerFunction_Women, [1, 1, 1, i_care.n_care_blocks]);
            AgeSpecificProgToCancer(:, :, i_sexes.M, :) = repmat(AgeSpecCancerFunction_Men,   [1, 1, 1, i_care.n_care_blocks]);

            Transitions.From(length(Transitions.From) + 1) = i_natural_hist.ImmTol;
            Transitions.To(length(Transitions.To) + 1) = i_natural_hist.HCC;
            Transitions.Values_withouttreat{length(Transitions.Values_withouttreat) + 1} = AgeSpecificProgToCancer;

            Transitions.From(length(Transitions.From) + 1) = i_natural_hist.ImmReact;
            Transitions.To(length(Transitions.To) + 1) = i_natural_hist.HCC;
            Transitions.Values_withouttreat{length(Transitions.Values_withouttreat) + 1} = 2 * AgeSpecificProgToCancer;

            Transitions.From(length(Transitions.From) + 1) = i_natural_hist.AsymptCarr;
            Transitions.To(length(Transitions.To) + 1) = i_natural_hist.HCC;
            Transitions.Values_withouttreat{length(Transitions.Values_withouttreat) + 1} = 0.5 * AgeSpecificProgToCancer;

            Transitions.From(length(Transitions.From) + 1) = i_natural_hist.Chronic;
            Transitions.To(length(Transitions.To) + 1) = i_natural_hist.HCC;
            Transitions.Values_withouttreat{length(Transitions.Values_withouttreat) + 1} = 2 * AgeSpecificProgToCancer;

            Transitions.From(length(Transitions.From) + 1) = i_natural_hist.CompCirr;
            Transitions.To(length(Transitions.To) + 1) = i_natural_hist.HCC;
            Transitions.Values_withouttreat{length(Transitions.Values_withouttreat) + 1} = 13 * AgeSpecificProgToCancer;


    

            years_vec_01yr = start_year:dt:end_year;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Load (WUENIC) data on BD and HepB3 coverage - we will use/modify these in the scenarios below.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % The code below loads WUENIC coverage data released in July
            % 2020 - this is what was used in the last round of fitting by
            % Margaret de Villers.
            % Extracts the years from the array HepB3_table (years where we have HepB3 coverage data from WUENIC 2020)
            yearsHepB3_coverage_available_wuenic2020 = cellfun(@(yyy) str2double(erase(yyy,"x")),HepB3_table.Properties.VariableNames);
            %InfantVacc_vec = HepB3_table; 
            %InfantVacc_vec = InfantVacc_vec(ISO,:);
            HepB3_wuenic2020 = HepB3_table(ISO,:);   %Extract the given country's BD coverage
            HepB3_wuenic2020 = HepB3_wuenic2020{:,:};  %Convert from a table to a vector
            assert(yearsHepB3_coverage_available_wuenic2020(1)==1980)

            % Extracts the years from the array BD_table (years where we have BD coverage data from WUENIC 2020)
            yearsBD_coverage_available_wuenic2020 = cellfun(@(yyy) str2double(erase(yyy,"x")),BD_table.Properties.VariableNames);
            %BirthDose_vec = BD_table;
            %BirthDose_vec = BirthDose_vec(ISO,:);
            BirthDose_wuenic2020 = BD_table(ISO,:);    %Extract the given country's (2020 WUENIC) BD coverage
            BirthDose_wuenic2020 = BirthDose_wuenic2020{:,:}; %Convert from a table to a vector
            
            last_available_coverage_data_year_wuenic2020 = 2019.0;
            index_last_available_year_WUENIC2020 = last_available_coverage_data_year_wuenic2020 - 1979;

            %% Check that the data contains the expected years (1980-2019):
            assert(isequal(size(HepB3_wuenic2020),[1 (2019-1979)]))
            assert(isequal(size(BirthDose_wuenic2020),[1 (2019-1979)]))
            assert(all(HepB3_wuenic2020>=0) && all(HepB3_wuenic2020<=1), "WUENIC 2020 HepB3 coverage must be between 0 and 1")
            assert(all(BirthDose_wuenic2020>=0) && all(BirthDose_wuenic2020<=1), "WUENIC 2020 BD coverage must be between 0 and 1")


            HepB3_wuenic2020 = HepB3_wuenic2020(1:index_last_available_year_WUENIC2020); % coverage of vaccination from 1980 to last_available_year
            BirthDose_wuenic2020 = BirthDose_wuenic2020(1:index_last_available_year_WUENIC2020); 

            assert(yearsBD_coverage_available_wuenic2020(1)==1980)
            assert(isequal(yearsBD_coverage_available_wuenic2020,yearsHepB3_coverage_available_wuenic2020))

            % 2019 is the last year of vaccination available from WUENIC2020

            %% Now load up the 2020-2024 BD and Hepb3 DATA
            %% For now we append the 2020-2024 values onto the 2019 WUENIC data (there are small differences in estimates up to 2019)
            %% So here we just extract the 2020-2024 data (2024 is the most recent year):
            if(ISO=="NIC")  %% Nicaragua currently not in WUENIC 2025 data!?
                BirthDose_wuenic2025 = [BirthDose_wuenic2020, ones(1,2024-last_available_coverage_data_year_wuenic2020)*BirthDose_wuenic2020(end)];
                HepB3_wuenic2025 = [HepB3_wuenic2020, ones(1,2024-last_available_coverage_data_year_wuenic2020)*HepB3_wuenic2020(end)];
                
            else
                firstyear_BD_coverage_available_wuenic2024 = min(table2array(WUENIC2024BDdata(strcmp(WUENIC2024BDdata.CODE,ISO), "YEAR")));
                firstyear_HepB3_coverage_available_wuenic2024 = min(table2array(WUENIC2024HepB3data(strcmp(WUENIC2024HepB3data.CODE,ISO), "YEAR")));
                
                %% Extract rows of the table corresponding to the country "ISO":
                BirthDose_wuenic2025_table_unsorted = WUENIC2024BDdata(strcmp(WUENIC2024BDdata.CODE,ISO), :);
                HepB3_wuenic2025_table_unsorted = WUENIC2024HepB3data(strcmp(WUENIC2024HepB3data.CODE,ISO), :);
                %% Sort the table by year (with year increasing)
                BirthDose_wuenic2025_table_sorted = sortrows(BirthDose_wuenic2025_table_unsorted, "YEAR");
                HepB3_wuenic2025_table_sorted = sortrows(HepB3_wuenic2025_table_unsorted, "YEAR");
                %% Now just pull out the coverage column, then convert from table to vector and divide by 100 (so proportion instead of percentage)
                BirthDose_wuenic2025 =  (table2array(BirthDose_wuenic2025_table_sorted(:,"COVERAGE"))')/100.0; 
                HepB3_wuenic2025 =  (table2array(HepB3_wuenic2025_table_sorted(:,"COVERAGE"))')/100.0; 
                %% Years where there is no coverage data are given as NA. Set these to be 0:
                BirthDose_wuenic2025(isnan(BirthDose_wuenic2025)) = 0;
                HepB3_wuenic2025(isnan(HepB3_wuenic2025)) = 0;

                BirthDose_wuenic2025 = [BirthDose_wuenic2020(1:(firstyear_BD_coverage_available_wuenic2024-yearsBD_coverage_available_wuenic2020(1))), BirthDose_wuenic2025];
                HepB3_wuenic2025 = [HepB3_wuenic2020(1:(firstyear_HepB3_coverage_available_wuenic2024-yearsHepB3_coverage_available_wuenic2020(1))), HepB3_wuenic2025];
                assert(isequal(length(BirthDose_wuenic2025),(length(BirthDose_wuenic2020)+5))) %% Check that the 2024 wuenic is now the same length as 2019 + 5
                assert(isequal(length(HepB3_wuenic2025),(length(HepB3_wuenic2020)+5))) %% Check that the 2024 wuenic is now the same length as 2019 + 5

                assert(all(HepB3_wuenic2025>=0) && all(HepB3_wuenic2025<=1), "WUENIC 2024 HepB3 coverage must be between 0 and 1")
                assert(all(BirthDose_wuenic2025>=0) && all(BirthDose_wuenic2025<=1), "WUENIC 2024 BD coverage must be between 0 and 1")

            end

       
            
            %% Get country-specific intervention data:
            Intervention_data_thiscountry = Countrylevel_intervention_params(strcmp(Countrylevel_intervention_params.ISO,ISO), :);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Set up scenarios:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% TUTAJ:
            
            

            
       
            %% For each scenario we determine which intervention "levers" are used.
            switch scenario_num
                case i_scenario_SQ     %% WUENIC 2025 BD+HepB3, 2016 treatment, no new interventions. Addresses - how have changes in BD+Hep B3 coverage impacted result?
                    scenario_BD = I_BD_WUENIC2025;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_SQ;
                    scenario_Treatment = I_TREAT.SQ;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp    %% Continued improvement - BD can increase (+ starts up in GAVI-approved countries). HepB3 can in crease
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusB3     %% Hep B3 increases to 90% 2026-2029 (T_INTERVENTION_START_HepB3-T_INTERVENTION_END_HepB3)
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_WHOtarget;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";                    
                case i_scenario_ContImp_plusBD_IF     %% BD increases - increasing OOF coverage
                    scenario_BD = I_BD_IFexpansion;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusBD_OOF     %% BD increases - increasing OOF coverage
                    scenario_BD = I_BD_OOFexpansion;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusBD_IF_OOF     %% BD increases - increasing IF+OOF coverage
                    scenario_BD = I_BD_IF_OOFexpansion;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusPAP_HVL_targeted     %% ContImp+ PAP for HVL only, capped at level of availability of VL testing.
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_targeted;
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusPAP_PoC     %% ContImp+ PAP eligibility through PoC test.
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_PoC;
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusPAP_all     %% ContImp+ PAP eligibility through PoC test.
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_all;  %% PAP available to all pregnant women regardless of VL
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";  
                
                %% ContImp with additional diagnosis among pregnant women during ANC (ANC-1 capped at HIV screening %)
                case i_scenario_ContImp_plusDx_ANCscreening
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "ANC screening";  
                 %% ContImp with additional diagnosis through birth cohort screening (those born within 5 years of country introduction of HepB3)
                case i_scenario_ContImp_plusDx_BirthCohort
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "Birth cohort screening";  
                case i_scenario_ContImp_plusDx_IFscreening
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.IFscreening;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% ContImp with additional diagnosis through (realistic) community screening similar to PROLIFICA.
                case i_scenario_ContImp_plusDx_CommunityScreening
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "Community screening";  
                %% ContImp with additional diagnosis through perfect community screening.
                case i_scenario_ContImp_plusDx_CommunityScreening_perfect
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.continuedimprovement;
                    scenario_AddScreenIntervention = "Perfect community screening";  
                %% CompInt with integrated services (TBD?)
                case i_scenario_ContImp_plusDx_IntegratedServices
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.IntegratedServices;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% CompInt with PoC treatment eligibility 
                case i_scenario_ContImp_plusTx_PoCeligibility
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.PoCeligibility;
                    scenario_AddScreenIntervention = "No additional screening";  
                case i_scenario_ContImp_plusTx_treatall
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.universal;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% Long-acting treatment available:
                case i_scenario_ContImp_plusTx_LA
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.LA;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% Decentralised testing and treatment:
                case i_scenario_ContImp_plusDecentralisedDxTx
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.decentralised;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% Bepi-like cure available:
                case i_scenario_ContImp_plusTx_cure_Bepi
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.cureBepi;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% Better-than-Bepi cure available:
                case i_scenario_ContImp_plusTx_cure_improved
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT.curev2;
                    scenario_AddScreenIntervention = "No additional screening";  
                otherwise
                    disp("Error - unknown scenario. Exiting")
                    return  %% Exit the script.
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% PARAMETERS FOR BD SCENARIOS:
            %% Here we specify what the maximum increase in BD coverage from MAP/CPAD would be:
            %%prop_accept_MAP = 0.9;   %% Placeholder assumption - easier to accept a patch than a needle
            %%prop_accept_CPAD = 0.85; %% Placeholder assumption

            %% In countries where BD is introduced in continued improvement scenario, this represents 
            %% Analysis is simply the median of WUENIC2024/(% in-facility births) for countries with BD - see BD.divided.by.infacility variable in process_polaris_timeseries_data.R
            BD_in_facility_acceptance_contimp = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'BD_in_facility_acceptance_contimp'),:).Value;

            %%prop_OOF_births = 1-GHO_infacilitybirthproportion_map(ISO);

            % Make model scenario birth dose coverage over time:
            %%last_BD_scaleup_year = 2030.0;  %% Assumption that BD plateaus after this time.
            
            switch scenario_BD
                case I_BD_WUENIC2025
                    disp("I_BD_WUENIC2025")
                    year_last_BD_data = 2024;
                    %% Update using WUENIC 2025: follow WUENIC2025 and after 2024 coverage remains at last (2024) value
                    coverage_BD_to_last_datapoint = BirthDose_wuenic2025;
                    future_xvals_vec = [2024.0, 2025.0, end_year];
                    future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), BirthDose_wuenic2025(end)];
                    %% No MAP or CPAD introduced:
                    scenario_BDcoverage_fromMAP = zeros(1,length(years_vec_01yr));
                    scenario_BDcoverage_fromCPAD = zeros(1,length(years_vec_01yr));
                case I_BD_contimp
                    year_last_BD_data = 2024;
                    disp("I_BD_contimp")
                    coverage_BD_to_last_datapoint = BirthDose_wuenic2025;
                    %% Non-GAVI eligible countries (as of 2026):
                    if(Intervention_data_thiscountry.HasBDorGAVIeligible==0)
                        future_xvals_vec = [2024.0, end_year];
                        future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end)];
                    else
                        %% BD increases at a slow rate in countries which already have BD; increases to BD_in_facility_acceptance_contimp% of in-facility births in 
                        %% Annual percentage point increase in Hep B3 (capped at 90%)
                        annual_BDimprovement = Intervention_data_thiscountry.ContImp_BD_annual_increase;
                        %% Check whether, at this rate of increase from 2025 onwards, if we ever go above in-facility births:
                        potential_bd_target = BirthDose_wuenic2025(end)+(end_year-2025)*annual_BDimprovement;
                        if(potential_bd_target>BirthDose_wuenic2025(end))
                            if(potential_bd_target>(GHO_infacilitybirthproportion_map(ISO)*BD_in_facility_acceptance_contimp))
                                year_reach_IF_target = floor(((GHO_infacilitybirthproportion_map(ISO)*BD_in_facility_acceptance_contimp)-BirthDose_wuenic2025(end))/annual_BDimprovement);
                                future_xvals_vec = [2024.0, 2025.0, year_reach_IF_target, end_year];
                                future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), GHO_infacilitybirthproportion_map(ISO), GHO_infacilitybirthproportion_map(ISO)];
                            else
                                future_xvals_vec = [2024.0, 2025.0, end_year];
                                future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), potential_bd_target];
                            end
                        else
                            future_xvals_vec = [2024.0, 2025.0, end_year];
                            future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), potential_bd_target];
                        end
                    end
                    %% No MAP or CPAD introduced:
                    scenario_BDcoverage_fromMAP = zeros(1,length(years_vec_01yr));
                    scenario_BDcoverage_fromCPAD = zeros(1,length(years_vec_01yr));
                case I_BD_IFexpansion  %% Optimisation of in-facility BD
                    year_last_BD_data = 2024;
                    disp("I_BD_IFexpansion")
                    coverage_BD_to_last_datapoint = BirthDose_wuenic2025;
                    %% Coverage up to % in-facility births, or current (2025 WUENIC) value - whichever is bigger.
                    max_in_facility_coverage = GHO_infacilitybirthproportion_map(ISO)*BD_in_facility_acceptance_contimp;
                    max_coverage = max(max_in_facility_coverage,BirthDose_wuenic2025(end));
                    assert(max_coverage<=1);
                    %% Currently 5 year scale-up of BD.
                    future_xvals_vec = [2024.0, T_INTERVENTION_START_BD, T_INTERVENTION_END_BD, end_year];
                    future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), max_coverage, max_coverage];
                    %% No MAP or CPAD introduced:
                    scenario_BDcoverage_fromMAP = zeros(1,length(years_vec_01yr));
                    scenario_BDcoverage_fromCPAD = zeros(1,length(years_vec_01yr));

                case I_BD_OOFexpansion
                    year_last_BD_data = 2024;
                    disp("I_BD_OOFexpansion")
                    coverage_BD_to_last_datapoint = BirthDose_wuenic2025;
                    %% Coverage up to % in-facility births, or current (2025 WUENIC) value - whichever is bigger.
                    max_in_facility_coverage = GHO_infacilitybirthproportion_map(ISO)*BD_in_facility_acceptance_contimp;

                    BD_oof_coverage = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'BD_oof_coverage'),:).Value;
                    max_OOF_coverage = (1.0-GHO_infacilitybirthproportion_map(ISO))*BD_oof_coverage;
                    
                    max_coverage = max(max_in_facility_coverage + max_OOF_coverage,BirthDose_wuenic2025(end));
                    assert(max_coverage<=1);
                    %% Currently 5 year scale-up of BD.
                    future_xvals_vec = [2024.0, T_INTERVENTION_START_BD, T_INTERVENTION_END_BD, end_year];
                    future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), max_coverage, max_coverage];
                    %% No MAP or CPAD introduced:
                    scenario_BDcoverage_fromMAP = zeros(1,length(years_vec_01yr));
                    scenario_BDcoverage_fromCPAD = zeros(1,length(years_vec_01yr));
                %% Combining in-facility and out-of-facility expansion:
                case I_BD_IF_OOFexpansion
                    year_last_BD_data = 2024;
                    disp("I_BD_OOFexpansion")
                    coverage_BD_to_last_datapoint = BirthDose_wuenic2025;
                    %% Coverage up to % in-facility births, or current (2025 WUENIC) value - whichever is bigger.
                    max_in_facility_coverage = GHO_infacilitybirthproportion_map(ISO)*BD_in_facility_acceptance_contimp;

                    BD_oof_coverage = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'BD_oof_coverage'),:).Value;
                    max_OOF_coverage = (1.0-GHO_infacilitybirthproportion_map(ISO))*BD_oof_coverage;
                    
                    max_coverage = max(max_in_facility_coverage + max_OOF_coverage,BirthDose_wuenic2025(end));
                    assert(max_coverage<=1);
                    %% Currently 5 year scale-up of BD.
                    future_xvals_vec = [2024.0, T_INTERVENTION_START_BD, T_INTERVENTION_END_BD, end_year];
                    future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), max_coverage, max_coverage];
                    %% No MAP or CPAD introduced:
                    scenario_BDcoverage_fromMAP = zeros(1,length(years_vec_01yr));
                    scenario_BDcoverage_fromCPAD = zeros(1,length(years_vec_01yr));

                    
                    % case I_BD_MAP  %% MAP introduced:
                %     disp("I_BD_MAP")
                %     year_last_BD_data = 2024;
                %     %% Follow WUENIC2025, then an extra (different efficacy) product increases overall BD coverage up to a level capped by out-of-facility deliveries.
                %     coverage_BD_to_last_datapoint = BirthDose_wuenic2025;
                %     %% This governs the coverage of the standard BD injection:
                %     future_xvals_vec = [2024.0, 2025.0, end_year];  
                %     future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), BirthDose_wuenic2025(end)];
                %     %% This governs the coverage of the additional MAP injection:
                %     future_xvals_vec_MAP = [2024.0, T_INTERVENTION_START, T_INTERVENTION_END, end_year];
                %     %% Increase in BD if introduce MAP (requires BD to currently be available):
                %     if(BirthDose_wuenic2025(end)>0)  %% BD already available
                %         %% Increase in BD is capped to not exceed the current proportion not getting BD.
                %         current_prop_not_getting_BD = (1-BirthDose_wuenic2025(end));
                %         BD_increase_from_MAP  = min(prop_OOF_births*prop_accept_MAP, current_prop_not_getting_BD);
                %     else                             %% BD not currently available
                %         BD_increase_from_MAP = 0;
                %     end
                % 
                %     future_yvals_vec_MAP = [0, 0, BD_increase_from_MAP, BD_increase_from_MAP];
                %     coverageMAP_to_present = zeros(1,length(BirthDose_wuenic2025));  
                %     disp("MAP1")
                %     scenario_BDcoverage_fromMAP = make_coverage_vec(start_year,num_year_divisions,dt,end_year,coverageMAP_to_present,future_xvals_vec_MAP, future_yvals_vec_MAP, year_last_BD_data);
                %     scenario_BDcoverage_fromCPAD = zeros(1,length(years_vec_01yr));
                 
                otherwise 
                    disp("Error - unknown scenario_BD. Exiting")
                    return  %% Exit the script.

            end  %% end switch scenario_BD
            
            %% Now get the full timetrend of BD coverage from start_year to end_year (note that MAP/CPAD coverage is stored separately in scenario_BDcoverage_fromMAP/CPAD)
            % disp([start_year,num_year_divisions,dt,end_year])
            % disp(coverage_BD_to_last_datapoint)
            % disp("A")
            % disp(future_xvals_vec)
            % disp("B")
            % disp(future_yvals_vec)
            % disp(year_last_BD_data)
            % disp("DONE")
            scenario_BDcoverage = make_coverage_vec(start_year,num_year_divisions,dt,end_year,coverage_BD_to_last_datapoint,future_xvals_vec,future_yvals_vec,year_last_BD_data);
            scenario_BDcoverage = min(1,scenario_BDcoverage);    % Ensure coverage is <=100% at every timestep:
            assert(isequal(size(scenario_BDcoverage),size(years_vec_01yr)))
            assert(isequal(size(scenario_BDcoverage_fromMAP),size(years_vec_01yr)))
            assert(isequal(size(scenario_BDcoverage_fromCPAD),size(years_vec_01yr)))
            %% Make sure scenario_BDcoverage_fromMAP/CPAD lies in range 0-1:
            assert(all(scenario_BDcoverage_fromMAP>=0) && all(scenario_BDcoverage_fromMAP<=1))
            assert(all(scenario_BDcoverage_fromCPAD>=0) && all(scenario_BDcoverage_fromCPAD<=1))
            %% Make sure total coverage of BD (normal BD/MAP/CPAD) is in range 0-1:
            %%mustBeBetween((scenario_BDcoverage+scenario_BDcoverage_fromMAP+scenario_BDcoverage_fromCPAD),0,1);
            assert(all((scenario_BDcoverage+scenario_BDcoverage_fromMAP+scenario_BDcoverage_fromCPAD)>=0) && all((scenario_BDcoverage+scenario_BDcoverage_fromMAP+scenario_BDcoverage_fromCPAD)<=1))

            %% Make sure that the total (normal BD + MAP/CPAD) BD coverage is 100% or less:
            %%scenario_BDcoverage_fromMAP_CPAD = min(scenario_BDcoverage_fromMAP_CPAD, 1 - scenario_BDcoverage);
            


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Now HepB3 coverage scenarios:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            switch scenario_HepB3
                case I_HEPB3_WUENIC2025
                    disp("I_HEPB3_WUENIC2025")
                    year_last_HepB3_data = 2024;
                    coverage_HepB3_to_last_datapoint = HepB3_wuenic2025;
                    future_xvals_vec = [2024.0, 2025.0, end_year];
                    future_yvals_vec = [HepB3_wuenic2025(end), HepB3_wuenic2025(end), HepB3_wuenic2025(end)];
                case I_HEPB3_contimp
                    disp("I_HEPB3_contimp")
                    year_last_HepB3_data = 2024;
                    coverage_HepB3_to_last_datapoint = HepB3_wuenic2025;
                    %% Annual percentage point increase in Hep B3 (capped at 90%)
                    annual_hepB3improvement = Intervention_data_thiscountry.ContImp_HepB3_annual_increase;
                    %% Check whether, at this rate of increase from 2025 onwards, if we ever go above the WHO target of 90%:
                    potential_hepb3_target = HepB3_wuenic2025(end)+(end_year-2025)*annual_hepB3improvement;
                    if(annual_hepB3improvement>0 && potential_hepb3_target>HepB3_WHO_target_coverage)
                        year_reach_whotarget = floor((HepB3_WHO_target_coverage-HepB3_wuenic2025(end))/annual_hepB3improvement);
                       
                        future_xvals_vec = [2024.0, 2025.0, year_reach_whotarget, end_year];
                        future_yvals_vec = [HepB3_wuenic2025(end), HepB3_wuenic2025(end), HepB3_WHO_target_coverage, HepB3_WHO_target_coverage];
                    else
                        future_xvals_vec = [2024.0, 2025.0, end_year];
                        future_yvals_vec = [HepB3_wuenic2025(end), HepB3_wuenic2025(end), potential_hepb3_target];
                    end
                case I_HEPB3_WHOtarget
                    disp("I_HEPB3_WHOtarget")
                    year_last_HepB3_data = 2024;
                    coverage_HepB3_to_last_datapoint = HepB3_wuenic2025;
                    %% Increase to HepB3_WHO_target_coverage (90%) (or current value if higher) from 2026 to 2029
                    hepb3_target = max(HepB3_wuenic2025(end),HepB3_WHO_target_coverage);
                    future_xvals_vec = [2024.0, T_INTERVENTION_START_HepB3 T_INTERVENTION_END_HepB3, end_year];
                    future_yvals_vec = [HepB3_wuenic2025(end), HepB3_wuenic2025(end), hepb3_target hepb3_target];
                otherwise
                    disp("Error: Unknown value for scenario_HepB3. Exiting")
                    return
                %if strcmp(sensitivity_analysis,'infant_100')
                % Ramp up coverage to 100% from first_expansion_year to
                % (first_expansion_year+0.1), and then keep it at 100%
                % until the end of the simulation.
                %future_xvals_vec = [2019.0 first_expansion_year (first_expansion_year+0.1) end_year];
                %future_yvals_vec = [HepB3_wuenic2020(end) HepB3_wuenic2020(end) 1 1];
                
            end  %% End switch scenario_HepB3

            % Using the above, create coverage vector that has coverage at each timestep:
          
            scenario_HepB3coverage = make_coverage_vec(start_year,num_year_divisions,dt,end_year,coverage_HepB3_to_last_datapoint,future_xvals_vec,future_yvals_vec,year_last_HepB3_data);
            
            % Ensure coverage is <=100% at every timestep:
            scenario_HepB3coverage = min(1,scenario_HepB3coverage);
            assert(isequal(size(scenario_HepB3coverage),size(years_vec_01yr)))
            %% Make sure scenario_HepB3coverage lies in range 0-1:
            assert(all(scenario_HepB3coverage>=0) && all(scenario_HepB3coverage<=1))
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% PAP coverage:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            PAP_cov_params = struct('current_cov_BDandPAP_EAgHighVL', 0,...
                    'current_cov_BDandPAP_SAgHighVL', 0,...
                    'current_cov_BDandPAP_EAgLowVL', 0,...
                    'current_cov_BDandPAP_SAgLowVL', 0,...
                    'current_cov_PAPonly_EAgHighVL', 0,...
                    'current_cov_PAPonly_SAgHighVL', 0,...
                    'current_cov_PAPonly_EAgLowVL', 0,...
                    'current_cov_PAPonly_SAgLowVL', 0,...
                    'max_cov_BDandPAP_EAgHighVL', 0,...  %% Coverage ceiling in future
                    'max_cov_BDandPAP_SAgHighVL', 0,...
                    'max_cov_BDandPAP_EAgLowVL', 0,...
                    'max_cov_BDandPAP_SAgLowVL', 0,...
                    'max_cov_PAPonly_EAgHighVL', 0,...
                    'max_cov_PAPonly_SAgHighVL', 0,...
                    'max_cov_PAPonly_EAgLowVL', 0,...
                    'max_cov_PAPonly_SAgLowVL', 0,...
                    'Past_TScaleup_PAP_start', 2024,... %% Dummy values:
                    'Past_TScaleup_PAP_end', 2025,...
                    'Intervention_TScaleup_PAP_start', 2026,...
                    'Intervention_TScaleup_PAP_end', 2030);
            %% This sets the past PAP coverage from Polaris estimates
            %% Dropbox_copy/Hepatits B/Data/Polaris/Polaris Database Query – CDA Foundation.xlsx
            switch ISO
                %%Bosnia - introduced 2018 at 43% overall
                case "BIH"
                    PAP_current_coverage = 0.43;
                    PAP_cov_params.Past_TScaleup_PAP_start = 2018;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2020;
                case "BLR" %% Belarus PAP introduced in 2019, 2025 coverage is 99%
                    PAP_current_coverage = 0.99;
                    PAP_cov_params.Past_TScaleup_PAP_start = 2019;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2021;                            
                case "CHN"
                    %% China 2019, 26% coverage.
                    PAP_current_coverage = 0.26;
                    PAP_cov_params.Past_TScaleup_PAP_start = 2019;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2021;
                case "CUB"
                    %% Cuba - first data in 2016, reached 100% in 2022.
                    PAP_current_coverage = 1.0;
                    PAP_cov_params.Past_TScaleup_PAP_start = 2016;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2022;                            
                case "ECU"
                    %% Ecuador introduced 2021, constant at 10% 2021-2025
                    PAP_current_coverage = 0.1;
                    PAP_cov_params.Past_TScaleup_PAP_start = 2021;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2023;
                case "EGY"
                    %% Egypt 45% in 2021 remaining roughly constant.
                    PAP_current_coverage = 0.45;
                    PAP_cov_params.Past_TScaleup_PAP_start = 2021;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2023;
                case "FSM"
                    % Micronesia: 50% since 2016
                    PAP_current_coverage = 0.50;
                    PAP_cov_params.Past_TScaleup_PAP_start = 2016;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2018;
                case "IRN"
                    % Iran: 10% in 2015, reaching 25% by 2021
                    PAP_current_coverage = 0.25;
                    PAP_cov_params.Past_TScaleup_PAP_start = 2015;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2021;
                case "KHM"
                    %% Cambodia 6% in 2022, reaching 23% by 2025
                    PAP_current_coverage = 0.23;
                    PAP_cov_params.Past_TScaleup_PAP_start = 2022;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2025;
                case "THA"
                    % Thailand 39% in 2020, reaching 41% by 2025
                    PAP_current_coverage = 0.41;
                    PAP_cov_params.Past_TScaleup_PAP_start = 2020;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2022;
                case "ZAF"
                    %% South Africa, 7% in 2015, reaching 11% in 2025.
                    PAP_current_coverage = 0.11;
                    PAP_cov_params.Past_TScaleup_PAP_start = 2015;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2019; 
                otherwise
                    %% No PAP currently available:
                    PAP_current_coverage = 0.0; 
                    PAP_cov_params.Past_TScaleup_PAP_start = 2019;
                    PAP_cov_params.Past_TScaleup_PAP_end = 2021;                            
            end
            
            PAP_cov_params.current_cov_BDandPAP_EAgHighVL = PAP_current_coverage;
            PAP_cov_params.current_cov_BDandPAP_SAgHighVL = PAP_current_coverage;
            PAP_cov_params.current_cov_BDandPAP_EAgLowVL  = 0;  
            PAP_cov_params.current_cov_BDandPAP_SAgLowVL  = 0;
            PAP_cov_params.current_cov_PAPonly_EAgHighVL = PAP_current_coverage;
            PAP_cov_params.current_cov_PAPonly_SAgHighVL = PAP_current_coverage;
            PAP_cov_params.current_cov_PAPonly_EAgLowVL  = 0;
            PAP_cov_params.current_cov_PAPonly_SAgLowVL  = 0;        

            %% Now set the future PAP coverage trend:
            PAP_cov_params.Intervention_TScaleup_PAP_start = 2026;
            PAP_cov_params.Intervention_TScaleup_PAP_end = 2030;
            switch scenario_PAP
                %% Scenarios where only HVL is eligible:
                case {I_PAP_SQ,I_PAP_HVL_contimp,I_PAP_HVL_targeted}
                    if (scenario_PAP==I_PAP_SQ)
                        PAP_max_coverage = PAP_current_coverage;
                    elseif (scenario_PAP==I_PAP_HVL_contimp)
                        %% 5% of HVL or current coverage (whichever is higher - this is done in the Excel file)
                        PAP_max_coverage = Intervention_data_thiscountry.ContImp_PAP_final_coverage;
                    elseif (scenario_PAP==I_PAP_HVL_targeted)
                        %% Currently a placeholder, but will be related to VL testing:
                        PAP_max_coverage = Intervention_data_thiscountry.PAP_HVL_targeted_coverage;
                    end
                    %% These scenarios all provide PAP only to HVL:
                    PAP_cov_params.max_cov_BDandPAP_EAgHighVL = PAP_max_coverage;
                    PAP_cov_params.max_cov_BDandPAP_SAgHighVL = PAP_max_coverage;
                    PAP_cov_params.max_cov_BDandPAP_EAgLowVL  = 0;  
                    PAP_cov_params.max_cov_BDandPAP_SAgLowVL  = 0;
                    PAP_cov_params.max_cov_PAPonly_EAgHighVL = PAP_max_coverage;
                    PAP_cov_params.max_cov_PAPonly_SAgHighVL = PAP_max_coverage;
                    PAP_cov_params.max_cov_PAPonly_EAgLowVL  = 0;
                    PAP_cov_params.max_cov_PAPonly_SAgLowVL  = 0;        

                case I_PAP_PoC
                    %% Coverage capped at ANC1
                    PoC_coverage = ANC_coverage_map(ISO);
                    PoC_sensitivity = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'PAP_PoC_sensitivity'),:).Value;
                    PoC_specificity = Global_intervention_params(strcmp(Global_intervention_params.Parameter,'PAP_PoC_specificity'),:).Value;

                    PAP_cov_params.max_cov_BDandPAP_EAgHighVL = PoC_coverage*PoC_sensitivity;
                    PAP_cov_params.max_cov_BDandPAP_SAgHighVL = PoC_coverage*PoC_sensitivity;
                    PAP_cov_params.max_cov_BDandPAP_EAgLowVL  = PoC_coverage*(1-PoC_specificity);  
                    PAP_cov_params.max_cov_BDandPAP_SAgLowVL  = PoC_coverage*(1-PoC_specificity);
                    PAP_cov_params.max_cov_PAPonly_EAgHighVL = PoC_coverage*PoC_sensitivity;
                    PAP_cov_params.max_cov_PAPonly_SAgHighVL = PoC_coverage*PoC_sensitivity;
                    PAP_cov_params.max_cov_PAPonly_EAgLowVL  = PoC_coverage*(1-PoC_specificity);
                    PAP_cov_params.max_cov_PAPonly_SAgLowVL  = PoC_coverage*(1-PoC_specificity);
                case I_PAP_all
                    %% Coverage up to ANC1
                    %% PAP is universal, capped at ANC-1 coverage.
                    ANC_coverage_level = ANC_coverage_map(ISO);
                    
                    PAP_cov_params.max_cov_BDandPAP_EAgHighVL = ANC_coverage_level; %% Everyone who wants PAP who got BD + is HVL gets it
                    PAP_cov_params.max_cov_BDandPAP_SAgHighVL = ANC_coverage_level;
                    PAP_cov_params.max_cov_BDandPAP_EAgLowVL  = ANC_coverage_level;  
                    PAP_cov_params.max_cov_BDandPAP_SAgLowVL  = ANC_coverage_level;
                    PAP_cov_params.max_cov_PAPonly_EAgHighVL = ANC_coverage_level;
                    PAP_cov_params.max_cov_PAPonly_SAgHighVL = ANC_coverage_level;
                    PAP_cov_params.max_cov_PAPonly_EAgLowVL  = ANC_coverage_level;
                    PAP_cov_params.max_cov_PAPonly_SAgLowVL  = ANC_coverage_level;        

                    %% In this code chunk, if ANC>scenario_BDcoverage(end) then we assume everyone who gets BD would also attend ANC.
                    % BD_coverage_level = scenario_BDcoverage(end);
                    % %% in_facility_births = GHO_infacilitybirthproportion_map(ISO);
                    % if(ANC_coverage_level>=BD_coverage_level)
                    %     PAP_coverage_withBD_ANCcap = 1.0; %% Everyone 
                    %     PAP_coverage_withoutBD_ANCcap = (ANC_coverage_level-BD_coverage_level)/(1.0-BD_coverage_level);
                    % else
                    %     %% Only those who got BD get PAP:
                    %     PAP_coverage_withBD_ANCcap = ANC_coverage_level/BD_coverage_level;
                    %     assert(PAP_coverage_withBD_ANCcap<=1)
                    %     PAP_coverage_withoutBD_ANCcap = 0; 
                    % end
                    
                otherwise
                    disp("Error: Unknown value for scenario_PAP. Exiting")
                    return
            end
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Now get time trends of PAP coverage (divided into those with/without BD, and by whether EAg+/SAg+ and high/low VL): 
            % Coverage of PAP among those with BD
            %% PAP code: cov_BirthDoseAndTDF_EAgHighVL_itt
            PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL = PAP_coverage_scaleup(start_year, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, end_year,...
                PAP_cov_params.current_cov_BDandPAP_EAgHighVL, PAP_cov_params.max_cov_BDandPAP_EAgHighVL, dt);

            PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL = PAP_coverage_scaleup(start_year, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, end_year,...
                PAP_cov_params.current_cov_BDandPAP_EAgLowVL, PAP_cov_params.max_cov_BDandPAP_EAgLowVL, dt);
            
            PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL = PAP_coverage_scaleup(start_year, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, end_year,...
                PAP_cov_params.current_cov_BDandPAP_SAgHighVL, PAP_cov_params.max_cov_BDandPAP_SAgHighVL, dt);
    
            PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL = PAP_coverage_scaleup(start_year, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, end_year,...
                PAP_cov_params.current_cov_BDandPAP_SAgLowVL, PAP_cov_params.max_cov_BDandPAP_SAgLowVL, dt);

            % Coverage of PAP among those not with BD
            PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgHighVL = PAP_coverage_scaleup(start_year, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, end_year,...
                PAP_cov_params.current_cov_PAPonly_EAgHighVL, PAP_cov_params.max_cov_PAPonly_EAgHighVL, dt);

            PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgLowVL = PAP_coverage_scaleup(start_year, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, end_year,...
                PAP_cov_params.current_cov_PAPonly_EAgLowVL, PAP_cov_params.max_cov_PAPonly_EAgLowVL, dt);
            
            PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgHighVL = PAP_coverage_scaleup(start_year, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, end_year,...
                PAP_cov_params.current_cov_PAPonly_SAgHighVL, PAP_cov_params.max_cov_PAPonly_SAgHighVL, dt);
    
            PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgLowVL = PAP_coverage_scaleup(start_year, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, end_year,...
                PAP_cov_params.current_cov_PAPonly_SAgLowVL, PAP_cov_params.max_cov_PAPonly_SAgLowVL, dt);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Diagnosis and treatment:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Dx_coverage_2016_thiscountry = 0.01; %% PLACEHOLDER
            annual_increase_Dx_past_thicountry = (Polaris_diagnosis_coverage_map(ISO) - Dx_coverage_2016_thiscountry) / (T_INTERVENTION_START-2016);
            %% Maybe change below to (*Note* - needs to be treat/diagnosis as we change Dx rates through interventions - this then filters through only if we use TxifDx rate)
            %%treatment_rate_params.annual_increase_TxifDx_past = (Polaris_treat_coverage_map(ISO) - HBsAg_treat_cov_all_ages) / (T_INTERVENTION_START-2016);

            %% For CHN:
            %%annual_increase_TxifDx_past = 0.0191;
            %%treatment_rate_params.Tx_coverage_2016 = 0.0413 (aka HBsAg_treat_cov_all_ages)
            %%Polaris_treat_coverage_map("CHN") = 0.30;
            %%Polaris_diagnosis_coverage_map("CHN") = 0.68; 
            treatment_rate_params = struct('Tx_coverage_2016', HBsAg_treat_cov_all_ages, ...
                    'Dx_coverage_2016', Dx_coverage_2016_thiscountry, ...
                    'annual_increase_TxifDx_past', stochas_params_mat(stochas_run_num,country_start_col+7),...
                    'annual_increase_Dx_past',annual_increase_Dx_past_thicountry,...        
                    't_treatment_scaleup_start', T_INTERVENTION_START,...  %% Treatment takes a few years to change from current rate of increase to new one.
                    't_treatment_scaleup_end', T_INTERVENTION_END,...
                    'annual_increase_Dx_future',0,...
                    'annual_increase_TxifDx_future',0);

            %% annual_increase_TxifDx_past is the annual rate of treatment increase from 2016 to current time
            switch scenario_Treatment
                case I_TREAT.SQ           %% Use current rates of treatment uptake and failure.
                    scenario_treat_elig = "Current treatment";                    
                    treatment_rate_params.annual_increase_Dx_future = 0;
                    treatment_rate_params.annual_increase_TxifDx_future = 0;                    
                case I_TREAT.continuedimprovement 
                    scenario_treat_elig = "Current treatment";
                    treatment_rate_params.annual_increase_Dx_future = Intervention_data_thiscountry.ContImp_Dx_annual_increase;
                    treatment_rate_params.annual_increase_TxifDx_future = Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                case I_TREAT.IFscreening
                    scenario_treat_elig = "Current treatment";
                    %% PLACEHOLDER - use data by country 
                    %prop_accessing_healthcare_F = [0,0,0,0.02,0.02,0.02,0.03,0.03,0.03,0.04,0.04,0.05,0.1,0.2,0.4,0.4,0.5,0.5,0.5,0.5];
                    %prop_accessing_healthcare_M = [0,0,0,0.02,0.02,0.02,0.03,0.03,0.03,0.04,0.04,0.05,0.1,0.2,0.4,0.4,0.5,0.5,0.5,0.5];
                    prop_accessing_healthcare_and_accepttest = Intervention_data_thiscountry.Dx_coverage_with_infacilitytesting;
                    
                    treatment_rate_params.annual_increase_Dx_future = Intervention_data_thiscountry.ContImp_Dx_annual_increase + prop_accessing_healthcare_and_accepttest;
                    treatment_rate_params.annual_increase_TxifDx_future = Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                %% PLACEHOLDER - does nothing:
                case I_TREAT.IntegratedServices
                    scenario_treat_elig = "Current treatment";
                    treatment_rate_params.annual_increase_Dx_future = Intervention_data_thiscountry.ContImp_Dx_annual_increase;
                    treatment_rate_params.annual_increase_TxifDx_future = Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                %% PLACEHOLDER - does nothing:
                case I_TREAT.PoCeligibility
                    scenario_treat_elig = "Current treatment";
                    treatment_rate_params.annual_increase_Dx_future = Intervention_data_thiscountry.ContImp_Dx_annual_increase;
                    treatment_rate_params.annual_increase_TxifDx_future = Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                %% PLACEHOLDER:
                case I_TREAT.universal
                    scenario_treat_elig = "Universal treatment";
                    treatment_rate_params.annual_increase_Dx_future = Intervention_data_thiscountry.ContImp_Dx_annual_increase;
                    treatment_rate_params.annual_increase_TxifDx_future = Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                %% PLACEHOLDER:
                case I_TREAT.LA
                    scenario_treat_elig = "Current treatment";
                    RR_LA = 1.0; %% PLACEHOLDER - Arbitrary increase in TxifDx if needed
                    treatment_rate_params.annual_increase_Dx_future = Intervention_data_thiscountry.ContImp_Dx_annual_increase;
                    treatment_rate_params.annual_increase_TxifDx_future = RR_LA*Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                case I_TREAT.decentralised
                    scenario_treat_elig = "Current treatment";
                    %% TO DO - DECIDE IF THERE SHOULD BE AN (IMMEDIATE?) INCREASE IN THOSE CURRENTLY ON TREATMENT BY THE SAME MULTIPLIER.
                    %% Proportion of country that is rural:
                    prop_rural = Intervention_data_thiscountry.Decentralisation_prop_rural;
                    increase_by_decentralisation = 1 + prop_rural/(1-prop_rural);
                    treatment_rate_params.annual_increase_Dx_future = increase_by_decentralisation*Intervention_data_thiscountry.ContImp_Dx_annual_increase;
                    treatment_rate_params.annual_increase_TxifDx_future = increase_by_decentralisation*Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                %% PLACEHOLDER
                case I_TREAT.cureBepi
                    scenario_treat_elig = "Current treatment";
                    treatment_rate_params.annual_increase_Dx_future = Intervention_data_thiscountry.ContImp_Dx_annual_increase;
                    treatment_rate_params.annual_increase_TxifDx_future = Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                %% PLACEHOLDER
                case I_TREAT.curev2
                    scenario_treat_elig = "Current treatment";
                    treatment_rate_params.annual_increase_Dx_future = Intervention_data_thiscountry.ContImp_Dx_annual_increase;
                    treatment_rate_params.annual_increase_TxifDx_future = Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                % case I_diag70percent
                %     scenario_treat_elig = "Current treatment";
                %     treatment_rate_params.prop_diagnosed_now = Polaris_diagnosis_coverage_map(ISO);
                %     treatment_rate_params.prop_treat_now = Polaris_treat_coverage_map(ISO);
                %     treatment_rate_params.annual_increase_Dx_future = max(0, (0.70-Polaris_diagnosis_coverage_map(ISO))/(T_INTERVENTION_END - T_INTERVENTION_START));
                %     treatment_rate_params.annual_increase_TxifDx_future = 0;

                otherwise
                    disp("Error: Unknown value for scenario_Treatment. Exiting")
                    return
            end

            %% This is now dealt with in HBVmodel.m:
            % switch scenario_AddScreenIntervention
            %     case I_NO_ADDITIONAL_SCREENING
            %         a=7;
            %     case I_BIRTHCOHORT_SCREENING
            %         a=8;
            %     otherwise
            %         disp("Error: Unknown value for scenario_AddScreenIntervention. Exiting")
            %         return
            % end  

            % switch sensitivity_analysis
            %     case {'default','infant_100'}
            %         params.PriorTDFTreatRate = stochas_params_mat(stochas_run_num,country_start_col+7);
            %         assert((params.PriorTDFTreatRate>=treatment_boundaries_vec(2)) && (params.PriorTDFTreatRate<=treatment_boundaries_vec(3)))
            %     case 'treat_medium'
            %         params.PriorTDFTreatRate = treatment_boundaries_vec(3); % 40%
            %     case 'treat_high'
            %         params.PriorTDFTreatRate = treatment_boundaries_vec(5); % 80%
            %     % otherwise
            %     %     disp("Error: Unknown value for sensitivity_analysis. Exiting")
            %     %     return
            % end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Now alter Transitions to account for treatment - this is the least bad way I can see to do this.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% RRprogress are the relative rates when on treatment (effective or non-adherent) at baseline (i.e. prior to any intervention, pre-2026) 
            %% and "future" - future could include long-acting treatment (which is dealt with here) or changes to guidelines (dealt with later in get_treatment_eligible_ageindices()).
            RRprogress_effective_treatment_nonCC_baseline = RRprogress_effective_TDFtreatment_nonCC;
            RRprogress_effective_treatment_CC_baseline = RRprogress_effective_TDFtreatment_CC;
            RRprogress_nonadherent_treatment_nonCC_baseline = RRprogress_nonadherent_TDFtreatment_nonCC;
            RRprogress_nonadherent_treatment_CC_baseline = RRprogress_nonadherent_TDFtreatment_CC;

            %% Now post-present (with changes to guidelines and treatments) relative rates:
            if(scenario_Treatment==I_TREAT.SQ || scenario_Treatment==I_TREAT.continuedimprovement || ...
                    scenario_Treatment==I_TREAT.IFscreening || scenario_Treatment==I_TREAT.IntegratedServices || ...
                    scenario_Treatment==I_TREAT.PoCeligibility || scenario_Treatment==I_TREAT.universal || ...
                    scenario_Treatment==I_TREAT.decentralised || ...
                    scenario_Treatment==I_TREAT.cureBepi || scenario_Treatment==I_TREAT.curev2)
                %% Tenofovir-based treatment:
                RRprogress_effective_treatment_nonCC_future = RRprogress_effective_TDFtreatment_nonCC;
                RRprogress_effective_treatment_CC_future = RRprogress_effective_TDFtreatment_CC;
                RRprogress_nonadherent_treatment_nonCC_future = RRprogress_nonadherent_TDFtreatment_nonCC;
                RRprogress_nonadherent_treatment_CC_future = RRprogress_nonadherent_TDFtreatment_CC;
            elseif(scenario_Treatment==I_TREAT.LA)
                %% Long-acting treatment:
                RRprogress_effective_treatment_nonCC_future = RRprogress_effective_LAtreatment_nonCC;
                RRprogress_effective_treatment_CC_future = RRprogress_effective_LAtreatment_CC;
                RRprogress_nonadherent_treatment_nonCC_future = RRprogress_nonadherent_LAtreatment_nonCC;
                RRprogress_nonadherent_treatment_CC_future = RRprogress_nonadherent_LAtreatment_CC;
            else
                disp("Error - unknown treatment scenario. Exiting")
                return
            end

            for i_transition = 1:length(Transitions.From)

                start_state = Transitions.From(i_transition);
                end_state = Transitions.To(i_transition);
                
                %% This sets the scenario's eligibility guidelines for baseline and future:
                treat_eligibility_ageindices_baseline = get_treatment_eligible_ageindices("Current treatment", start_state, i_natural_hist, ages);
                treat_eligibility_ageindices_future = get_treatment_eligible_ageindices(scenario_treat_elig, start_state, i_natural_hist, ages);
                
                %% Now modify the rate of progression of temparray for any age groups that can be in treatment:
                
                %% We use different relative rate multiplier when on treatment (vs not) for HCC/death versus other disease progressions:
                if(~isempty(treat_eligibility_ageindices_baseline))  %% Checks if any age groups can be in treatment for this natural history state
                    %% Different relative rate multiplier when on treatment (vs not) for HCC/death versus other disease progressions:
                    if(end_state==i_natural_hist.HCC || end_state==i_natural_hist.HBVdeath)
                        thisRR_effective_treatment_baseline = RRprogress_effective_treatment_CC_baseline;
                        thisRR_nonadherent_treatment_baseline = RRprogress_nonadherent_treatment_CC_baseline;
                    elseif(start_state==i_natural_hist.ImmReact || start_state==i_natural_hist.Chronic ...
                            || start_state==i_natural_hist.CompCirr || start_state==i_natural_hist.DecompCirr)
                        %% These are the "non-CC" transitions:
                        thisRR_effective_treatment_baseline = RRprogress_effective_treatment_nonCC_baseline;
                        thisRR_nonadherent_treatment_baseline = RRprogress_nonadherent_treatment_nonCC_baseline;
                    else
                        %% No effect (because not on treatment):
                        thisRR_effective_treatment_baseline = 1;
                        thisRR_nonadherent_treatment_baseline = 1;
                    end 
                end
                %% Now repeat for future (i.e. when we may have different treatment guidelines/treatment types):
                if(~isempty(treat_eligibility_ageindices_future))  %% Checks if any age groups can be in treatment for this natural history state
                    %% Different relative rate multiplier when on treatment (vs not) for HCC/death versus other disease progressions:
                    if(end_state==i_natural_hist.HCC || end_state==i_natural_hist.HBVdeath)
                        thisRR_effective_treatment_future = RRprogress_effective_treatment_CC_future;
                        thisRR_nonadherent_treatment_future = RRprogress_nonadherent_treatment_CC_future;
                    elseif(start_state==i_natural_hist.ImmReact || start_state==i_natural_hist.Chronic ...
                            || start_state==i_natural_hist.CompCirr || start_state==i_natural_hist.DecompCirr)
                        %% These are the "non-CC" transitions:
                        thisRR_effective_treatment_future = RRprogress_effective_treatment_nonCC_future;
                        thisRR_nonadherent_treatment_future = RRprogress_nonadherent_treatment_nonCC_future;
                    else
                        %% No effect (because not on treatment):
                        thisRR_effective_treatment_future = 1;
                        thisRR_nonadherent_treatment_future = 1;
                    end 
                end
                %% Transitions.Values_withouttreat{i_transition} is the transition matrix in the absence of treatment.
                %% As the transition matrix is age-dependent, we need to make a pre-intervention ("baseline") transition matrix and a post-intervention ("future") one. 
                %% Specifically for "Universal treatment" - the reason is that we need the current guidelines (up to 2026 in the simulation), and then we need the new guidelines.
                temparray_baseline = Transitions.Values_withouttreat{i_transition}; %% Pre-2026
                temparray_future = Transitions.Values_withouttreat{i_transition};   %% Post-2026
                    
                temparray_baseline(:,treat_eligibility_ageindices_baseline,:,i_care.appropriate_management) = ... 
                    thisRR_effective_treatment_baseline*temparray_baseline(:,treat_eligibility_ageindices_baseline,:,i_care.appropriate_management);
                temparray_baseline(:,treat_eligibility_ageindices_baseline,:,i_care.incare_nonadherent) = ... 
                    thisRR_nonadherent_treatment_baseline*temparray_baseline(:,treat_eligibility_ageindices_baseline,:,i_care.incare_nonadherent);

                temparray_future(:,treat_eligibility_ageindices_future,:,i_care.appropriate_management) = ... 
                    thisRR_effective_treatment_future*temparray_future(:,treat_eligibility_ageindices_future,:,i_care.appropriate_management);
                temparray_future(:,treat_eligibility_ageindices_future,:,i_care.incare_nonadherent) = ... 
                    thisRR_nonadherent_treatment_future*temparray_future(:,treat_eligibility_ageindices_future,:,i_care.incare_nonadherent);

                %% Store the updated matrices (for baseline and future):
                Transitions.Values_baseline{i_transition} = temparray_baseline;
                Transitions.Values_future{i_transition} = temparray_future;
                %disp(append('i_transition1 =',num2str(i_transition),' ',num2str(size(Transitions.Values{i_transition}))))
                
            end  %% End for loop.




            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Now call HBVmodel.m:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            num_year_1980_2100 = 2100 - 1980 + 1;

            %% MP: added so we just store the CSV from the default scenario for now to save on storage. Can be changed as needed.
            if(strcmp(sensitivity_analysis,'default')==1)
                store_results_as_text = 1;
            else
                store_results_as_text = 0;
            end
            %% Run scenarios:
            lastrun = HBVmodel(source_HBsAg,...
                num_year_divisions,dt,ages,num_age_steps,i_natural_hist,i_sexes,i_care,...
                start_year,num_years_simul,...
                theta,ECofactor,treatment_rate_params, treatment_start_year-dt, ...
                params, PAP_VL_params, PAP_cov_params, ...
                Global_intervention_params, Intervention_data_thiscountry, ...
                p_ChronicCarriage,Prog_scenario,Transitions,......
                scenario_BDcoverage, scenario_BDcoverage_fromMAP,...
                scenario_BDcoverage_fromCPAD, scenario_HepB3coverage, ...
                scenario_Treatment, I_TREAT,...
                scenario_treat_elig, scenario_data_ANCHBVtestingbyage_thiscountry, ...
                ISO, scenario_num, scenario_AddScreenIntervention, ...
                num_year_1980_2100, life_expectancy, ...
                stochas_run_str, sensitivity_analysis, basedir, store_results_as_text);
            lastrun.DALYPerYear = make_daly_mat(lastrun,num_years_simul,num_year_1980_2100,life_expectancy);
            assert(isequal(size(lastrun.DALYPerYear),[100 num_years_simul + 1]))

            assert(isequal(size(lastrun.Time),[1 (num_years_simul + 1)]))
            i1980 = find(lastrun.Time>=1980, 1);
            i2100 = find(lastrun.Time>=2100, 1);
            num_cols_out = i2100 - i1980 + 1; % every output should have entries for the years 1980 to 2100
            i5y = 6;
            index_under_5y = 1:5;
            num_cols_in = i2100 - i1980 + 1;
            index_last_year = i2100;
            
            lastrun = rmfield(lastrun,'Prev_Deaths_1yr');
            lastrun = rmfield(lastrun,'yld_1yr');
            lastrun_fields = sort(fields(lastrun));
            fields_of_interest = {'Time',...
                'Tot_Pop_1yr','num_births_1yr',...
                'Incid_chronic_all_1yr_approx',...
                'NumSAg_1yr','NumSAg_chronic_1yr',...
                ...%'Prev_TDF_treat_1yr','Prev_Immune_Reactive_1yr','Prev_Chronic_Hep_B_1yr','Prev_Comp_Cirr_1yr','Prev_Decomp_Cirr_1yr',...
                'Prev_TDF_treat_1yr','Prev_treatment_eligible_1yr',...
                'Incid_Deaths_1yr_approx',...
                'DALYPerYear',...
                'HBVPregnantWomenNeedToEvaluate', 'NewChronicInfectionRate', 'NewChronicInfectionRate_NeonatesOnly',...
                'NumDecompCirr', 'NumEAg_chronic_1yr', 'NumEAg_chronic_acute_1yr', 'NumLiverCancer', ...
                'PeripartumTreatment_HbEAg_HighVL_approx', 'PeripartumTreatment_HbEAg_LowVL_approx', ...
                'PeripartumTreatment_HbSAg_HighVL_approx', 'PeripartumTreatment_HbSAg_LowVL_approx', ...
                'PregnantWomenNeedToScreen', 'PrevEAg', 'RateBirthDoseVacc', 'RateInfantVacc', 'RatePeripartumTreatment', ...
                'beta_U5', 'num_births_1yr_approx', 'num_births_chronic_HbEAgWomenHVL_1yr_approx', ...
                'num_births_chronic_HbEAgWomenLVL_1yr_approx', 'num_births_chronic_HbSAgWomenHVL_1yr_approx', ...
                'num_births_chronic_HbSAgWomenLVL_1yr_approx', 'num_births_toHbEAgWomenHVL_1yr_approx', ...
                'num_births_toHbEAgWomenLVL_1yr_approx', 'num_births_toHbSAgWomenHVL_1yr_approx', ...
                'num_births_toHbSAgWomenLVL_1yr_approx', 'p_HbSAg_av', 'p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc', ...
                'p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP', 'p_VerticalTransmission_HbEAgHighVL_NoIntv', ...
                'p_VerticalTransmission_HbEAgHighVL_PAP', 'p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc', ...
                'p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP', 'p_VerticalTransmission_HbEAgLowVL_NoIntv', ...
                'p_VerticalTransmission_HbEAgLowVL_PAP', 'p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL', ...
                'p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc', 'p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP', ...
                'p_VerticalTransmission_HbSAgHighVL_NoIntv', 'p_VerticalTransmission_HbSAgHighVL_PAP', ...
                'p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc', 'p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP', ...
                'p_VerticalTransmission_HbSAgLowVL_NoIntv', 'p_VerticalTransmission_HbSAgLowVL_PAP', ...
                'p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL',...
                'num_starting_treatment_as_eligible'...
                };
            assert(all(ismember(fields_of_interest,lastrun_fields)))
            assert(all(ismember(lastrun_fields,fields_of_interest)))

            % lastrun.Prev_treatment_eligible_1yr = ...
            %     lastrun.Prev_Immune_Reactive_1yr + lastrun.Prev_Chronic_Hep_B_1yr + lastrun.Prev_Comp_Cirr_1yr + lastrun.Prev_Decomp_Cirr_1yr;
                %%lastrun.Prev_TDF_treat_1yr;
            % lastrun = rmfield(lastrun,'Prev_Immune_Reactive_1yr');
            % lastrun = rmfield(lastrun,'Prev_Chronic_Hep_B_1yr');
            % lastrun = rmfield(lastrun,'Prev_Comp_Cirr_1yr');
            % lastrun = rmfield(lastrun,'Prev_Decomp_Cirr_1yr');
            %% MAGIC NUMBERS - 2 (NUMBER OF SEXES?) AND 100 (NUMBER OF 1-YEAR AGE GROUPS?)
            assert(isequal(size(lastrun.Incid_Deaths_1yr_approx),[2 100 (num_years_simul + 1)]))
            birth_cohorts_fun = @(yy) arrayfun(@(xx) sum(diag(yy,xx)), 0:(num_cols_out-1));
            assert(isequal(size(squeeze(sum(lastrun.Incid_Deaths_1yr_approx(:,:,i1980:index_last_year),1))),[100 num_cols_in]))
            lastrun.Incid_Deaths_1yr_approx_birth_cohorts = birth_cohorts_fun(squeeze(sum(lastrun.Incid_Deaths_1yr_approx(:,:,i1980:index_last_year),1))); % 2 x 100 x (num_years_simul + 1)
            assert(isequal(size(lastrun.Incid_Deaths_1yr_approx_birth_cohorts),[1 num_cols_out]))
            lastrun.Time = lastrun.Time(i1980:i2100); % 1 x (num_years_simul + 1)
            lastrun.Tot_Pop_1yr_5_year_olds = squeeze(sum(sum(lastrun.Tot_Pop_1yr(:,i5y,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)
            lastrun.Tot_Pop_1yr_under_5_year_olds = squeeze(sum(sum(lastrun.Tot_Pop_1yr(:,index_under_5y,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)
            lastrun.NumSAg_1yr_5_year_olds = squeeze(sum(sum(lastrun.NumSAg_1yr(:,i5y,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)
            lastrun.NumSAg_1yr_under_5_year_olds = squeeze(sum(sum(lastrun.NumSAg_1yr(:,index_under_5y,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)
            lastrun.NumSAg_chronic_1yr_5_year_olds = squeeze(sum(sum(lastrun.NumSAg_chronic_1yr(:,i5y,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)            
            lastrun.NumSAg_chronic_1yr_under_5_year_olds = squeeze(sum(sum(lastrun.NumSAg_chronic_1yr(:,index_under_5y,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)            
            lastrun.Tot_Pop_1yr = squeeze(sum(sum(lastrun.Tot_Pop_1yr(:,:,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)
            lastrun.num_births_1yr = lastrun.num_births_1yr(i1980:i2100); % 1 x (num_years_simul + 1)
            lastrun.Incid_chronic_all_1yr_approx = squeeze(sum(sum(lastrun.Incid_chronic_all_1yr_approx(:,:,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)
            lastrun.NumSAg_1yr = squeeze(sum(sum(lastrun.NumSAg_1yr(:,:,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)
            lastrun.NumSAg_chronic_1yr = squeeze(sum(sum(lastrun.NumSAg_chronic_1yr(:,:,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)
            lastrun.Prev_TDF_treat_1yr = squeeze(sum(sum(lastrun.Prev_TDF_treat_1yr(:,:,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)
            lastrun.Prev_treatment_eligible_1yr = squeeze(sum(sum(lastrun.Prev_treatment_eligible_1yr(:,:,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)
            lastrun.Incid_Deaths_1yr_approx = squeeze(sum(sum(lastrun.Incid_Deaths_1yr_approx(:,:,i1980:i2100),1),2))'; % 2 x 100 x (num_years_simul + 1)
            lastrun.DALYPerYear = squeeze(sum(lastrun.DALYPerYear(:,i1980:i2100),1)); % 100 x (num_years_simul + 1)
            assert(isequal(size(lastrun.Time),[1 num_cols_out]))
            assert(isequal(size(lastrun.Tot_Pop_1yr_5_year_olds),[1 num_cols_out]))
            assert(isequal(size(lastrun.Tot_Pop_1yr_under_5_year_olds),[1 num_cols_out]))
            assert(isequal(size(lastrun.Tot_Pop_1yr),[1 num_cols_out]))
            assert(isequal(size(lastrun.DALYPerYear),[1 num_cols_out]))
            assert(all(lastrun.Tot_Pop_1yr>=lastrun.NumSAg_1yr))
            assert(all(lastrun.NumSAg_1yr>=lastrun.NumSAg_chronic_1yr))
            disp([lastrun.NumSAg_chronic_1yr;lastrun.Prev_treatment_eligible_1yr])
            assert(all(lastrun.NumSAg_chronic_1yr>=lastrun.Prev_treatment_eligible_1yr))
            assert(all(lastrun.Prev_treatment_eligible_1yr>=lastrun.Prev_TDF_treat_1yr))
            assert(all(lastrun.Tot_Pop_1yr>=lastrun.Incid_Deaths_1yr_approx))
            assert(all(lastrun.Tot_Pop_1yr>=lastrun.Tot_Pop_1yr_5_year_olds))
            assert(all(lastrun.Tot_Pop_1yr>=lastrun.Tot_Pop_1yr_under_5_year_olds))
            assert(all(lastrun.Tot_Pop_1yr_5_year_olds>=lastrun.NumSAg_1yr_5_year_olds))
            assert(all(lastrun.Tot_Pop_1yr_under_5_year_olds>=lastrun.NumSAg_1yr_under_5_year_olds))
            assert(all(lastrun.NumSAg_1yr_5_year_olds>=lastrun.NumSAg_chronic_1yr_5_year_olds))
            assert(all(lastrun.NumSAg_1yr_under_5_year_olds>=lastrun.NumSAg_chronic_1yr_under_5_year_olds))


            if strcmp(stochas_run_str,'1')
                lastrun.country_name = params.country_name;
                %%lastrun.scenario = scenario;

                i1980 = find(years_vec_01yr>=1980,1);
                i2100 = find(years_vec_01yr>=2100,1);
                lastrun.InfantVacc = scenario_HepB3coverage(i1980:i2100);
                lastrun.BirthDoseVacc = scenario_BDcoverage(i1980:i2100);
            end


            countryMap(ISO) = lastrun;

        end % end for country_num loop

        outMap(num2str(scenario_num)) = countryMap;
        %% MP: I removed this because I never use these files:
        %%save(fullfile(basedir,'outputs',filename_results),'outMap') 
        %%if strcmp(stochas_run_str,'1') && strcmp(sensitivity_analysis,'default')
        %%    save(fullfile(basedir,'outputs','scenarios_array.mat'),'label_array') % only version saved after last scenario is correct
        %%end

        time_taken_for_scenario = datetime('now') - begin_time_scenario;
        scenario_hours_vec(scenario_num) = time_taken_for_scenario;
        %%assert(all(scenario_hours_vec(1:scenario_num)>0))
        %%average_time_per_scenario = mean(scenario_hours_vec(1:scenario_num));
        %% min_time_per_scenario = min(scenario_hours_vec(1:scenario_num)); %% MP: not used.
        %% max_time_per_scenario = max(scenario_hours_vec(1:scenario_num)); %% MP: not used.
        %%num_scenarios_left = num_scenarios - scenario_num;
        %%mean_time_left = num_scenarios_left * average_time_per_scenario;
        %% min_time_left = num_scenarios_left * min_time_per_scenario; %% MP: not used.
        %% max_time_left = num_scenarios_left * max_time_per_scenario; %% MP: not used.
        %%if num_scenarios_left>0
        %%    disp(append('There are ',num2str(num_scenarios_left),' scenarios left for run number ',stochas_run_str,' (',sensitivity_analysis,'), which will take about ',char(mean_time_left),' hh:mm:ss.'));
        %%end
    

    end % end for scenario_num loop
    disp("BBB")
    %%assert(length(scenario_hours_vec)==num_scenarios)
    %%end_time_run_num = datetime('now');
    %%disp(end_time_run_num)
    %%time_taken_for_run = end_time_run_num - begin_time_run_num;
    %%disp(['The duration of run number ' stochas_run_str ' (of ' num2str(num_stochas_runs) ') was ' char(time_taken_for_run) ' hours.\n\n'])
    %%if stochas_run_num<num_stochas_runs
    %%    num_runs_left = num_stochas_runs - stochas_run_num;
    %%    approximate_time_left = num_runs_left * time_taken_for_run;
    %%    disp(['There are ' num2str(num_runs_left) ' runs left for sensitivity analysis "' sensitivity_analysis '", which will take about ' char(approximate_time_left) ' hours.'])
    %%end

end % end function country_level_analyses








function out = make_coverage_vec(start_year,num_year_divisions,dt,end_year,yleft,xright,yright,last_available_year)
% yleft contains Montagu coverage values from 1980 to last_available_year
% xright contains boundary years from (last_available_year+dt) onwards
% yright contains coverage values from (last_available_year+dt) to 2101 for boundary years in xright

% expand yleft
% expand yright
% join them
    assert(yleft(1)==0)
    x_vec_before_1980 = start_year:dt:(1980-dt);
    num_time_steps_before_1980 = length(x_vec_before_1980);
    y_vec_before_1980 = repmat(yleft(1),1,num_time_steps_before_1980);

    x_vec_1980_last_available_year = 1980:dt:(last_available_year-dt);
    num_time_steps_1980_last_available_year = length(x_vec_1980_last_available_year);
    y_vec_1980_last_available_year = repmat(yleft(1:(end-1)),num_year_divisions,1);
    %%size(y_vec_1980_last_available_year)
    assert(isequal(size(y_vec_1980_last_available_year),[num_year_divisions (last_available_year-1)-1979]))
    y_vec_1980_last_available_year = y_vec_1980_last_available_year(:);
    assert(isequal(size(y_vec_1980_last_available_year),[num_time_steps_1980_last_available_year 1]))

    assert(xright(1)==last_available_year)
    outright = interp1(xright, yright, last_available_year:dt:end_year, 'linear');

    out = [y_vec_before_1980 y_vec_1980_last_available_year' outright];
    assert(isequal(size(out),size(start_year:dt:end_year)))

end


%% Extract % from 'Status quo infant & BD expansion to 25%' etc
% function pc = extract_percent_from_BDexpansion_scenario_label(s)
%     string_array = strsplit(s);
%     percent_with_symbol = string_array{end};
%     pc = str2double(replace(percent_with_symbol,'%',''));
% end
% 
% %% Extract % from 'Status quo infant & BD drop 5 2020' etc
% function pc = extract_percent_from_BDdrop_scenario_label(s)
%     string_array = strsplit(s);
%     pc_cell=string_array(7);
%     pc = str2double(cell2mat(pc_cell));
% end


% function coverage = PAP_coverage_scaleup(start_year_simul, TScaleup_PAP_start, TScaleup_PAP_end,...
%     last_year_run, dt, PAP_coverage_thissubgroup)
%     xvals_vec = [start_year_simul TScaleup_PAP_start TScaleup_PAP_end last_year_run];
%     % Scales up linearly from 0 to PAP_coverage_thissubgroup over the period
%     % (TScaleup_PAP-1) to TScaleup_PAP
%     yvals_vec = [0 0 PAP_coverage_thissubgroup PAP_coverage_thissubgroup];
% 
%     TimeSteps = start_year_simul:dt:last_year_run; % 1 x 2101 double; [1890 1890.1 1890.2 ... 2099.8 2099.9 2100 2100.1 ... 2100.8 2100.9 2101]
% 
% 
%     coverage = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
%     % Ensure coverage is capped at 100%:
%     coverage = min(1,coverage); 
% end
% 

function coverage = PAP_coverage_scaleup(start_year_simul, Past_TScaleup_PAP_start, Past_TScaleup_PAP_end,...
    Intervention_TScaleup_PAP_start, Intervention_TScaleup_PAP_end, last_year_run, ...
    Current_PAP_coverage_thissubgroup, Max_PAP_coverage_thissubgroup, dt)

    xvals_vec = [start_year_simul Past_TScaleup_PAP_start Past_TScaleup_PAP_end Intervention_TScaleup_PAP_start Intervention_TScaleup_PAP_end last_year_run];
    % Scales up linearly from 0 to PAP_coverage_thissubgroup over the period
    % (TScaleup_PAP-1) to TScaleup_PAP
    yvals_vec = [0 0 Current_PAP_coverage_thissubgroup Current_PAP_coverage_thissubgroup Max_PAP_coverage_thissubgroup Max_PAP_coverage_thissubgroup];

    TimeSteps = start_year_simul:dt:last_year_run; % 1 x 2101 double; [1890 1890.1 1890.2 ... 2099.8 2099.9 2100 2100.1 ... 2100.8 2100.9 2101]


    coverage = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
    % Ensure coverage is capped at 100%:
    coverage = min(1,coverage); 
end


%% Dx and TxifDx time-trends. Treatment (and diagnosis) is assumed to begin in t0_treatment at coverage level "coverage_t0" 
%% (t0_treatment is set as 2016 in the main code).
function coverage = Dx_and_Tx_coverage_scaleup(start_year_simul, t0_treatment, treatment_intervention_start,...
    last_year_run, coverage_t0, historic_annual_increase_coverage, ...
    intervention_annual_increase_coverage, ceiling_coverage, dt)

    xvals_vec = [start_year_simul (t0_treatment-dt) t0_treatment treatment_intervention_start last_year_run];

    % Scales up linearly from 0 to PAP_coverage_thissubgroup over the period
    % (TScaleup_PAP-1) to TScaleup_PAP
    coverage_now = coverage_t0 + (treatment_intervention_start-t0_treatment)*historic_annual_increase_coverage;
    %% We allow this coverage to be >1 (i.e. above 100%) - we cap the coverage later on.
    coverage_max = coverage_now + (last_year_run-treatment_intervention_start)*intervention_annual_increase_coverage;
    yvals_vec = [0 0 coverage_t0 coverage_now coverage_max];

    TimeSteps = start_year_simul:dt:last_year_run; % 1 x 2101 double; [1890 1890.1 1890.2 ... 2099.8 2099.9 2100 2100.1 ... 2100.8 2100.9 2101]
    coverage = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
    
    assert(ceiling_coverage<=1)
    %% Here we ensure that coverage saturates (at a value <100%):
    coverage = min(ceiling_coverage,coverage); 
    
end



% %% Function is used as part of the code to set up the Transitions matrices (disease progression).
% %% Eligibility is determined separately (in get_treatment_eligible_ageindices()).
% %% This is a code snippet that determines whether the current transition
% function RRprogression = get_RRprogression_type(start_state, end_state, i_natural_hist, ...
%     RRprogress_effective_treatment_CC, RRprogress_nonadherent_treatment_CC, ...
%     RRprogress_effective_treatment_nonCC, RRprogress_nonadherent_treatment_nonCC)
%     %% Different relative rate multiplier when on treatment (vs not) for CC/death versus other disease progressions:
%     if(end_state==i_natural_hist.HCC || end_state==i_natural_hist.HBVdeath)
%         thisRR_effective_treatment = RRprogress_effective_treatment_CC;
%         thisRR_nonadherent_treatment = RRprogress_nonadherent_treatment_CC;
%     elseif(start_state==i_natural_hist.ImmReact || start_state==i_natural_hist.Chronic ...
%             || start_state==i_natural_hist.CompCirr || start_state==i_natural_hist.DecompCirr)
%         %% These are the "non-CC" transitions:
%         thisRR_effective_treatment = RRprogress_effective_treatment_nonCC;
%         thisRR_nonadherent_treatment = RRprogress_nonadherent_treatment_nonCC;
%     else
%         %% No effect (because not on treatment):
%         thisRR_effective_treatment = 1;
%         thisRR_nonadherent_treatment = 1;
%     end 
% 
%     RRprogression = [thisRR_effective_treatment, thisRR_nonadherent_treatment];
% end
