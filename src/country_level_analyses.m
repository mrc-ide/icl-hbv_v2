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
    basedir,...
    num_states,num_year_divisions,dt,ages,num_age_steps,start_year,num_years_simul,end_year,...
    theta,CFR_Acute,rate_6months,ECofactor,p_ChronicCarriage,life_expectancy, Prog)

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
    
    % TUTAJ:
    num_scenarios = 11;
    %start_scenario = 17;
    start_scenario = 6;

    assert(ismember(sensitivity_analysis,{'default','infant_100','treat_medium','treat_high'}))

    stochas_run_num = str2double(stochas_run_str);

    filename_results = ['results_countries_', sensitivity_analysis, '_stochastic_run_', stochas_run_str, '.mat'];


    %%begin_time_run_num = datetime('now');
    %disp(append('Run number ',stochas_run_str,' (of ',num2str(num_stochas_runs),') started at ',string(begin_time_run_num)));

    outMap = containers.Map; 
    % for this particle (stochas_run_num), outMap contains the num_scenario scenarios, each of which contains countryMap, which contains the model results (lastrun) for 110 countries
    label_array = cell(1,num_scenarios);
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
                %% treatment_rate_params.Treatmentrate_2016 = stochas_params_mat(stochas_run_num,country_start_col+7))
                HBsAg_treat_cov_all_ages = pop_size_HBsAg_treatment_2016_vec(4);
                assert(HBsAg_treat_cov_all_ages>0 && HBsAg_treat_cov_all_ages<1, "HBsAg_treat_cov_all_ages must be between 0 and 1")
                assert(in_treatment_2016_CDA==pop_size_HBsAg_treatment_2016_vec(5))
            else
                HBsAg_treat_cov_all_ages = 0;
            end

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
           
            

            Prog_scenario(8, 11) = params.CancerDeathRate;  % HCC to HBV death.


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


            %% TAM: assign PAP/VL parameters:
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
            [transactions_from, transactions_to] = find(Prog_scenario > 0);

            % Load these into a data-structure:
            Transactions = [];
            for tr = 1:length(transactions_from)
                Transactions.From(tr) = transactions_from(tr);
                Transactions.To(tr) = transactions_to(tr);
                Transactions.Values(tr) = {repmat(Prog_scenario(transactions_from(tr), transactions_to(tr)), [1, num_age_steps, 2, 2])}; 
                % Transaction.Values is identical across age, gender, treatment
                % For each to-from pair, form a 1 x num_age_steps x 2 x 2 double containing the progression parameter for that to-from transition
                % Arrange this sequence of matrices in a cell array called Transactions.Values, which is contained in Transactions
            end


            % Add age-specific progression acute severe/non-severe to immune tolerant or immune (14, 15)-->(2, 9)
            Transactions.From(length(Transactions.From) + 1) = 14;
            Transactions.To(length(Transactions.To) + 1) = 2;
            Transactions.Values(length(Transactions.Values) + 1) = {p_ChronicCarriage * rate_6months};

            Transactions.From(length(Transactions.From) + 1) = 14;
            Transactions.To(length(Transactions.To) + 1) = 9;
            Transactions.Values(length(Transactions.Values) + 1) = {(1 - p_ChronicCarriage) * rate_6months};

            Transactions.From(length(Transactions.From) + 1) = 15;
            Transactions.To(length(Transactions.To) + 1) = 2;
            Transactions.Values(length(Transactions.Values) + 1) = {p_ChronicCarriage * (1 - CFR_Acute) * rate_6months};

            Transactions.From(length(Transactions.From) + 1) = 15;
            Transactions.To(length(Transactions.To) + 1) = 9;
            Transactions.Values(length(Transactions.Values) + 1) = {(1 - p_ChronicCarriage) * (1 - CFR_Acute) * rate_6months};



            % Add age-specific progression immune tolerant -> immune reactive -> asymptomatic (2,3) and (3,4)
            AgeSpecELossFunction=params.SpeedUpELoss_F*exp(-params.SpeedUpELoss_Beta*ages);

            Transactions.From(length(Transactions.From)+1)=2;
            Transactions.To(length(Transactions.To)+1)=3;
            Transactions.Values(length(Transactions.Values)+1)={repmat(Prog_scenario(2,3)*AgeSpecELossFunction,[1 1 2 2])};

            Transactions.From(length(Transactions.From)+1)=3;
            Transactions.To(length(Transactions.To)+1)=4;
            Transactions.Values(length(Transactions.Values)+1)={repmat(Prog_scenario(3,4)*AgeSpecELossFunction,[1 1 2 2])};



            % Modify immune reactive to chronic (3,5)  to an age-specific progression
            shortcut_rate_age = 20;
            indicator_vec = [zeros(1,shortcut_rate_age*10) ones(1,num_age_steps - shortcut_rate_age*10)];
            tmp_pos = find((Transactions.From == 3) & (Transactions.To == 5));
            assert(isscalar(tmp_pos))
            tmp_mat = Transactions.Values(tmp_pos);
            tmp_mat = tmp_mat{1};
            assert(min(min(min(min(tmp_mat))))==max(max(max(max(tmp_mat)))))
            Transactions.Values(tmp_pos) = {repmat(Prog_scenario(3,5)*indicator_vec,[1 1 2 2])}; % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair

            trans_rate_age = 25; % generic setting - at this age the rate of transition from chronic Hep B to comp cirrhosis is minimised (it's 0).
            trans_rate_by_age = (params.cirrh_rate_coeff*(ages - trans_rate_age)).^2; 
            trans_rate_by_age = trans_rate_by_age .* [zeros(1,trans_rate_age*10) ones(1,num_age_steps - trans_rate_age*10)];

            % Modify Chronic Hep B to Comp Cirrhosis (5,6) to an
            % age-specific progression:
            trans_rate_by_age_5_6 = Prog_scenario(5,6)*trans_rate_by_age;
            tmp_pos = find((Transactions.From == 5) & (Transactions.To == 6));
            assert(isscalar(tmp_pos))
            tmp_mat = Transactions.Values(tmp_pos);
            tmp_mat = tmp_mat{1};
            assert(min(min(min(min(tmp_mat))))==max(max(max(max(tmp_mat)))))
            tmp_mat = repmat(trans_rate_by_age_5_6,[1 1 2 2]); % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
            tmp_mat(:,:,1,:) = tmp_mat(:,:,1,:)*params.CirrhosisRate_WomenCoFactor;
            tmp_mat(:,:,2,:) = tmp_mat(:,:,2,:)*params.CirrhosisRate_MenCoFactor;
            tmp_mat = min(5, tmp_mat);
            Transactions.Values(tmp_pos) = {tmp_mat};

            % Modify Immune Reactive to Comp Cirrhosis (3,6) to an age-specific progression
            trans_rate_by_age_3_6 = Prog_scenario(3,6)*trans_rate_by_age;
            tmp_pos = find((Transactions.From == 3) & (Transactions.To == 6));
            assert(isscalar(tmp_pos))
            tmp_mat = Transactions.Values(tmp_pos);
            tmp_mat = tmp_mat{1};
            assert(min(min(min(min(tmp_mat))))==max(max(max(max(tmp_mat)))))
            tmp_mat = repmat(trans_rate_by_age_3_6,[1 1 2 2]); % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
            tmp_mat(:,:,1,:) = tmp_mat(:,:,1,:)*params.CirrhosisRate_WomenCoFactor;
            tmp_mat(:,:,2,:) = tmp_mat(:,:,2,:)*params.CirrhosisRate_MenCoFactor;
            tmp_mat = min(5, tmp_mat);
            Transactions.Values(tmp_pos) = {tmp_mat};


            % Add sex-specific co-factor to clearance (asympt -> immune) (4, 9)
            tmp_pos = find((Transactions.From == 4) & (Transactions.To == 9));
            assert(isscalar(tmp_pos))
            tmp = Transactions.Values{tmp_pos}; % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
            tmp(:, :, 1, :) = tmp(:, :, 1, :) * params.ClearanceRateWomenCoFactor;
            Transactions.Values(tmp_pos) = {tmp};


            % Add age-specific progresion to HCC (2, 3, 4, 5, 6)-->8
            cancer_rate_age = 10;
            cancer_rate_by_age = (params.cancer_rate_coeff*(ages - cancer_rate_age)).^2; 
            cancer_rate_by_age = cancer_rate_by_age .* [zeros(1,cancer_rate_age*10) ones(1,num_age_steps - cancer_rate_age*10)];
            AgeSpecCancerFunction_Women = min(1, params.CancerRate_WomenCoFactor * cancer_rate_by_age);
            AgeSpecCancerFunction_Men = min(1, params.CancerRate_MenCoFactor * cancer_rate_by_age);

            AgeSpecificProgToCancer = zeros(1, num_age_steps, 2, 2);
            AgeSpecificProgToCancer(:, :, 1, :) = repmat(AgeSpecCancerFunction_Women, [1, 1, 1, 2]);
            AgeSpecificProgToCancer(:, :, 2, :) = repmat(AgeSpecCancerFunction_Men, [1, 1, 1, 2]);

            Transactions.From(length(Transactions.From) + 1) = 2;
            Transactions.To(length(Transactions.To) + 1) = 8;
            Transactions.Values(length(Transactions.Values) + 1) = {AgeSpecificProgToCancer};

            Transactions.From(length(Transactions.From) + 1) = 3;
            Transactions.To(length(Transactions.To) + 1) = 8;
            Transactions.Values(length(Transactions.Values) + 1) = {2 * AgeSpecificProgToCancer};

            Transactions.From(length(Transactions.From) + 1) = 4;
            Transactions.To(length(Transactions.To) + 1) = 8;
            Transactions.Values(length(Transactions.Values) + 1) = {0.5 * AgeSpecificProgToCancer};

            Transactions.From(length(Transactions.From) + 1) = 5;
            Transactions.To(length(Transactions.To) + 1) = 8;
            Transactions.Values(length(Transactions.Values) + 1) = {2 * AgeSpecificProgToCancer};

            Transactions.From(length(Transactions.From) + 1) = 6;
            Transactions.To(length(Transactions.To) + 1) = 8;
            Transactions.Values(length(Transactions.Values) + 1) = {13 * AgeSpecificProgToCancer};




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
            %% Here we define the indices for each scenario we are looking at:
            %%i_scenario_BASE2020notreat = 1;
            %%i_scenario_BASE2025notreat = 1;   %% WUENIC 2025 BD+HepB3, no treatment, no new interventions. Introduced because existing treatment scenarios have treatment coverage increasing to ~70% by 2100.
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

            % i_scenario_ContImp_plusdiag = 6;  %% As i_scenario_ContImp, but maximising PAP+BD coverage (so no overlap if possible)
            % i_scenario_ContImp_plustreat = 7;  %% As i_scenario_ContImp, but maximising PAP+BD coverage (so no overlap if possible)
            % i_scenario_ContImp_plusB3_BD = 8;  
            % i_scenario_ContImp_plusB3_BD_PAP = 9;  %% ContImp + BD to 60% in non-facility births (
            % i_scenario_ContImp_plusB3_BD_PAP_diag = 10;  %% ContImp + 90% PAP coverage of high VL
            % i_scenario_ContImp_plusB3_BD_PAP_diag_treat = 11;  %% As i_scenario_ContImp, but maximising PAP+BD coverage (so no overlap if possible)
            

            % i_scenario_INFACILITYBD_VLPAPtoANC = 4;  %% As i_scenario_INFACILITYBD, but with HVL PAP to ANC coverage levels.
            % i_scenario_INFACILITYBD_UniPAPtoBD = 5;  %% As i_scenario_INFACILITYBD, but with universal PAP to BD coverage levels.
            % i_scenario_INFACILITYBD_UniPAPtoANC = 6;  %% As i_scenario_INFACILITYBD, but with universal PAP to ANC coverage levels.
            % i_scenario_BD75percent = 7;       %% As i_scenario_BASE2025notreat, but BD reaches 75% coverage by T_INTERVENTION_END (needs introduction of new tech like MAP/CPAD - but for now we just take an overall effective coverage).
            % i_scenario_BDWHOtarget = 8;       %% As i_scenario_BASE2025notreat, but BD reaches 90% coverage by T_INTERVENTION_END (needs introduction of new tech like MAP/CPAD - but for now we just take an overall effective coverage).
            % i_scenario_BDWHOtarget_VLPAPtoBD = 9;   %% As i_scenario_BDWHOtarget, but with HVL PAP to BD coverage levels.
            % i_scenario_BDWHOtarget_VLPAPtoANC = 10;   %% As i_scenario_BDWHOtarget, but with HVL PAP to ANC coverage levels.
            % i_scenario_BDWHOtarget_UniPAPtoBD = 11;   %% As i_scenario_BDWHOtarget, but with universal PAP to BD coverage levels.
            % i_scenario_BDWHOtarget_UniPAPtoANC = 12;   %% As i_scenario_BDWHOtarget, but with universal PAP to ANC coverage levels.
            % i_scenario_HepB3WHOtarget = 13;    %% As i_scenario_BDWHOtarget_VLPAPtoANC but also HepB3 reach 90% coverage by T_INTERVENTION_END
            % i_scenario_Treatlink45 = 14;    %% As i_scenario_HepB3WHOtarget, but treatment linkage reach 45% coverage (no increase in diagnosis though). Treatment rate scales up over period T_INTERVENTION_START to T_INTERVENTION_END
            % i_scenario_Treatlink80 = 15;    %% As i_scenario_HepB3WHOtarget, but treatment linkage reach 80% coverage (no increase in diagnosis though). Treatment rate scales up over period T_INTERVENTION_START to T_INTERVENTION_END
            % i_scenario_Diag30 = 16;    %% As i_scenario_Treatlink80, and diagnosis reaches 30%, scaling up over period T_INTERVENTION_START to T_INTERVENTION_END
            % i_scenario_Diag70 = 17;    %% As i_scenario_Treatlink80, and diagnosis reaches 70% (similar to China), scaling up over period T_INTERVENTION_START to T_INTERVENTION_END
            
            
            % %%i_scenario_TreatWHOtarget = 7; %% Treatment reach 80% coverage. Treatment rate scales up over period T_INTERVENTION_START to T_INTERVENTION_END
            % i_scenario_WHOtarget = 8;     %% BD+HepB3 reach 90% coverage by T_INTERVENTION_END, treatment reaches 80% by T_INTERVENTION_END
            % i_scenario_MAP = 9; %% WUENIC 2025 BD+HepB3, 2016 treatment, Microarray patch introduced in T_INTERVENTION_START (increase BD coverage, but lower efficacy).
            % i_scenario_CPAD = 10; %% WUENIC 2025 BD+HepB3, 2016 treatment, CPAD patch introduced (increase BD but lower eff and different cost to MAP).
            % i_scenario_BD2025_birthcohorttest = 11;   %% WUENIC 2025 BD+HepB3, 2016 treatment, Thai-B-type testing of pre-BD birth cohort on top of existing testing (cap so cannot test >100% of any age stratum).
            % i_scenario_PAP_TREAThighVL = 12;   %% PAP for High VL pregnant women
            % i_scenario_PAP_TREATeAgpos = 13;             %% eAg+
            % i_scenario_PAP_TREAT_highVL_or_eAgpos = 14;  %% Either high VL or eAg+ (or both)
            % i_scenario_BASE2020_WITHTREAT = 15; %% The 'default' scenario - WUENIC 2019 BD+HepB3, 2016 treatment, no new interventions.
            % i_scenario_BASE2025_WITHTREAT = 16;     %% WUENIC 2025 BD+HepB3, 2016 treatment, no new interventions. Addresses - how have changes in BD+Hep B3 coverage impacted result?

            %%i_scenario_BD2025_LA_TDF = 6; %% WUENIC 2025 BD+HepB3, 2016 treatment, long-acting treatment introduced (increases coverage of TDF treatment).
            %%i_scenario_BD2025_PoC_ALT_HBcrAg = 7;    %% WUENIC 2025 BD+HepB3, 2016 treatment, PoC ALT and HBcrAg introduced - higher treatment coverage, also some people on treatment who don't need it.
            %%i_scenario_BD2025_cure = 9; %% WUENIC 2025 BD+HepB3, 2016 treatment, (hypothetical) cure replaces treatment at current test rates.
            
        
            %% Index values for scenario_BD: Governs BD coverage time trends, introduction of different BD devices (MAP, CPAD).
            I_BD_WUENIC2025 = 100;  %% Follow WUENIC2025 and after 2024 coverage remains at last (2024) value
            I_BD_contimp = 101;  %% Possible increase beyond WUENIC2025 at a slow rate (currently 0). Introduction of BD to GAVI-approved countries
            I_BD_IFexpansion = 102; %% BD increases to current in-facilility % (or current coverage if higher) in all countries with BD/GAVI-approved.
            I_BD_OOFexpansion = 103; %% I_BD_IFexpansion with additional intervention to reach BD_oof_coverage=60% of OOF births (e.g. CHW, MAP). BD introduced in countries without BD.
            I_BD_IF_OOFexpansion = 104; %% Combining IF expansion + OOF 

            %% Index values for scenario_HepB3: Hep B3 scenarios. 
            I_HEPB3_WUENIC2025 = 201;  %% Fixed at WUENIC 2025
            I_HEPB3_contimp = 202;      %% Possible increase beyond WUENIC2025 at a slow rate (currently 0).
            I_HEPB3_WHOtarget = 203;        %% Increase to WHO target.
            
            
            
            %% Index values for scenario_PAP: peripartum antiviral prophylaxis (PAP) treatment for HBsAg+ mothers (treat all, treat high VL etc).
            I_PAP_SQ = 301;      %% No PAP unless already present.
            I_PAP_HVL_contimp = 302; %% Increase PAP to 20% of HVL unless already present.
            I_PAP_HVL_targeted = 303; %% Increase PAP to X% of HVL unless already present (X% placeholder but will be based on VL testing availability).
            I_PAP_PoC = 304;      %% Eligibility based on PoC test with given sensitivity, specificity and (globally constant) coverage.
            I_PAP_all = 305; %% All eligible (no VL criteria). Uptake to national ANC-1 value.

            
            %% scenario_Treatment: governs how treatment happens:
            % Modified by treat-all, introduction of PoC HBcrAg, PoC ALT tests, 
    
            I_TREAT_SQ = 1;         %% Current treatment (Capped at most recent treatemnt data - currently Polaris 2025).
            I_TREAT_continuedimprovement = 2;
            I_TREAT_ANCscreening = 3;
            I_TREAT_IFscreeing = 4;
            I_TREAT_IntegratedServices = 5; 
            I_TREAT_PoCeligibility = 6;
            I_TREAT_universal = 7;
            I_TREAT_LA = 8;
            I_TREAT_decentralised = 9;
            I_TREAT_cureBepi = 10;
            I_TREAT_curev2 = 11;
            

            
       
            %% For each scenario we determine which intervention "levers" are used.
            switch scenario_num
                case i_scenario_SQ     %% WUENIC 2025 BD+HepB3, 2016 treatment, no new interventions. Addresses - how have changes in BD+Hep B3 coverage impacted result?
                    scenario_BD = I_BD_WUENIC2025;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_SQ;
                    scenario_Treatment = I_TREAT_SQ;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp    %% Continued improvement - BD can increase (+ starts up in GAVI-approved countries). HepB3 can in crease
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusB3     %% Hep B3 increases to 90% 2026-2029 (T_INTERVENTION_START_HepB3-T_INTERVENTION_END_HepB3)
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_WHOtarget;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";                    
                case i_scenario_ContImp_plusBD_IF     %% BD increases - increasing OOF coverage
                    scenario_BD = I_BD_IFexpansion;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusBD_OOF     %% BD increases - increasing OOF coverage
                    scenario_BD = I_BD_OOFexpansion;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusBD_IF_OOF     %% BD increases - increasing IF+OOF coverage
                    scenario_BD = I_BD_IF_OOFexpansion;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusPAP_HVL_targeted     %% ContImp+ PAP for HVL only, capped at level of availability of VL testing.
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_targeted;
                    scenario_Treatment = I_TREAT_continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusPAP_PoC     %% ContImp+ PAP eligibility through PoC test.
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_PoC;
                    scenario_Treatment = I_TREAT_continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";
                case i_scenario_ContImp_plusPAP_all     %% ContImp+ PAP eligibility through PoC test.
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_all;  %% PAP available to all pregnant women regardless of VL
                    scenario_Treatment = I_TREAT_continuedimprovement;
                    scenario_AddScreenIntervention = "No additional screening";  
                WYBOR:
                %% ContImp with additional diagnosis among pregnant women during ANC (ANC-1 capped at HIV screening %)
                case i_scenario_ContImp_plusDx_ANCscreening
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_ANCscreening;
                    scenario_AddScreenIntervention = "ANC screening";  
                 %% ContImp with additional diagnosis through birth cohort screening (those born within 5 years of country introduction of HepB3)
                case i_scenario_ContImp_plusDx_BirthCohort
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_continuedimprovement;
                    scenario_AddScreenIntervention = "Birth cohort screening";  
                case i_scenario_ContImp_plusDx_IFscreening
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_IFscreeing;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% ContImp with additional diagnosis through (realistic) community screening similar to PROLIFICA.
                case i_scenario_ContImp_plusDx_CommunityScreening
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_continuedimprovement;
                    scenario_AddScreenIntervention = "Community screening";  
                %% ContImp with additional diagnosis through perfect community screening.
                case i_scenario_ContImp_plusDx_CommunityScreening_perfect
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_continuedimprovement;
                    scenario_AddScreenIntervention = "Perfect community screening";  
                %% CompInt with integrated services (TBD?)
                case i_scenario_ContImp_plusDx_IntegratedServices
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_IntegratedServices;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% CompInt with PoC treatment eligibility 
                case i_scenario_ContImp_plusTx_PoCeligibility
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_PoCeligibility;
                    scenario_AddScreenIntervention = "No additional screening";  
                case i_scenario_ContImp_plusTx_treatall
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_universal;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% Long-acting treatment available:
                case i_scenario_ContImp_plusTx_LA
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_LA;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% Decentralised testing and treatment:
                case i_scenario_ContImp_plusDecentralisedDxTx
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_decentralised;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% Bepi-like cure available:
                case i_scenario_ContImp_plusTx_cure_Bepi
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_cureBepi;
                    scenario_AddScreenIntervention = "No additional screening";  
                %% Better-than-Bepi cure available:
                case i_scenario_ContImp_plusTx_cure_improved
                    scenario_BD = I_BD_contimp;
                    scenario_HepB3 = I_HEPB3_contimp;
                    scenario_PAP = I_PAP_HVL_contimp;
                    scenario_Treatment = I_TREAT_curev2;
                    scenario_AddScreenIntervention = "No additional screening";  
                % case i_scenario_ContImp_plusdiag    %% Diagnosis increases to 70%     
                % case i_scenario_ContImp_plustreat    %% Treatment (of diagnosed+eligible) increases to 45%     
                                
                % case i_scenario_INFACILITYBD_VLPAPtoANC     %% WUENIC 2025 HepB3, 2016 treatment, BD introduced in countries where it is not already present - coverage capped at in-facility birth coverage.
                %     scenario_PAP = I_PAP_highVL_ANCcoverage;
                % case i_scenario_INFACILITYBD_UniPAPtoANC     %% WUENIC 2025 HepB3, 2016 treatment, BD introduced in countries where it is not already present - coverage capped at in-facility birth coverage.
                %     scenario_PAP = I_PAP_universal_ANCcoverage;
                % case i_scenario_BDWHOtarget    %% BD reaches 90% coverage target
                %     scenario_BD = I_BD_WHOtarget;
                % case i_scenario_BDWHOtarget_VLPAPtoBD    %% BD reaches 90% coverage target, PAP for HVL to BD coverage levels
                %     scenario_BD = I_BD_WHOtarget;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_highVL_BDcoverage;
                %     scenario_Treatment = I_TREAT_SQ;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_BDWHOtarget_VLPAPtoANC    %% BD reaches 90% coverage target, PAP for HVL to ANC coverage levels
                %     scenario_BD = I_BD_WHOtarget;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_highVL_ANCcoverage;
                %     scenario_Treatment = I_TREAT_SQ;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_BDWHOtarget_UniPAPtoBD    %% BD reaches 90% coverage target, universal PAP to BD coverage levels
                %     scenario_BD = I_BD_WHOtarget;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_universal_BDcoverage;
                %     scenario_Treatment = I_TREAT_SQ;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_BDWHOtarget_UniPAPtoANC    %% BD reaches 90% coverage target, universal PAP to ANC coverage levels
                %     scenario_BD = I_BD_WHOtarget;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_universal_ANCcoverage;
                %     scenario_Treatment = I_TREAT_SQ;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_HepB3WHOtarget  %% HepB3 reaches 90% coverage target
                %     scenario_BD = I_BD_WHOtarget;
                %     scenario_HepB3 = I_HEPB3_WHOtarget;
                %     scenario_PAP = I_PAP_highVL_ANCcoverage;
                %     scenario_Treatment = I_TREAT_SQ;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_Treatlink45  %% 45% of those diag+elig start treatment
                %     scenario_BD = I_BD_WHOtarget;
                %     scenario_HepB3 = I_HEPB3_WHOtarget;
                %     scenario_PAP = I_PAP_highVL_ANCcoverage;
                %     scenario_Treatment = I_TREATlink45;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_Treatlink80  %% 80% of those diag+elig start treatment
                %     scenario_BD = I_BD_WHOtarget;
                %     scenario_HepB3 = I_HEPB3_WHOtarget;
                %     scenario_PAP = I_PAP_highVL_ANCcoverage;
                %     scenario_Treatment = I_TREATlink80;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_Diag30  %% 30% diagnosed. Treatment linkage as for I_TREATlink80.
                %     scenario_BD = I_BD_WHOtarget;
                %     scenario_HepB3 = I_HEPB3_WHOtarget;
                %     scenario_PAP = I_PAP_highVL_ANCcoverage;
                %     scenario_Treatment = I_diag30percent;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_Diag70  %% 70% diagnosed. Treatment linkage as for I_TREATlink80.
                %     scenario_BD = I_BD_WHOtarget;
                %     scenario_HepB3 = I_HEPB3_WHOtarget;
                %     scenario_PAP = I_PAP_highVL_ANCcoverage;
                %     scenario_Treatment = I_diag70percent;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_TreatWHOtarget  %% Treatment reaches 80% target
                %     disp("Treatment target scenario")
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_NOTREAT;
                %     scenario_Treatment = I_TREAT_WHOtarget;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_WHOtarget       %% BD, HepB3 and treatment reach targets
                %     disp("WHO target scenario")
                %     scenario_BD = I_BD_WHOtarget;
                %     scenario_HepB3 = I_HEPB3_WHOtarget;
                %     scenario_PAP = I_PAP_highVL_ANCcoverage;
                %     scenario_Treatment = I_TREAT_WHOtarget;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_MAP %% WUENIC 2025 BD+HepB3, 2016 treatment, Microarray patch introduced in T_INTERVENTION_START (increase BD coverage, but lower efficacy).
                %     scenario_BD = I_BD_MAP;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_NOTREAT;
                %     scenario_Treatment = I_NOTREAT;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_CPAD %% WUENIC 2025 BD+HepB3, 2016 treatment, CPAD patch introduced (increase BD but lower eff and different cost to MAP).
                %     scenario_BD = I_BD_CPAD;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_NOTREAT;
                % case i_scenario_BD2025_birthcohorttest   %% WUENIC 2025 BD+HepB3, 2016 treatment, Thai-B-type testing of pre-BD birth cohort on top of existing testing (cap so cannot test >100% of any age stratum).
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_NOTREAT;
                %     scenario_Treatment = I_NOTREAT;
                %     scenario_AddScreenIntervention = I_BIRTHCOHORT_SCREENING;     
                % case i_scenario_PAP_TREAThighVL
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_TREAThighVL;
                %     scenario_Treatment = I_NOTREAT;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_PAP_TREATeAgpos
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_TREATeAgpos;
                %     scenario_Treatment = I_NOTREAT;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_PAP_TREAT_highVL_or_eAgpos
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_TREAT_highVL_or_eAgpos;
                %     scenario_Treatment = I_NOTREAT;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_BASE2020_WITHTREAT   %%case 'Status quo infant & BD'
                %     scenario_BD = I_BD_WUENIC2020;
                %     scenario_HepB3 = I_HEPB3_WUENIC2020;
                %     scenario_PAP = I_PAP_NOTREAT;
                %     scenario_Treatment = I_TREAT_SQ;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_BASE2025_WITHTREAT     %% WUENIC 2025 BD+HepB3, 2016 treatment, no new interventions. Addresses - how have changes in BD+Hep B3 coverage impacted result?
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_NOTREAT;
                %     scenario_Treatment = I_TREAT_SQ;
                %     scenario_AddScreenIntervention = "No additional screening";

                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_BD2025_LA_TDF %% WUENIC 2025 BD+HepB3, 2016 treatment, long-acting treatment introduced (increases coverage of TDF treatment).
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_NOTREAT;
                %     scenario_Treatment = I_TREAT_INIT_LA;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_BD2025_PoC_ALT_HBcrAg   %% WUENIC 2025 BD+HepB3, 2016 treatment, PoC ALT and HBcrAg introduced - higher treatment coverage, also some people on treatment who don't need it.
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_NOTREAT;
                %     scenario_Treatment = I_TREAT_INIT_POC_cr_ALT;
                %     scenario_AddScreenIntervention = "No additional screening";
                % case i_scenario_BD2025_cure
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_NOTREAT;
                %     scenario_Treatment = I_CURE;
                %     scenario_AddScreenIntervention = "No additional screening";
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
                case I_BD_IF_OOFexpansion
                    BIDOOF
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
            disp([start_year,num_year_divisions,dt,end_year])
            disp(coverage_BD_to_last_datapoint)
            disp("A")
            disp(future_xvals_vec)
            disp("B")
            disp(future_yvals_vec)
            disp(year_last_BD_data)
            disp("DONE")
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
                    if(potential_hepb3_target>HepB3_WHO_target_coverage)
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
                    PAP_current_coverage = 0.26
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
                    PoC_coverage = Intervention_data_thiscountry.PAP_PoC_coverage;
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
        

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Diagnosis and treatment:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            treatment_rate_params = struct('Treatmentrate_2016', stochas_params_mat(stochas_run_num,country_start_col+7),...
                    'Treatmentrate_final', 0,...
                    't_treatment_scaleup_start', T_INTERVENTION_START,...
                    't_treatment_scaleup_end', T_INTERVENTION_END,...
                    't_remove_treatment_barriers', T_INTERVENTION_START,...
                    'prop_diagnosed_t0',0,...
                    'prop_treatifdiag_t0',0,...
                    'annual_increase_diagnosis',0,...
                    'annual_increase_treatifdiag',0);
                    %%'prop_wouldseektreat_t0',0,...
                    %%'prop_wouldseektreat_tchange',0.5);
            switch scenario_Treatment
                case I_TREAT_SQ           %% Use current rates of treatment uptake and failure.
                    %%params.PriorTDFTreatRate = stochas_params_mat(stochas_run_num,country_start_col+7);
                    treatment_rate_params.Treatmentrate_final = treatment_rate_params.Treatmentrate_2016;
                    treatment_rate_params.prop_diagnosed_t0 = Polaris_diagnosis_coverage_map(ISO);
                    treatment_rate_params.prop_treatifdiag_t0 = Polaris_treat_coverage_map(ISO);
                    treatment_rate_params.annual_increase_diagnosis = 0;
                    treatment_rate_params.annual_increase_treatifdiag = 0;
                    max_treatment_coverage = treatment_rate_params.prop_diagnosed_t0*treatment_rate_params.prop_treatifdiag_t0;
                    %%treatment_rate_params.prop_wouldseektreat_tchange = Polaris_diagnosis_coverage_map(ISO)*Polaris_treat_coverage_map(ISO);
                    treatment_rate_params.t_remove_treatment_barriers = 9999; % Dummy value at time beyond any simulation.
                    HAS_TREATMENT = 1; % Treatment introduced in 2016 at rate from MJdV.
                % case I_NOTREAT
                %     treatment_rate_params.Treatmentrate_final = 0; % no treatment
                %     treatment_rate_params.prop_wouldseektreat_t0 = 0;
                %     treatment_rate_params.prop_wouldseektreat_tchange = 0;
                %     treatment_rate_params.t_remove_treatment_barriers = 9999; % Dummy value at time beyond any simulation.
                %     HAS_TREATMENT = 0;  % Treatment never introduced
                case I_TREAT_continuedimprovement
                    treatment_rate_params.annual_increase_diagnosis = Intervention_data_thiscountry.ContImp_Dx_annual_increase;
                    treatment_rate_params.annual_increase_treatifdiag = Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                    treatment_rate_params.Treatmentrate_final = treatment_boundaries_vec(5); % 80% scenario from MJDV code - corresponds to about 0.15/yr for ETH/GMB.
                    treatment_rate_params.prop_diagnosed_t0 = Polaris_diagnosis_coverage_map(ISO);
                    treatment_rate_params.prop_treatifdiag_t0 = Polaris_treat_coverage_map(ISO);
                    max_treatment_coverage = 0.45*0.7;
                    %new_treatlink_prop = max(0.80,Polaris_treat_coverage_map(ISO));
                    %new_diag_prop = max(0.70,Polaris_diagnosis_coverage_map(ISO));
                    %treatment_rate_params.prop_wouldseektreat_tchange = new_diag_prop*new_treatlink_prop;
                    treatment_rate_params.t_remove_treatment_barriers = T_INTERVENTION_START;
                    HAS_TREATMENT = 1;
                case I_TREAT_ANCscreening  %% as Cont Imp treatment except for pregnant women.
                    %% Ceiling = number ANC-1 by age x acceptability of test in a country (use HIV test?)
                    treatment_rate_params.annual_increase_diagnosis = Intervention_data_thiscountry.ContImp_Dx_annual_increase;
                    treatment_rate_params.annual_increase_treatifdiag = Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                    treatment_rate_params.Treatmentrate_final = treatment_boundaries_vec(5); % 80% scenario from MJDV code - corresponds to about 0.15/yr for ETH/GMB.
                    treatment_rate_params.prop_diagnosed_t0 = Polaris_diagnosis_coverage_map(ISO);
                    treatment_rate_params.prop_treatifdiag_t0 = Polaris_treat_coverage_map(ISO);
                    max_treatment_coverage = 0.45*0.7;
                    %new_treatlink_prop = max(0.80,Polaris_treat_coverage_map(ISO));
                    %new_diag_prop = max(0.70,Polaris_diagnosis_coverage_map(ISO));
                    %treatment_rate_params.prop_wouldseektreat_tchange = new_diag_prop*new_treatlink_prop;
                    treatment_rate_params.t_remove_treatment_barriers = T_INTERVENTION_START;
                    HAS_TREATMENT = 1;
                case I_TREAT_IFscreeing
                    treatment_rate_params.annual_increase_diagnosis = Intervention_data_thiscountry.ContImp_Dx_annual_increase;
                    treatment_rate_params.annual_increase_treatifdiag = Intervention_data_thiscountry.ContImp_TxifDx_annual_increase;
                    treatment_rate_params.Treatmentrate_final = treatment_boundaries_vec(5); % 80% scenario from MJDV code - corresponds to about 0.15/yr for ETH/GMB.
                    treatment_rate_params.prop_diagnosed_t0 = Polaris_diagnosis_coverage_map(ISO);
                    treatment_rate_params.prop_treatifdiag_t0 = Polaris_treat_coverage_map(ISO);
                    max_treatment_coverage = 0.45*0.7;
                    %new_treatlink_prop = max(0.80,Polaris_treat_coverage_map(ISO));
                    %new_diag_prop = max(0.70,Polaris_diagnosis_coverage_map(ISO));
                    %treatment_rate_params.prop_wouldseektreat_tchange = new_diag_prop*new_treatlink_prop;
                    treatment_rate_params.t_remove_treatment_barriers = T_INTERVENTION_START;
                    HAS_TREATMENT = 1;
                case I_TREATlink45
                    treatment_rate_params.Treatmentrate_final = treatment_boundaries_vec(5); % 80% scenario from MJDV code - corresponds to about 0.15/yr for ETH/GMB.
                    treatment_rate_params.prop_diagnosed_t0 = Polaris_diagnosis_coverage_map(ISO);
                    treatment_rate_params.prop_treatifdiag_t0 = Polaris_treat_coverage_map(ISO);
                    treatment_rate_params.annual_increase_diagnosis = 0;
                    treatment_rate_params.annual_increase_treatifdiag = max(0, (0.45-Polaris_treat_coverage_map(ISO))/(T_INTERVENTION_END - T_INTERVENTION_START));
                    max_treatment_coverage = treatment_rate_params.prop_diagnosed_t0;
                    %%treatment_rate_params.prop_wouldseektreat_t0 = Polaris_diagnosis_coverage_map(ISO)*Polaris_treat_coverage_map(ISO);
                    %%new_treatlink_prop = max(0.45,Polaris_treat_coverage_map(ISO));
                    %%treatment_rate_params.prop_wouldseektreat_tchange = Polaris_diagnosis_coverage_map(ISO)*new_treatlink_prop;
                    treatment_rate_params.t_remove_treatment_barriers = T_INTERVENTION_START;
                    HAS_TREATMENT = 1;
                case I_diag70percent
                    treatment_rate_params.Treatmentrate_final = treatment_boundaries_vec(5); % 80% scenario from MJDV code - corresponds to about 0.15/yr for ETH/GMB.
                    treatment_rate_params.prop_diagnosed_t0 = Polaris_diagnosis_coverage_map(ISO);
                    treatment_rate_params.prop_treatifdiag_t0 = Polaris_treat_coverage_map(ISO);
                    treatment_rate_params.annual_increase_diagnosis = max(0, (0.70-Polaris_diagnosis_coverage_map(ISO))/(T_INTERVENTION_END - T_INTERVENTION_START));
                    treatment_rate_params.annual_increase_treatifdiag = 0;
                    max_treatment_coverage = treatment_rate_params.prop_treatifdiag_t0;
                    %%treatment_rate_params.prop_wouldseektreat_t0 = Polaris_diagnosis_coverage_map(ISO)*Polaris_treat_coverage_map(ISO);
                    %%new_treatlink_prop = max(0.80,Polaris_treat_coverage_map(ISO));
                    %%treatment_rate_params.prop_wouldseektreat_tchange = Polaris_diagnosis_coverage_map(ISO)*new_treatlink_prop;
                    treatment_rate_params.t_remove_treatment_barriers = T_INTERVENTION_START;
                    HAS_TREATMENT = 1;
                case I_TREAT_PLUS
                    treatment_rate_params.Treatmentrate_final = treatment_boundaries_vec(5); % 80% scenario from MJDV code - corresponds to about 0.15/yr for ETH/GMB.
                    treatment_rate_params.prop_diagnosed_t0 = Polaris_diagnosis_coverage_map(ISO);
                    treatment_rate_params.prop_treatifdiag_t0 = Polaris_treat_coverage_map(ISO);
                    treatment_rate_params.annual_increase_diagnosis = max(0, (0.70-Polaris_diagnosis_coverage_map(ISO))/(T_INTERVENTION_END - T_INTERVENTION_START));
                    treatment_rate_params.annual_increase_treatifdiag = max(0, (0.45-Polaris_treat_coverage_map(ISO))/(T_INTERVENTION_END - T_INTERVENTION_START));
                    max_treatment_coverage = max(0.45,Polaris_treat_coverage_map(ISO)) * max(0.7,Polaris_diagnosis_coverage_map(ISO));
                    %%treatment_rate_params.prop_wouldseektreat_t0 = Polaris_diagnosis_coverage_map(ISO)*Polaris_treat_coverage_map(ISO);
                    %%new_treatlink_prop = max(0.45,Polaris_treat_coverage_map(ISO));
                    %%treatment_rate_params.prop_wouldseektreat_tchange = Polaris_diagnosis_coverage_map(ISO)*new_treatlink_prop;
                    treatment_rate_params.t_remove_treatment_barriers = T_INTERVENTION_START;
                    HAS_TREATMENT = 1;
                % case I_diag70percent
                %     treatment_rate_params.Treatmentrate_final = treatment_boundaries_vec(5); % 80% scenario from MJDV code - corresponds to about 0.15/yr for ETH/GMB.
                %     treatment_rate_params.prop_diagnosed_t0 = Polaris_diagnosis_coverage_map(ISO);
                %     treatment_rate_params.prop_treatifdiag_t0 = Polaris_treat_coverage_map(ISO);
                %     %%treatment_rate_params.prop_wouldseektreat_t0 = Polaris_diagnosis_coverage_map(ISO)*Polaris_treat_coverage_map(ISO);
                %     %%new_treatlink_prop = max(0.80,Polaris_treat_coverage_map(ISO));
                %     %%new_diag_prop = max(0.7,Polaris_diagnosis_coverage_map(ISO));
                %     %%treatment_rate_params.prop_wouldseektreat_tchange = new_diag_prop*new_treatlink_prop;
                %     treatment_rate_params.t_remove_treatment_barriers = T_INTERVENTION_START;
                %     HAS_TREATMENT = 1;

                % case I_TREAT_SQ           %% Use current rates of treatment uptake and failure.
                %     %%params.PriorTDFTreatRate = stochas_params_mat(stochas_run_num,country_start_col+7);
                %     treatment_rate_params.Treatmentrate_final = treatment_rate_params.Treatmentrate_2016;
                %     HAS_TREATMENT = 1; % Treatment introduced in 2016 at rate from MJdV.
                % case I_TREAT_WHOtarget
                %     treatment_rate_params.Treatmentrate_final = treatment_boundaries_vec(5); % 80% scenario from MJDV code
                %     HAS_TREATMENT = 1;
                % case I_TREAT_40percent
                %     treatment_rate_params.Treatmentrate_final = treatment_boundaries_vec(3); % 40% scenario from MJDV code
                %     HAS_TREATMENT = 1;
                
                % case I_TREAT_INIT_POC_cr_ALT   %% Introduce PoC tests for HBcrAg and ALT - increase rate of treatment initiation in eligible groups.
                %     treatment_rate_params.Treatmentrate_final = treatment_rate_params.Treatmentrate_2016;
                % case I_TREAT_INIT_LA           %% Introduce long-acting treatment. SQ treatment failure rate is very low (0.001), so we take the "TDF treatment" group to be "On treatment, adherent and not going to drop out". LA treatment then just increases the proportion of people in this compartment (either by improving adherence, preventing dropout, or offering a more convenient/preffered option).
                %     treatment_rate_params.Treatmentrate_final = treatment_rate_params.Treatmentrate_2016;
                %%case I_CURE                    %% Cure replaces treatment - in this case 
                  
                %     case 'treat_medium'
                %         params.PriorTDFTreatRate = treatment_boundaries_vec(3); % 40%
                %     case 'treat_high'
                %         params.PriorTDFTreatRate = treatment_boundaries_vec(5); % 80%
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


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Now get time trends of PAP coverage (divided into those with/without BD, and by whether EAg+/SAg+ and high/low VL): 
            % Coverage of PAP among those with BD
            %% TAM cov_BirthDoseAndTDF_EAgHighVL_itt
            start_year_simul = 1890;
            last_year_run = end_year;
            PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL = PAP_coverage_scaleup(start_year_simul, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, last_year_run,...
                PAP_cov_params.current_cov_BDandPAP_EAgHighVL, PAP_cov_params.max_cov_BDandPAP_EAgHighVL, dt);

            PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL = PAP_coverage_scaleup(start_year_simul, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, last_year_run,...
                PAP_cov_params.current_cov_BDandPAP_EAgLowVL, PAP_cov_params.max_cov_BDandPAP_EAgLowVL, dt);
            

            PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL = PAP_coverage_scaleup(start_year_simul, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, last_year_run,...
                PAP_cov_params.current_cov_BDandPAP_SAgHighVL, PAP_cov_params.max_cov_BDandPAP_SAgHighVL, dt);
    
            PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL = PAP_coverage_scaleup(start_year_simul, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, last_year_run,...
                PAP_cov_params.current_cov_BDandPAP_SAgLowVL, PAP_cov_params.max_cov_BDandPAP_SAgLowVL, dt);

            % Coverage of PAP among those not with BD
            PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgHighVL = PAP_coverage_scaleup(start_year_simul, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, last_year_run,...
                PAP_cov_params.current_cov_PAPonly_EAgHighVL, PAP_cov_params.max_cov_PAPonly_EAgHighVL, dt);

            PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgLowVL = PAP_coverage_scaleup(start_year_simul, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, last_year_run,...
                PAP_cov_params.current_cov_PAPonly_EAgLowVL, PAP_cov_params.max_cov_PAPonly_EAgLowVL, dt);
            

            PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgHighVL = PAP_coverage_scaleup(start_year_simul, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, last_year_run,...
                PAP_cov_params.current_cov_PAPonly_SAgHighVL, PAP_cov_params.max_cov_PAPonly_SAgHighVL, dt);
    
            PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgLowVL = PAP_coverage_scaleup(start_year_simul, PAP_cov_params.Past_TScaleup_PAP_start, PAP_cov_params.Past_TScaleup_PAP_end,...
                PAP_cov_params.Intervention_TScaleup_PAP_start, PAP_cov_params.Intervention_TScaleup_PAP_end, last_year_run,...
                PAP_cov_params.current_cov_PAPonly_SAgLowVL, PAP_cov_params.max_cov_PAPonly_SAgLowVL, dt);


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
                num_states,num_year_divisions,dt,ages,num_age_steps,start_year,num_years_simul,...
                theta,ECofactor,treatment_rate_params, treatment_start_year-dt, ...
                HBsAg_treat_cov_all_ages, ...
                params, PAP_VL_params, PAP_cov_params, ...
                Global_intervention_params, ...
                p_ChronicCarriage,Prog_scenario,Transactions,......
                scenario_BDcoverage, scenario_BDcoverage_fromMAP,...
                scenario_BDcoverage_fromCPAD, scenario_HepB3coverage, ...
                HAS_TREATMENT, max_treatment_coverage, ...
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
                'Prev_TDF_treat_1yr','Prev_Immune_Reactive_1yr','Prev_Chronic_Hep_B_1yr','Prev_Comp_Cirr_1yr','Prev_Decomp_Cirr_1yr',...
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
                'p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL'...
                };
            assert(all(ismember(fields_of_interest,lastrun_fields)))
            assert(all(ismember(lastrun_fields,fields_of_interest)))

            lastrun.Prev_treatment_eligible_1yr = ...
                lastrun.Prev_Immune_Reactive_1yr + lastrun.Prev_Chronic_Hep_B_1yr + lastrun.Prev_Comp_Cirr_1yr + lastrun.Prev_Decomp_Cirr_1yr + ...
                lastrun.Prev_TDF_treat_1yr;
            lastrun = rmfield(lastrun,'Prev_Immune_Reactive_1yr');
            lastrun = rmfield(lastrun,'Prev_Chronic_Hep_B_1yr');
            lastrun = rmfield(lastrun,'Prev_Comp_Cirr_1yr');
            lastrun = rmfield(lastrun,'Prev_Decomp_Cirr_1yr');

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