function country_level_analyses(sensitivity_analysis,...
    stochas_run_str,...
    ListOfISOs,countries_to_run,...
    BD_table,HepB3_table,...
    num_in_treatment_2016_map,pop_size_HBsAg_treatment_map,treatment_rates_map,...
    country_s_e_HCCdeaths_map,...
    params_map,stochas_params_mat,country_start_cols,...
    WUENIC2024BDdata, WUENIC2024HepB3data, ...
    GHO_infacilitybirthproportion_map, ...
    basedir,...
    num_states,num_year_divisions,dt,ages,num_age_steps,start_year,num_years_simul,end_year,...
    theta,CFR_Acute,rate_6months,ECofactor,p_ChronicCarriage,life_expectancy, Prog)

    if nargin < 1
       error('No input')
    end
    
    T_INTERVENTION_START = 2026.0;
    T_INTERVENTION_END = 2029.0;

    % TUTAJ:
    num_scenarios = 11;
    %start_scenario = 14;
    start_scenario = 11;

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
    %%for scenario_num = 10:10
        disp("Rnning scenario")
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

       
            



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Set up scenarios:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% TUTAJ:
            %% Here we define the indices for each scenario we are looking at:
            iscenario_BASE2020notreat = 1;
            iscenario_BASE2025notreat = 2;   %% WUENIC 2025 BD+HepB3, no treatment, no new interventions. Introduced because existing treatment scenarios have treatment coverage increasing to ~70% by 2100.
            iscenario_INFACILITYBD = 3;  %% WUENIC 2025 HepB3, 2016 treatment, BD introduced in countries where it is not already present - coverage capped at in-facility birth coverage.
            iscenario_BDWHOtarget = 4;   %% BD reaches 90% coverage by T_INTERVENTION_END
            iscenario_HepB3WHOtarget = 5; %% HepB3 reach 90% coverage by T_INTERVENTION_END
            iscenario_Treat40percent = 6; %% Treatment reach 40% coverage. Treatment rate scales up over period T_INTERVENTION_START to T_INTERVENTION_END
            iscenario_TreatWHOtarget = 7; %% Treatment reach 80% coverage. Treatment rate scales up over period T_INTERVENTION_START to T_INTERVENTION_END
            iscenario_WHOtarget = 8;     %% BD+HepB3 reach 90% coverage by T_INTERVENTION_END, treatment reaches 80% by T_INTERVENTION_END
            iscenario_MAP = 9; %% WUENIC 2025 BD+HepB3, 2016 treatment, Microarray patch introduced in T_INTERVENTION_START (increase BD coverage, but lower efficacy).
            iscenario_CPAD = 10; %% WUENIC 2025 BD+HepB3, 2016 treatment, CPAD patch introduced (increase BD but lower eff and different cost to MAP).
            iscenario_BD2025_birthcohorttest = 11;   %% WUENIC 2025 BD+HepB3, 2016 treatment, Thai-B-type testing of pre-BD birth cohort on top of existing testing (cap so cannot test >100% of any age stratum).
            iscenario_PAP_TREAThighVL = 12;   %% PAP for High VL pregnant women
            iscenario_PAP_TREATeAgpos = 13;             %% eAg+
            iscenario_PAP_TREAT_highVL_or_eAgpos = 14;  %% Either high VL or eAg+ (or both)
            iscenario_BASE2020_WITHTREAT = 15; %% The 'default' scenario - WUENIC 2019 BD+HepB3, 2016 treatment, no new interventions.
            iscenario_BASE2025_WITHTREAT = 16;     %% WUENIC 2025 BD+HepB3, 2016 treatment, no new interventions. Addresses - how have changes in BD+Hep B3 coverage impacted result?

            %%iscenario_BD2025_LA_TDF = 6; %% WUENIC 2025 BD+HepB3, 2016 treatment, long-acting treatment introduced (increases coverage of TDF treatment).
            %%iscenario_BD2025_PoC_ALT_HBcrAg = 7;    %% WUENIC 2025 BD+HepB3, 2016 treatment, PoC ALT and HBcrAg introduced - higher treatment coverage, also some people on treatment who don't need it.
            %%iscenario_BD2025_cure = 9; %% WUENIC 2025 BD+HepB3, 2016 treatment, (hypothetical) cure replaces treatment at current test rates.
            
        
            %% Index values for scenario_BD: Governs BD coverage time trends, introduction of different BD devices (MAP, CPAD).
            I_BD_WUENIC2020 = 1;  %% Follow WUENIC2020 (existing scenario) and after 2019 coverage remains at last (2019) value
            I_BD_WUENIC2025 = 2;  %% Follow WUENIC2025 and after 2024 coverage remains at last (2024) value
            I_BD_INFACILITY_INTRODUCTION = 3;  %% Follow WUENIC2025. For countries without BD, introduce BD in T_INTERVENTION_START and scale up to be some % of the in-facility births. For countries with BD already this will be identical to I_BD_WUENIC2025.
            I_BD_MAP = 4;  %% Follow WUENIC2025, then an extra (different efficacy) product increases overall BD coverage up to a level capped by out-of-facility deliveries.
            I_BD_CPAD = 5; %% Follow WUENIC2025, then an extra (different efficacy) product increases overall BD coverage up to a level capped by out-of-facility deliveries.
            I_BD_WHOtarget = 6; %% 90% BD coverage by T_INTERVENTION_END.

            %% Index values for scenario_HepB3: Hep B3 scenarios. Currently just use 2020 or 2025 WUENIC data.
            I_HEPB3_WUENIC2020 = 1;
            I_HEPB3_WUENIC2025 = 2;
            I_HEPB3_WHOtarget = 3; %% 90% coverage.
            
            %% Index values for scenario_PAP: peripartum antiviral prophylaxis (PAP) treatment for HBsAg+ mothers (treat all, treat high VL etc).
            %%I_PAP_SQ = 0;      %% No PAP unless already present.
            I_PAP_NOTREAT = 1;
            %% Current WHO recommendation is PAP for women HBsAg+ with HBV DNA levels >=200,000 IU/ml or if HBeAg+. 
            I_PAP_TREAThighVL = 2;             %% High VL
            I_PAP_TREATeAgpos = 3;             %% eAg+
            I_PAP_TREAT_highVL_or_eAgpos = 4;  %% Either high VL or eAg+ (or both)

            %% scenario_Treatment: governs how treatment happens:
            % Modified by treat-all, introduction of PoC HBcrAg, PoC ALT tests, 
            I_TREAT_INIT_SQ = 1;         %% Use current rates of treatment uptake and failure.
            I_TREAT_WHOtarget = 2;       %% rate to get 80% of eligibles on treatment (TDF).
            I_TREAT_40percent = 3;       %% rate to get 40% of eligibles on treatment (TDF).
            I_NOTREAT = 4;       %% no treatment ever.
            %%I_TREAT_INIT_POC_cr_ALT = 4; %% Introduce PoC tests for HBcrAg and ALT - increase rate of treatment initiation in eligible groups.
            %%I_TREAT_INIT_LA = 5;         % Introduce long-acting treatment. SQ treatment failure rate is very low (0.001), so we take the "TDF treatment" group to be "On treatment, adherent and not going to drop out". LA treatment then just increases the proportion of people in this compartment (either by improving adherence, preventing dropout, or offering a more convenient/preffered option).
            %%I_CURE = 6;                  % Cure replaces treatment - in this case 


            % %% scenario_TreatmentRetention: governs introduction of long-acting treatment (TDF).
            % I_TREAT_RET_SQ = 1;
            % I_TREAT_RET_LA_TREAT = 2;
            
            
            %% scenario_CohortTesting: governs whether we have testing of a birth cohort (like Thai-B). 
            I_NO_COHORT_TEST = 1;
            I_COHORT_TEST = 2;
       
            %% For each scenario we determine which intervention "levers" are used.
            switch scenario_num
                case iscenario_BASE2020notreat   %%case 'Status quo infant & BD'
                    scenario_BD = I_BD_WUENIC2020;
                    scenario_HepB3 = I_HEPB3_WUENIC2020;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_NOTREAT;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_BASE2025notreat     %% WUENIC 2025 BD+HepB3, 2016 treatment, no new interventions. Addresses - how have changes in BD+Hep B3 coverage impacted result?
                    scenario_BD = I_BD_WUENIC2025;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_NOTREAT;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_INFACILITYBD     %% WUENIC 2025 HepB3, 2016 treatment, BD introduced in countries where it is not already present - coverage capped at in-facility birth coverage.
                    scenario_BD = I_BD_INFACILITY_INTRODUCTION;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_NOTREAT;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_BDWHOtarget    %% BD reaches 90% coverage target
                    disp("BD target scenario")
                    scenario_BD = I_BD_WHOtarget;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_NOTREAT;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_HepB3WHOtarget  %% HepB3 reaches 90% coverage target
                    disp("HepB3 target scenario")
                    scenario_BD = I_BD_WUENIC2025;
                    scenario_HepB3 = I_HEPB3_WHOtarget;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_NOTREAT;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_Treat40percent  %% Treatment reaches 40% (half of WHO target)
                    disp("Treatment target scenario")
                    scenario_BD = I_BD_WUENIC2025;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_TREAT_40percent;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_TreatWHOtarget  %% Treatment reaches 80% target
                    disp("Treatment target scenario")
                    scenario_BD = I_BD_WUENIC2025;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_TREAT_WHOtarget;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_WHOtarget       %% BD, HepB3 and treatment reach targets
                    disp("WHO target scenario")
                    scenario_BD = I_BD_WHOtarget;
                    scenario_HepB3 = I_HEPB3_WHOtarget;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_TREAT_WHOtarget;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_MAP %% WUENIC 2025 BD+HepB3, 2016 treatment, Microarray patch introduced in T_INTERVENTION_START (increase BD coverage, but lower efficacy).
                    scenario_BD = I_BD_MAP;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_NOTREAT;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_CPAD %% WUENIC 2025 BD+HepB3, 2016 treatment, CPAD patch introduced (increase BD but lower eff and different cost to MAP).
                    scenario_BD = I_BD_CPAD;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_NOTREAT;
                case iscenario_BD2025_birthcohorttest   %% WUENIC 2025 BD+HepB3, 2016 treatment, Thai-B-type testing of pre-BD birth cohort on top of existing testing (cap so cannot test >100% of any age stratum).
                    scenario_BD = I_BD_WUENIC2025;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_NOTREAT;
                    scenario_CohortTesting = I_COHORT_TEST;     
                case iscenario_PAP_TREAThighVL
                    scenario_BD = I_BD_WUENIC2025;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_TREAThighVL;
                    scenario_Treatment = I_NOTREAT;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_PAP_TREATeAgpos
                    scenario_BD = I_BD_WUENIC2025;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_TREATeAgpos;
                    scenario_Treatment = I_NOTREAT;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_PAP_TREAT_highVL_or_eAgpos
                    scenario_BD = I_BD_WUENIC2025;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_TREAT_highVL_or_eAgpos;
                    scenario_Treatment = I_NOTREAT;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_BASE2020_WITHTREAT   %%case 'Status quo infant & BD'
                    scenario_BD = I_BD_WUENIC2020;
                    scenario_HepB3 = I_HEPB3_WUENIC2020;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_TREAT_INIT_SQ;
                    scenario_CohortTesting = I_NO_COHORT_TEST;
                case iscenario_BASE2025_WITHTREAT     %% WUENIC 2025 BD+HepB3, 2016 treatment, no new interventions. Addresses - how have changes in BD+Hep B3 coverage impacted result?
                    scenario_BD = I_BD_WUENIC2025;
                    scenario_HepB3 = I_HEPB3_WUENIC2025;
                    scenario_PAP = I_PAP_NOTREAT;
                    scenario_Treatment = I_TREAT_INIT_SQ;
                    scenario_CohortTesting = I_NO_COHORT_TEST;

                %     scenario_CohortTesting = I_NO_COHORT_TEST;
                % case iscenario_BD2025_LA_TDF %% WUENIC 2025 BD+HepB3, 2016 treatment, long-acting treatment introduced (increases coverage of TDF treatment).
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_NOTREAT;
                %     scenario_Treatment = I_TREAT_INIT_LA;
                %     scenario_CohortTesting = I_NO_COHORT_TEST;
                % case iscenario_BD2025_PoC_ALT_HBcrAg   %% WUENIC 2025 BD+HepB3, 2016 treatment, PoC ALT and HBcrAg introduced - higher treatment coverage, also some people on treatment who don't need it.
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_NOTREAT;
                %     scenario_Treatment = I_TREAT_INIT_POC_cr_ALT;
                %     scenario_CohortTesting = I_NO_COHORT_TEST;
                % case iscenario_BD2025_cure
                %     scenario_BD = I_BD_WUENIC2025;
                %     scenario_HepB3 = I_HEPB3_WUENIC2025;
                %     scenario_PAP = I_PAP_NOTREAT;
                %     scenario_Treatment = I_CURE;
                %     scenario_CohortTesting = I_NO_COHORT_TEST;
                otherwise
                    disp("Error - unknown scenario. Exiting")
                    return  %% Exit the script.
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% PARAMETERS FOR BD SCENARIOS:
            %% Here we specify what the maximum increase in BD coverage from MAP/CPAD would be:
            prop_accept_MAP = 0.9;   %% Placeholder assumption - easier to accept a patch than a needle
            prop_accept_CPAD = 0.85; %% Placeholder assumption

            %% Out-of-facility births: %% Placeholders
            % if(strcmp(ISO,"GMB"))
            %     prop_OOF_births = 0.837;
            % elseif(strcmp(ISO,"ETH"))
            %     prop_OOF_births = 0.475;
            % else 
            %     prop_OOF_births = 0.95; %% Placeholder - update to something better... 
            % end
            in_facility_BD_acceptance = 0.9;  %% For now assume up to 90% of in-facility births will give BD.
            prop_OOF_births = 1-GHO_infacilitybirthproportion_map(ISO);

            % Make model scenario birth dose coverage over time:
            %%last_BD_scaleup_year = 2030.0;  %% Assumption that BD plateaus after this time.
            
            switch scenario_BD
                case I_BD_WUENIC2020   
                    disp("I_BD_WUENIC2020")
                    coverage_BD_to_last_datapoint = BirthDose_wuenic2020;
                    year_last_BD_data = 2019;
                    %% Previously called 'Status quo infant & BD'. Follow WUENIC2020 (existing scenario) and after 2019 coverage remains at last (2019) value
                    future_xvals_vec = [2019.0, 2020.0, end_year];
                    future_yvals_vec = [BirthDose_wuenic2020(end), BirthDose_wuenic2020(end), BirthDose_wuenic2020(end)];
                    %% No MAP or CPAD introduced:
                    scenario_BDcoverage_fromMAP = zeros(1,length(years_vec_01yr));
                    scenario_BDcoverage_fromCPAD = zeros(1,length(years_vec_01yr));
                    %%label_array{scenario_num} = 'Status quo HepB3 & Hep2B-BD (baseline)';
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
                case I_BD_INFACILITY_INTRODUCTION
                    year_last_BD_data = 2024;
                    disp("I_BD_INFACILITY_INTRODUCTION")
                    coverage_BD_to_last_datapoint = BirthDose_wuenic2025;
                    %% Coverage up to 90% of the in-facility births, or current (2025 WUENIC) value - whichever is bigger.
                    max_in_facility_coverage = max(GHO_infacilitybirthproportion_map(ISO)*in_facility_BD_acceptance,BirthDose_wuenic2025(end));
                    future_xvals_vec = [2024.0, T_INTERVENTION_START, T_INTERVENTION_END, end_year];
                    future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), max_in_facility_coverage, max_in_facility_coverage];
                    %% No MAP or CPAD introduced:
                    scenario_BDcoverage_fromMAP = zeros(1,length(years_vec_01yr));
                    scenario_BDcoverage_fromCPAD = zeros(1,length(years_vec_01yr));
                case I_BD_MAP  %% MAP introduced:
                    disp("I_BD_MAP")
                    year_last_BD_data = 2024;
                    %% Follow WUENIC2025, then an extra (different efficacy) product increases overall BD coverage up to a level capped by out-of-facility deliveries.
                    coverage_BD_to_last_datapoint = BirthDose_wuenic2025;
                    %% This governs the coverage of the standard BD injection:
                    future_xvals_vec = [2024.0, 2025.0, end_year];  
                    future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), BirthDose_wuenic2025(end)];
                    %% This governs the coverage of the additional MAP injection:
                    future_xvals_vec_MAP = [2024.0, T_INTERVENTION_START, T_INTERVENTION_END, end_year];
                    %% Increase in BD if introduce MAP (requires BD to currently be available):
                    if(BirthDose_wuenic2025(end)>0)  %% BD already available
                        %% Increase in BD is capped to not exceed the current proportion not getting BD.
                        current_prop_not_getting_BD = (1-BirthDose_wuenic2025(end));
                        BD_increase_from_MAP  = min(prop_OOF_births*prop_accept_MAP, current_prop_not_getting_BD);
                    else                             %% BD not currently available
                        BD_increase_from_MAP = 0;
                    end

                    future_yvals_vec_MAP = [0, 0, BD_increase_from_MAP, BD_increase_from_MAP];
                    coverageMAP_to_present = zeros(1,length(BirthDose_wuenic2025));  
                    disp("MAP1")
                    scenario_BDcoverage_fromMAP = make_coverage_vec(start_year,num_year_divisions,dt,end_year,coverageMAP_to_present,future_xvals_vec_MAP, future_yvals_vec_MAP, year_last_BD_data);
                    scenario_BDcoverage_fromCPAD = zeros(1,length(years_vec_01yr));
                case I_BD_CPAD  %% CPAD introduced:
                    disp("I_BD_CPAD")
                    year_last_BD_data = 2024;
                    %% Follow WUENIC2025, then an extra (different efficacy) product increases overall BD coverage up to a level capped by out-of-facility deliveries.
                    coverage_BD_to_last_datapoint = BirthDose_wuenic2025;
                    %% This governs the coverage of the standard BD injection:
                    future_xvals_vec = [2024.0, 2025.0, end_year];  
                    future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), BirthDose_wuenic2025(end)];
                    
                    %% This governs the coverage of the additional CPAD injection:
                    future_xvals_vec_CPAD = [2024.0, T_INTERVENTION_START, T_INTERVENTION_END, end_year];
                    %% Increase in BD if introduce CPAD (requires BD to currently be available):
                    if(BirthDose_wuenic2025(end)>0)  %% BD already available
                        %% Increase in BD is capped to not exceed the current proportion not getting BD.
                        current_prop_not_getting_BD = (1-BirthDose_wuenic2025(end));
                        BD_increase_from_CPAD = min(prop_OOF_births*prop_accept_CPAD, current_prop_not_getting_BD);
                    else                             %% BD not currently available
                        BD_increase_from_CPAD = 0;
                    end

                    future_yvals_vec_CPAD = [0, 0, BD_increase_from_CPAD, BD_increase_from_CPAD];                    
                    coverageCPAD_to_present = zeros(1,length(BirthDose_wuenic2025));  
                    disp("CPAD1")                    
                    scenario_BDcoverage_fromCPAD = make_coverage_vec(start_year, num_year_divisions, dt, end_year, coverageCPAD_to_present, future_xvals_vec_CPAD, future_yvals_vec_CPAD, year_last_BD_data);
                    scenario_BDcoverage_fromMAP = zeros(1,length(years_vec_01yr));
                case I_BD_WHOtarget
                    disp("I_BD_WHOtarget")
                    year_last_BD_data = 2024;
                    %% Update using WUENIC 2025: follow WUENIC2025 and after 2024 coverage remains at last (2024) value
                    coverage_BD_to_last_datapoint = BirthDose_wuenic2025;
                    future_xvals_vec = [2024.0, T_INTERVENTION_START, T_INTERVENTION_END, end_year];
                    future_yvals_vec = [BirthDose_wuenic2025(end), BirthDose_wuenic2025(end), 0.9, 0.9];
                    %% No MAP or CPAD introduced:
                    scenario_BDcoverage_fromMAP = zeros(1,length(years_vec_01yr));
                    scenario_BDcoverage_fromCPAD = zeros(1,length(years_vec_01yr));
                otherwise 
                    disp("Error - unknown scenario_BD. Exiting")
                    return  %% Exit the script.

            end  %% end switch scenario_BD
            
            %% Now get the full timetrend of BD coverage from start_year to end_year (note that MAP/CPAD coverage is stored separately in scenario_BDcoverage_fromMAP/CPAD)
            disp("BD1")
            scenario_BDcoverage = make_coverage_vec(start_year,num_year_divisions,dt,end_year,coverage_BD_to_last_datapoint,future_xvals_vec,future_yvals_vec,year_last_BD_data);
            disp("BD2")
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
            


            %     %%case {'Status quo infant & BD expansion to 25%','Status quo infant & BD expansion to 50%','Status quo infant & BD expansion to 75%','Status quo infant & BD expansion to 90%'}
            %         X = extract_percent_from_BDexpansion_scenario_label(scenario);
            %         % X in {25%, 50%, 75%, 90%}
            %         % Any country that is <X%, goes to X% by last_BD_scaleup_year (linear expansion); others that are >X% stay at that level.
            %         max_BD_cov_val = max([BirthDose_wuenic2020(end) X/100]);
            %         future_xvals_vec = [2019.0 first_expansion_year last_BD_scaleup_year end_year];
            %         future_yvals_vec = [BirthDose_wuenic2020(end) BirthDose_wuenic2020(end) max_BD_cov_val max_BD_cov_val];
            % 
            %         label_array_string = replace(scenario, 'Status quo infant & BD expansion', 'HepB-BD scale-up');
            %         label_array{scenario_num} = label_array_string;
            % 
            %     case {'Status quo infant & BD drop 5 2020', 'Status quo infant & BD drop 10 2020', 'Status quo infant & BD drop 15 2020', 'Status quo infant & BD drop 20 2020'}
            %         % Follows status quo but, during the year 2020, birth dose vaccination drops by X% in relative terms.
            %         % X in {5%, 10%, 15%, 20%}
            %         X = extract_percent_from_BDdrop_scenario_label(scenario);
            % 
            %         future_xvals_vec = [2019.0 2019.9 2020.0 2020.9 2021.0 2030 end_year];
            %         %% Reduce BD by X/100, so the remaining BD coverage is (1-X/100) times original coverage.
            %         pc_reduction = 1 - X/100;
            %         future_yvals_vec = [BirthDose_wuenic2020(end) BirthDose_wuenic2020(end) pc_reduction*BirthDose_wuenic2020(end) pc_reduction*BirthDose_wuenic2020(end) BirthDose_wuenic2020(end) BirthDose_wuenic2020(end) BirthDose_wuenic2020(end)];
            % 
            %         %% label_array_string should look like 'HepB-BD disruptions 5% in 2020'.
            %         label_array_string = replace(replace(scenario, 'Status quo infant & BD drop', 'HepB-BD disruptions'),' 2020','% in 2020');
            %         label_array{scenario_num} = label_array_string;
            %     case {'Status quo infant & BD delayed expansion 2023 to 2030','Status quo infant & BD delayed expansion 2023 to 2033','Status quo infant & BD delayed expansion 2025 to 2040'}
            % 
            %         max_BD_cov_val_90 = max([BirthDose_wuenic2020(end) 0.9]);                
            % 
            %         if(scenario=='Status quo infant & BD delayed expansion 2023 to 2030')
            %             % Planned expansion of birth-dose vaccination is
            %             % postponed by 3 years, finishing in 2030
            %             label_array{scenario_num} = 'HepB-BD delayed & fast scale-up 2023 to 2030';
            %             start_delay = 3;   % Delay in starting (after 2020)
            %             end_delay = 0;     % Delay in finishing (after 2030)
            %         elseif(scenario=='Status quo infant & BD delayed expansion 2023 to 2033')
            %             % Planned expansion of birth-dose vaccination is
            %             % postponed by 3 years, still taking 10 years
            %             label_array{scenario_num} = 'HepB-BD delayed & normal scale-up 2023 to 2033';
            %             start_delay = 3;   % Delay in starting (after 2020)
            %             end_delay = 3;     % Delay in finishing (after 2030)
            %         elseif(scenario=='Status quo infant & BD delayed expansion 2025 to 2040')
            %             % Planned expansion of birth-dose vaccination is postponed by 5 years and takes longer
            %             label_array{scenario_num} = 'HepB-BD delayed & slow scale-up 2025 to 2040';
            %             start_delay = 5;   % Delay in starting (after 2020)
            %             end_delay = 10;     % Delay in finishing (after 2030)
            %         else
            %             fprintf('Unknown scenario %s. Exiting',scenario)
            %             return
            %         end
            % 
            %         future_xvals_vec = [2019.0 first_expansion_year first_expansion_year+start_delay last_BD_scaleup_year+end_delay end_year];
            %         future_yvals_vec = [BirthDose_wuenic2020(end) BirthDose_wuenic2020(end) BirthDose_wuenic2020(end) max_BD_cov_val_90 max_BD_cov_val_90];                  
            % end
            % switch scenario_BD
            %     case I_BD_WUENIC2020
            %         future_xvals_vec = [2019.0 first_expansion_year (first_expansion_year+0.1) end_year];
            %     future_yvals_vec = [HepB3_wuenic2020(end) HepB3_wuenic2020(end) 1 1];
            % 


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Now HepB3 coverage scenarios:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            switch scenario_HepB3
                case I_HEPB3_WUENIC2020
                    disp("I_HEPB3_WUENIC2020")
                    year_last_HepB3_data = 2019;
                    %% Previously called 'Status quo infant & BD'. Follow WUENIC2020 (existing scenario) and after 2019 coverage remains at last (2019) value
                    coverage_HepB3_to_last_datapoint = HepB3_wuenic2020;
                    future_xvals_vec = [2019.0 2020.0 end_year];
                    future_yvals_vec = [HepB3_wuenic2020(end) HepB3_wuenic2020(end) HepB3_wuenic2020(end)];
                    %%label_array{scenario_num} = 'Status quo HepB3 & Hep2B-BD (baseline)';
                case I_HEPB3_WUENIC2025
                    disp("I_HEPB3_WUENIC2025")
                    year_last_HepB3_data = 2024;
                    coverage_HepB3_to_last_datapoint = HepB3_wuenic2025;
                    future_xvals_vec = [2024.0, 2025.0, end_year];
                    future_yvals_vec = [HepB3_wuenic2025(end), HepB3_wuenic2025(end), HepB3_wuenic2025(end)];
                case I_HEPB3_WHOtarget
                    disp("I_HEPB3_WHOtarget")
                    year_last_HepB3_data = 2024;
                    coverage_HepB3_to_last_datapoint = HepB3_wuenic2025;
                    future_xvals_vec = [2024.0, T_INTERVENTION_START T_INTERVENTION_END, end_year];
                    future_yvals_vec = [HepB3_wuenic2025(end), HepB3_wuenic2025(end), 0.9 0.9];
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
            
            %% Dead code:
            %% Now copy infant_vacc_vec into params (from now on we use params.InfantVacc):
            %% params.InfantVacc = infant_vacc_vec;
                    
            
 
            %% Now make the BirthDose coverage vector and store it in params:
            % params.scenario_BirthDose_coverage = make_coverage_vec(start_year,num_year_divisions,dt,end_year,BirthDose_wuenic2020,future_xvals_vec,future_yvals_vec);
            % params.scenario_BirthDose_coverage = min(1,params.scenario_BirthDose_coverage);
            % assert(isequal(size(params.scenario_BirthDose_coverage),size(years_vec_01yr)))

            %% Of women who received BD, the % of PAP-eligible who get PAP:
            PAP_coverage_withBD = 1.0;
            %% Of women who DID NOT received BD, the % of PAP-eligible who get PAP:
            PAP_coverage_withoutBD = 1.0;
            PAP_scaleup_end_year = 2030;
            PAP_cov_params = struct('max_cov_BDandPAP_EAgHighVL', 0,...
                    'max_cov_BDandPAP_SAgHighVL', 0,...
                    'max_cov_BDandPAP_EAgLowVL', 0,...
                    'max_cov_BDandPAP_SAgLowVL', 0,...
                    'max_cov_PAPonly_EAgHighVL', 0,...
                    'max_cov_PAPonly_SAgHighVL', 0,...
                    'max_cov_PAPonly_EAgLowVL', 0,...
                    'max_cov_PAPonly_SAgLowVL', 0,...
                    'TScaleup_PAP', 2020);
            switch scenario_PAP 
                case I_PAP_NOTREAT
                    % No PAP
                    PAP_cov_params.max_cov_BDandPAP_EAgHighVL = 0;
                    PAP_cov_params.max_cov_BDandPAP_SAgHighVL = 0;
                    PAP_cov_params.max_cov_BDandPAP_EAgLowVL  = 0;
                    PAP_cov_params.max_cov_BDandPAP_SAgLowVL  = 0;
                    % Fraction of those who do not get BD that do get PAP:
                    PAP_cov_params.max_cov_PAPonly_EAgHighVL = 0;   
                    PAP_cov_params.max_cov_PAPonly_SAgHighVL = 0;
                    PAP_cov_params.max_cov_PAPonly_EAgLowVL  = 0;
                    PAP_cov_params.max_cov_PAPonly_SAgLowVL  = 0;
                    PAP_cov_params.TScaleup_PAP = PAP_scaleup_end_year;
                case I_PAP_TREAThighVL
                    % PAP to (PAP_coverage)% of the pregnant women who have High VL (whatever eAg status) (irrespective of BD)
                    PAP_cov_params.max_cov_BDandPAP_EAgHighVL = PAP_coverage_withBD;
                    PAP_cov_params.max_cov_BDandPAP_SAgHighVL = PAP_coverage_withBD;
                    PAP_cov_params.max_cov_BDandPAP_EAgLowVL  = 0;  
                    PAP_cov_params.max_cov_BDandPAP_SAgLowVL  = 0;
                    PAP_cov_params.max_cov_PAPonly_EAgHighVL = PAP_coverage_withoutBD;   % coverage among those who missed BD
                    PAP_cov_params.max_cov_PAPonly_SAgHighVL = PAP_coverage_withoutBD;
                    PAP_cov_params.max_cov_PAPonly_EAgLowVL  = 0;
                    PAP_cov_params.max_cov_PAPonly_SAgLowVL  = 0;
                    PAP_cov_params.TScaleup_PAP = PAP_scaleup_end_year;
                case I_PAP_TREATeAgpos
                    % PAP to (PAP_coverage)% of the pregnant women who are HBeAg+ (whatever VL status) (irrespective of BD)
                    PAP_cov_params.max_cov_BDandPAP_EAgHighVL = PAP_coverage_withBD;
                    PAP_cov_params.max_cov_BDandPAP_SAgHighVL = 0;
                    PAP_cov_params.max_cov_BDandPAP_EAgLowVL = PAP_coverage_withBD;  
                    PAP_cov_params.max_cov_BDandPAP_SAgLowVL = 0;
                    PAP_cov_params.max_cov_PAPonly_EAgHighVL = PAP_coverage_withoutBD;   % coverage among those who missed BD
                    PAP_cov_params.max_cov_PAPonly_SAgHighVL = 0;
                    PAP_cov_params.max_cov_PAPonly_EAgLowVL  = PAP_coverage_withoutBD;
                    PAP_cov_params.max_cov_PAPonly_SAgLowVL  = 0;
                    PAP_cov_params.TScaleup_PAP = PAP_scaleup_end_year;
                case I_PAP_TREAT_highVL_or_eAgpos
                    % PAP to (PAP_coverage)% of the pregnant women who are HBeAg+ or high VL (or both) (irrespective of BD)
                    PAP_cov_params.max_cov_BDandPAP_EAgHighVL = PAP_coverage_withBD;
                    PAP_cov_params.max_cov_BDandPAP_SAgHighVL = PAP_coverage_withBD;
                    PAP_cov_params.max_cov_BDandPAP_EAgLowVL = PAP_coverage_withBD;  
                    PAP_cov_params.max_cov_BDandPAP_SAgLowVL = 0;
                    PAP_cov_params.max_cov_PAPonly_EAgHighVL = PAP_coverage_withoutBD;   % coverage among those who missed BD
                    PAP_cov_params.max_cov_PAPonly_SAgHighVL = PAP_coverage_withoutBD;
                    PAP_cov_params.max_cov_PAPonly_EAgLowVL  = PAP_coverage_withoutBD;
                    PAP_cov_params.max_cov_PAPonly_SAgLowVL  = 0;
                    PAP_cov_params.TScaleup_PAP = PAP_scaleup_end_year;
                case I_PAP_TREAT_sAgpos
                    % PAP to (PAP_coverage)% of the pregnant women who are HBeAg+ or high VL (or both) (irrespective of BD)
                    PAP_cov_params.max_cov_BDandPAP_EAgHighVL = PAP_coverage_withBD;
                    PAP_cov_params.max_cov_BDandPAP_SAgHighVL = PAP_coverage_withBD;
                    PAP_cov_params.max_cov_BDandPAP_EAgLowVL = PAP_coverage_withBD;
                    PAP_cov_params.max_cov_BDandPAP_SAgLowVL = PAP_coverage_withBD;
                    PAP_cov_params.max_cov_PAPonly_EAgHighVL = PAP_coverage_withoutBD;   % coverage among those who missed BD
                    PAP_cov_params.max_cov_PAPonly_SAgHighVL = PAP_coverage_withoutBD;
                    PAP_cov_params.max_cov_PAPonly_EAgLowVL  = PAP_coverage_withoutBD;
                    PAP_cov_params.max_cov_PAPonly_SAgLowVL  = PAP_coverage_withoutBD;
                    PAP_cov_params.TScaleup_PAP = PAP_scaleup_end_year;
                otherwise
                    disp("Error: Unknown value for scenario_PAP. Exiting")
                    return
            end

            
            treatment_rate_params = struct('Treatmentrate_2016', stochas_params_mat(stochas_run_num,country_start_col+7),...
                    'Treatmentrate_final', 0,...
                    't_treatment_scaleup_start', T_INTERVENTION_START,...
                    't_treatment_scaleup_end', T_INTERVENTION_END);
            switch scenario_Treatment
                case I_TREAT_INIT_SQ           %% Use current rates of treatment uptake and failure.
                    %%params.PriorTDFTreatRate = stochas_params_mat(stochas_run_num,country_start_col+7);
                    treatment_rate_params.Treatmentrate_final = treatment_rate_params.Treatmentrate_2016;
                    HAS_TREATMENT = 1; % Treatment introduced in 2016 at rate from MJdV.
                case I_TREAT_WHOtarget
                    treatment_rate_params.Treatmentrate_final = treatment_boundaries_vec(5); % 80% scenario from MJDV code
                    HAS_TREATMENT = 1;
                case I_TREAT_40percent
                    treatment_rate_params.Treatmentrate_final = treatment_boundaries_vec(3); % 40% scenario from MJDV code
                    HAS_TREATMENT = 1;
                case I_NOTREAT
                    treatment_rate_params.Treatmentrate_final = 0; % no treatment
                    HAS_TREATMENT = 0;  % Treatment never introduced
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
            % switch scenario_CohortTesting
            %     case I_NO_COHORT_TEST
            %         a=7;
            %     case I_COHORT_TEST
            %         a=8;
            %     otherwise
            %         disp("Error: Unknown value for scenario_CohortTesting. Exiting")
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
            start_year_vacc_coverage = 1890;
            last_year_run = end_year;
            PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgHighVL = PAP_coverage_scaleup(start_year_vacc_coverage, PAP_cov_params.TScaleup_PAP, ...
                last_year_run, dt, PAP_cov_params.max_cov_BDandPAP_EAgHighVL);

            PAP_cov_params.scenario_PAPcoverage_BDandPAP_EAgLowVL = PAP_coverage_scaleup(start_year_vacc_coverage, PAP_cov_params.TScaleup_PAP, ...
                last_year_run, dt, PAP_cov_params.max_cov_BDandPAP_EAgLowVL);

            PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgHighVL = PAP_coverage_scaleup(start_year_vacc_coverage, PAP_cov_params.TScaleup_PAP, ...
                last_year_run, dt, PAP_cov_params.max_cov_BDandPAP_SAgHighVL);
    
            PAP_cov_params.scenario_PAPcoverage_BDandPAP_SAgLowVL = PAP_coverage_scaleup(start_year_vacc_coverage, PAP_cov_params.TScaleup_PAP, ...
                last_year_run, dt, PAP_cov_params.max_cov_BDandPAP_SAgLowVL);

            % Coverage of PAP among those not with BD
            %%xvals_vec = [start_year_vacc_coverage TScaleup_PAP-1 TScaleup_PAP last_year_run];
            PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgHighVL = PAP_coverage_scaleup(start_year_vacc_coverage, PAP_cov_params.TScaleup_PAP, ...
                last_year_run, dt, PAP_cov_params.max_cov_PAPonly_EAgHighVL);

            PAP_cov_params.scenario_PAPcoverage_PAPonly_EAgLowVL = PAP_coverage_scaleup(start_year_vacc_coverage, PAP_cov_params.TScaleup_PAP, ...
                last_year_run, dt, PAP_cov_params.max_cov_PAPonly_EAgLowVL);

            PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgHighVL = PAP_coverage_scaleup(start_year_vacc_coverage, PAP_cov_params.TScaleup_PAP, ...
                last_year_run, dt, PAP_cov_params.max_cov_PAPonly_SAgHighVL);
   
            PAP_cov_params.scenario_PAPcoverage_PAPonly_SAgLowVL = PAP_coverage_scaleup(start_year_vacc_coverage, PAP_cov_params.TScaleup_PAP, ...
                last_year_run, dt, PAP_cov_params.max_cov_PAPonly_SAgLowVL);


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
                p_ChronicCarriage,Prog_scenario,Transactions,......
                scenario_BDcoverage, scenario_BDcoverage_fromMAP,...
                scenario_BDcoverage_fromCPAD, scenario_HepB3coverage, ...
                HAS_TREATMENT,...
                ISO, scenario_num, scenario_CohortTesting, ...
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
        save(fullfile(basedir,'outputs',filename_results),'outMap') 
        if strcmp(stochas_run_str,'1') && strcmp(sensitivity_analysis,'default')
            save(fullfile(basedir,'outputs','scenarios_array.mat'),'label_array') % only version saved after last scenario is correct
        end

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


function coverage = PAP_coverage_scaleup(start_year_vacc_coverage, TScaleup_PAP, ...
    last_year_run, dt, PAP_coverage_thissubgroup)
    xvals_vec = [start_year_vacc_coverage TScaleup_PAP-1 TScaleup_PAP last_year_run];
    % Scales up linearly from 0 to PAP_coverage_thissubgroup over the period
    % (TScaleup_PAP-1) to TScaleup_PAP
    yvals_vec = [0 0 PAP_coverage_thissubgroup PAP_coverage_thissubgroup];

    TimeSteps = start_year_vacc_coverage:dt:last_year_run; % 1 x 2101 double; [1890 1890.1 1890.2 ... 2099.8 2099.9 2100 2100.1 ... 2100.8 2100.9 2101]


    coverage = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
    % Ensure coverage is capped at 100%:
    coverage = min(1,coverage); 
end