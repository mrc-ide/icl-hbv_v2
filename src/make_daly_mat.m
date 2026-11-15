function outmat = make_daly_mat(model_run,num_years_simul,num_year_1980_2100,life_expectancy)


    yll_spread = squeeze(sum(model_run.Prev_Deaths_1yr(:,1:(life_expectancy-1),:),1)); % sum over gender
    % prevalence of deaths because one wants to count total deaths
    assert(isequal(size(yll_spread),[(life_expectancy-1) num_years_simul+1]))
    yll_spread = [yll_spread; zeros(100-(life_expectancy-1),num_years_simul+1)];
    assert(isequal(size(yll_spread),[100 num_years_simul+1]))
    
    %% MP: as stated below, the yld_1yr is calculated in HBVmodel.m as the DALY weight multiplied by the number of people in each compartment.
    %% Current DALY weights are for severe acute, decompensated cirrhosis and HCC only (death is dealth with separately):
    %% weights: susc=0, ImmTol=0, ImmReact=0, AsymptCarr=0, Chronic=0,
    % CompCirr=0, DecompCirr=0.178, HCC=0.54, Immune=0, TDFtreat_LEGACY=0, death(dealt with separately)=0,
    % i_funct_cure_path=0, i_funct_cure=0, NonSevAcute=0, SevereAcute=0.133
        
    yld_spread = squeeze(sum(model_run.yld_1yr,1)); % years living with the disease
    % Wikipedia: YLD = I x DW x L, where I = number of incident cases in the population, DW = disability weight of specific condition, and L = average duration of the case until remission or death (years)
    % I x L = P therefore YLD = P x DW, which is done in the model
    % model_run.yld_1yr is a 2 x 100 x num_year_1980_2100 matrix 
    % sum over genders to get a 100 x 101 matrix of 100 age groups (0 to 99) versus 101 years (2000 to 2100)
    assert(isequal(size(yld_spread),[100 num_years_simul+1]))

    outmat = yll_spread + yld_spread;
    % DALY = YLL + YLD
end % end function make_daly_mat

