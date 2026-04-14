%% Function previously two functions "make_effparams" and "augment_effparams".
%% These were local functions in country_level_analyses.m for the PAP model 
%% The PAP model is the "mjdv/PAP_scripts" branch of the icl-hbv github repo.
%% That model is also in the subdirectory icl-hbv-mjdv-PAP_scripts/ on my laptop.
function PAP_VL_params = assign_PAP_VL_params(params)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Function makes a structure "PAP_VL_params" that contains parameters relevant to
    %% the PAP part of the model and/or the VL division (pregnancy transmission is split into a
    %% portion due to high VL and a portion due to low VL, so that interventions - including PAP - that
    %% target women with high VL can be modelled.
    %% Function takes the parameter params.p_VerticalTransmission_HbSAg_NoIntv
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PAP_VL_params.FracEPosHighVL = 0.90; % proportion of HBsAg+ HBeAg+ pregnant women that have high viral load
    PAP_VL_params.FracSPosHighVL = 0.05; % proportion of HBsAg+ HBeAg- pregnant women that have high viral load
    %% We need p_VerticalTransmission_HbSAg_NoIntv to be in params (it was used in the non-PAP code branch too), so keep it.
    p_HbSAg_av = params.p_VerticalTransmission_HbSAg_NoIntv;
    p_HbEAg_av = params.p_VerticalTransmission_HbEAg_NoIntv;
    %%p_HbEAg_av = 0.9; % p_VerticalTransmission_HbEAg_NoIntv


    % -------
    % ** WAY OF INTERPRETING AVAILABLE DATA #1: Transmission parameters option:
    % just differentiated by VL **
    % [Assume that transmission rate only depends on VL ..[ and make the rates 
    % by e/s-status and VL consistent with the given average of transmsision 
    % among e and s-status]

    % Given average of HBsAg is p_VerticalTransmission_HbSAg_NoIntv)
    % Given average of HBeAg is 0.90 (using in previous fitting)
    
    % get acceptable values for 'FracSPosHighVL'
    range_for_FracSPosHighVL = get_range_of_FracSPosHighVL(PAP_VL_params.FracEPosHighVL, p_HbSAg_av, p_HbEAg_av);

    assert(PAP_VL_params.FracSPosHighVL>=range_for_FracSPosHighVL(1))
    assert(PAP_VL_params.FracSPosHighVL<=range_for_FracSPosHighVL(2))

    % Work out the transmission rates for low and high VL (for this
    % case of no differentiation by E/S state).
    % Use the constraints of the average E and S rates, and the fact that
    % transmission rate must not vary within VL categories
    beta_low_vl = max(1e-10,( p_HbEAg_av * PAP_VL_params.FracSPosHighVL - p_HbSAg_av * PAP_VL_params.FracEPosHighVL ) / ( PAP_VL_params.FracSPosHighVL - PAP_VL_params.FracEPosHighVL ) );
    beta_high_vl = (p_HbSAg_av - (1-PAP_VL_params.FracSPosHighVL)*beta_low_vl)/PAP_VL_params.FracSPosHighVL;
    assert(beta_low_vl>0); assert(beta_low_vl<1); assert(beta_high_vl>0); assert(beta_high_vl<1); assert(beta_high_vl>beta_low_vl)
    
   

    PAP_VL_params.p_VertTrans_HbSAgLowVL_NoIntv  = beta_low_vl; 
    PAP_VL_params.p_VertTrans_HbSAgHighVL_NoIntv = min(1,beta_high_vl);
    PAP_VL_params.p_VertTrans_HbEAgLowVL_NoIntv  = beta_low_vl;
    PAP_VL_params.p_VertTrans_HbEAgHighVL_NoIntv = min(1,beta_high_vl);

    %% These were the previous versions (from Shevanthi's PAP paper) but they are mathematically identical.
    %%ratio_high_to_low = beta_high_vl  / beta_low_vl;
    %%PAP_VL_params.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL=ratio_high_to_low;
    %%PAP_VL_params.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL=ratio_high_to_low; 

    %%p_VerticalTransmission_HbSAgLowVL_NoIntv=p_HbSAg_av / ( (1-PAP_VL_params.FracSPosHighVL) + PAP_VL_params.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL*PAP_VL_params.FracSPosHighVL ); 
    %%p_VerticalTransmission_HbSAgHighVL_NoIntv=min(1,PAP_VL_params.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL * p_VerticalTransmission_HbSAgLowVL_NoIntv);
    
    %%p_VerticalTransmission_HbEAgLowVL_NoIntv=p_HbEAg_av / ( (1-PAP_VL_params.FracEPosHighVL) + PAP_VL_params.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL*PAP_VL_params.FracEPosHighVL ); 
    %%p_VerticalTransmission_HbEAgHighVL_NoIntv=min(1,PAP_VL_params.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL * p_VerticalTransmission_HbEAgLowVL_NoIntv);
    
    


    % average transmission rate among S+ correct:
    assert(abs(p_HbSAg_av - (PAP_VL_params.p_VertTrans_HbSAgLowVL_NoIntv*(1-PAP_VL_params.FracSPosHighVL) + PAP_VL_params.p_VertTrans_HbSAgHighVL_NoIntv*PAP_VL_params.FracSPosHighVL))<0.001)
    % p_HbSAg_av - (p_VerticalTransmission_HbSAgLowVL_NoIntv*(1-PAP_VL_params.FracSPosHighVL) + p_VerticalTransmission_HbSAgHighVL_NoIntv*PAP_VL_params.FracSPosHighVL)
    % = p_HbSAg_av - (p_VerticalTransmission_HbSAgLowVL_NoIntv*(1-PAP_VL_params.FracSPosHighVL) + PAP_VL_params.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL*p_VerticalTransmission_HbSAgLowVL_NoIntv*PAP_VL_params.FracSPosHighVL)
    % = p_HbSAg_av - (p_VerticalTransmission_HbSAgLowVL_NoIntv - p_VerticalTransmission_HbSAgLowVL_NoIntv*PAP_VL_params.FracSPosHighVL + PAP_VL_params.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL*p_VerticalTransmission_HbSAgLowVL_NoIntv*PAP_VL_params.FracSPosHighVL)
    % = p_HbSAg_av - (p_VerticalTransmission_HbSAgLowVL_NoIntv * (1 - PAP_VL_params.FracSPosHighVL + PAP_VL_params.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL*PAP_VL_params.FracSPosHighVL))
    % = p_HbSAg_av - p_HbSAg_av
    % = 0
    
    % average transmission rate among E+ correct:
    assert(abs(p_HbEAg_av - (PAP_VL_params.p_VertTrans_HbEAgLowVL_NoIntv*(1-PAP_VL_params.FracEPosHighVL) + PAP_VL_params.p_VertTrans_HbEAgHighVL_NoIntv*PAP_VL_params.FracEPosHighVL))<0.001)
    % p_HbEAg_av - (p_VerticalTransmission_HbEAgLowVL_NoIntv*(1-PAP_VL_params.FracEPosHighVL) + p_VerticalTransmission_HbEAgHighVL_NoIntv*PAP_VL_params.FracEPosHighVL)
    % = p_HbEAg_av - (p_VerticalTransmission_HbEAgLowVL_NoIntv*(1-PAP_VL_params.FracEPosHighVL) + PAP_VL_params.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL*p_VerticalTransmission_HbEAgLowVL_NoIntv*PAP_VL_params.FracEPosHighVL)
    % = p_HbEAg_av - (p_VerticalTransmission_HbEAgLowVL_NoIntv - p_VerticalTransmission_HbEAgLowVL_NoIntv*PAP_VL_params.FracEPosHighVL + PAP_VL_params.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL*p_VerticalTransmission_HbEAgLowVL_NoIntv*PAP_VL_params.FracEPosHighVL)
    % = p_HbEAg_av - (p_VerticalTransmission_HbEAgLowVL_NoIntv * (1 - PAP_VL_params.FracEPosHighVL + PAP_VL_params.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL*PAP_VL_params.FracEPosHighVL))
    % = p_HbEAg_av - p_HbEAg_av

    % equal transmssion rates within a VL category
    %%assert(abs(PAP_VL_params.p_VerticalTransmission_HbSAgLowVL_NoIntv-PAP_VL_params.p_VerticalTransmission_HbEAgLowVL_NoIntv)<0.001);

    %%assert(abs(PAP_VL_params.p_VerticalTransmission_HbSAgHighVL_NoIntv-PAP_VL_params.p_VerticalTransmission_HbEAgHighVL_NoIntv)<0.001);
    % p_VerticalTransmission_HbSAgHighVL_NoIntv-p_VerticalTransmission_HbEAgHighVL_NoIntv
    % = PAP_VL_params.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL * p_VerticalTransmission_HbSAgLowVL_NoIntv - (PAP_VL_params.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL * p_VerticalTransmission_HbEAgLowVL_NoIntv)
    % = ratio_high_to_low * p_VerticalTransmission_HbSAgLowVL_NoIntv - (ratio_high_to_low * p_VerticalTransmission_HbEAgLowVL_NoIntv)
    % = ratio_high_to_low * (p_VerticalTransmission_HbSAgLowVL_NoIntv - p_VerticalTransmission_HbEAgLowVL_NoIntv)
    % = ratio_high_to_low * 0 % from the previous assert
    % = 0



    %% End of code from make_effparams()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% Start of code from augment_effparams()
    %% Note that some of this  code was taken from country_level_analyses.m (in PAP script):
    % if do_geometric_median % if doing the geometric median, use default Efficacy_BirthDoseVacc_HbEAg and Efficacy_InfantVacc values
    %     demog.pr_VerticalTransmission_HbSAgLowVL_PAP = mean([0 0.05]); % 0.025
    % else % each of the 200 stochastic particles has its own Efficacy_BirthDoseVacc_HbEAg and Efficacy_InfantVacc values
    %     demog.pr_VerticalTransmission_HbSAgLowVL_PAP = stochas_params_mat(stochas_run_num,end);
    % end
    %% It was then moved into PAP_VL_params in the (now defunct) function "augment_effparams()".
    %% PAP_VL_params.pr_VerticalTransmission_HbSAgLowVL_PAP = demog.pr_VerticalTransmission_HbSAgLowVL_PAP; % probability ratio relative to the corresponding group's _NoIntv level; 2022/12/02
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% "pr_" and "pRatio_" both mean a probability ratio (pRatio is MP's updated version). "p_" means a probability.
    %% probability ratios relative to the corresponding (SAg/Eag and High/LowVL group) group's _NoIntv level; 2022/03/01
    %% Note that probability ratio = 1 - Effectiveness.
    
    %% Low VL, SAg+ (EAg-) mothers:
    PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_BD = 0.0;              % if birth dose 
    PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_PAP = mean([0 0.05]);  % PAP given; 0.025
    PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_CPAD = 0.01;            % if CPAD
    PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_MAP = 0.02;            % if MAP
    PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_BD_PAP = 0.0;          % BD and PAP
    PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_CPAD_PAP = 0.0;          % CPAD and PAP
    PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_MAP_PAP = 0.0;          % MAP and PAP
    
    %%PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_Treat = 0.0;           % mother on treatment (TDF/interferon)
    %%PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_BD_Treat = 0.0;        % BD and treatment
    %%PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_CPAD_Treat = 0.0;        % CPAD and treatment
    %%PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_MAP_Treat = 0.0;        % MAP and treatment
    
    %% High VL, SAg+ (EAg-) mothers:
    PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_BD = 0.20;            % BD for SAg+ (EAg-) high VL (relative to no intervention SAg+, EAg-, high VL group)
    PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_PAP = 0.06;           % probability ratio relative to the corresponding group's _NoIntv level
    PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_CPAD = 0.3;            % if CPAD
    PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_MAP = 0.35;            % if MAP
    PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_BD_PAP = 0.01;        % probability ratio relative to the corresponding group's _NoIntv level
    PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_CPAD_PAP = 0.02;          % CPAD and PAP
    PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_MAP_PAP = 0.02;          % MAP and PAP
    
    %%PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_Treat = 0.02;         % mother on treatment (TDF/interferon). Non-VL model has treatment 98% effective at stopping MTCT (regardless of SAg/EAg).
    %%PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_BD_Treat = 0.01;       % BD and treatment (mother on TDF/interferon)
    %%PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_CPAD_Treat = 0.02;        % CPAD and treatment
    %%PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_MAP_Treat = 0.02;        % MAP and treatment
    
    
    %% Low VL, EAg+ (also SAg+) mothers:
    PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_BD        = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_BD;
    PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_PAP       = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_PAP;
    PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_CPAD      = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_CPAD;
    PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_MAP       = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_MAP;
    PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_BD_PAP    = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_BD_PAP;
    PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_CPAD_PAP  = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_CPAD_PAP;
    PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_MAP_PAP   = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_MAP_PAP;
    
    %%PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_Treat     = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_Treat;
    %%PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_BD_Treat  = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_BD_Treat;
    %%PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_CPAD_Treat= PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_CPAD_Treat;
    %%PAP_VL_params.pRatio_VertTrans_HbEAgLowVL_MAP_Treat = PAP_VL_params.pRatio_VertTrans_HbSAgLowVL_MAP_Treat;
    
    %% High VL, EAg+ (also SAg+) mothers:
    PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_BD       = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_BD;
    PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_PAP      = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_PAP;
    PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_CPAD     = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_CPAD;
    PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_MAP      = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_MAP;
    PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_BD_PAP   = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_BD_PAP;
    PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_CPAD_PAP = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_CPAD_PAP;
    PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_MAP_PAP = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_MAP_PAP;
    
    %%PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_Treat    = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_Treat;
    %%PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_BD_Treat = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_BD_Treat;
    %%PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_CPAD_Treat = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_CPAD_Treat;
    %%PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_MAP_Treat = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_MAP_Treat;

    % mother on treatment (TDF/interferon). Non-VL model has treatment 98% effective at stopping MTCT (regardless of SAg/EAg).
    PAP_VL_params.pRatio_VertTrans_Treat      = 0.02;           
    PAP_VL_params.pRatio_VertTrans_Treat_BD   = 0.004;      % BD and treatment (mother on TDF/interferon). PLACEHOLDER - Ball-park BD ~80% effective so 0.02 * (1-0.8) = 0.004
    PAP_VL_params.pRatio_VertTrans_Treat_CPAD = 0.006;      % CPAD and treatment - PLACEHOLDER
    PAP_VL_params.pRatio_VertTrans_Treat_MAP  = 0.008;      % MAP and treatment  - PLACEHOLDER



    %PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_BD = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_BD;
    %PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_BD_PAP = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_BD_PAP;
    %PAP_VL_params.pRatio_VertTrans_HbEAgHighVL_PAP = PAP_VL_params.pRatio_VertTrans_HbSAgHighVL_PAP;

end % end function augment_PAP_VL_params
