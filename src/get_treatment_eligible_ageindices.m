function treat_eligibility_ageindices = get_treatment_eligible_ageindices(scenario_treat_elig, current_natural_history_state, ...
    i_natural_hist, ages)

    num_age_steps = length(ages);

    %% Current treatment eligibility:
    if(strcmp(scenario_treat_elig,"Current treatment"))
        if(current_natural_history_state==i_natural_hist.ImmTol)
            i30y = find(ages >= 30, 1);       %% Age boundary for i_natural_hist.ImmTol state
            treat_eligibility_ageindices = i30y:num_age_steps;
        elseif(current_natural_history_state==i_natural_hist.ImmReact || current_natural_history_state==i_natural_hist.Chronic || ...
                current_natural_history_state==i_natural_hist.CompCirr || current_natural_history_state==i_natural_hist.DecompCirr)
            treat_eligibility_ageindices = 1:num_age_steps;
        else
            treat_eligibility_ageindices = [];
        end
    %% Universal treatment
    elseif(strcmp(scenario_treat_elig,"Universal treatment"))
        if(current_natural_history_state==i_natural_hist.ImmTol || current_natural_history_state==i_natural_hist.ImmReact || ...
                current_natural_history_state==i_natural_hist.AsymptCarr || current_natural_history_state==i_natural_hist.Chronic || ...
                current_natural_history_state==i_natural_hist.CompCirr || current_natural_history_state==i_natural_hist.DecompCirr)
            treat_eligibility_ageindices = 1:num_age_steps;
        else
            treat_eligibility_ageindices = [];
        end
    else
        disp("Error - unknown treatment eligbility scenario in get_treatment_eligible_ageindices(). Exiting")
        %% Hopefully returning without a value will cause an error.
        return
    end
end
