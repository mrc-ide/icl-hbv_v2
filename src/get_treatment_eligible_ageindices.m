function treat_eligibility_ageindices = get_treatment_eligible_ageindices(ages, current_natural_history_state, ...
    i_natural_hist)

    num_age_steps = length(ages);
    if(current_natural_history_state==i_natural_hist.ImmTol)
        i30y = find(ages >= 30, 1);       %% Age boundary for i_natural_hist.ImmTol state
        treat_eligibility_ageindices = i30y:num_age_steps;
    elseif(current_natural_history_state==i_natural_hist.ImmReact || current_natural_history_state==i_natural_hist.Chronic || ...
            current_natural_history_state==i_natural_hist.CompCirr || current_natural_history_state==i_natural_hist.DecompCirr)
        treat_eligibility_ageindices = 1:num_age_steps;
    else
        treat_eligibility_ageindices = [];
    end
end
