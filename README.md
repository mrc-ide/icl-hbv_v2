# icl-hbv

This is an updated version of the HBV model developed at Imperial College London, developed from the version described in the paper ["The impact of the timely birth dose vaccine on the global elimination of hepatitis B"](https://www.nature.com/articles/s41467-021-26475-6).

All of the scripts have been run in MATLAB (version R2025b).

In the folder `src`, the main files used are HBVmodel.m, country_level_analyses.m and main_script.m. 

`main_script.m` runs `country_level_analyses.m`, which in turn runs `HBVmodel.m` (the latter is in the folder `model`). These can generate up to 200 stochastic runs per country, across up to 110 countries, with multiple scenarios and sensitivity analyses.

The following 
<ol> 
  <li> </li>num_stochas_runs </li> 
<li> num_scenarios - determines number of scenarios to run (from ListOfScenarios) </li>
<li> i_start_country and i_end_country - indices for the start and end country to loop over. i_start_country=1 and i_end_country=110 would run all countries. i_start_country=1 and i_start_country=1 would run only the first country (Afghanistan). </li>
<li> sensitivity_analysis_list is a list of the sensitivity analyses to be run (for each stochastic run, country, and scenario). /li>
</ol>

At present there are additional Matlab scripts from the above paper which have not been updated (though should still run if the model is run for all countries and scenarios):  `summarize_stochastic_runs.m` (summarizes the results from the 800 results files into results for each country, each WHO region, as well as global results) and  `draw_figures.m`( draws Figures 1e, 1f, 2a, 2d, 3a, 3b, 4a, 4b, 5a, 5b, Supplementary figures 1, 2a, 2b and 3 and Supplementary Tables 3 and 4 (default analyses)).  Finally The folder `data_and_script_for_figures` contains the script `draw_figs.m` for drawing Figures 1f, 2b, 3a, 3b, 4a&ndash;4b, 5a&ndash;5b, Supplementary Figures 1, 2b and 3 and Supplementary Tables 3 and 4 (default analyses). Note that the Dropbox folder containing the data structure for running this code no longer exists, so for now `draw_figs.m` is unlikely to work.
