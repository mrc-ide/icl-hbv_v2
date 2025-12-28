function write_parameter_txt_files(ListOfISOs, BD_table, HepB3_table, country_s_e_HCCdeaths_map, dwvec, ...
    params_map, stochas_params_mat, num_in_treatment_2016_map, treatment_rates_map, num_countries)


    allcountry_ISOs = "";
    for country_num = 1:num_countries
        ISO = ListOfISOs{country_num};
        %%thiscountry_ISO = jsonencode(ListOfISOs(ISO),PrettyPrint=true);
        thiscountry_ISO = string(country_num)+ ": "+ISO+newline;
        allcountry_ISOs = strcat(allcountry_ISOs,thiscountry_ISO);
    end
    writelines(allcountry_ISOs,"raw_params/ListOfISOs.txt");

    %% BD_table and HepB3_table are 194 (all countries including those not in the modelling exercise) x40 (years?) tables. 
    writetable(BD_table,"raw_params/BD_table.csv",'WriteRowNames',true); 
    writetable(HepB3_table,"raw_params/HepB3_table.csv",'WriteRowNames',true); 


    allcountry_hccdeaths = "";
    for country_num = 1:num_countries
        ISO = ListOfISOs{country_num};
        thiscountry_hccdeaths = jsonencode(country_s_e_HCCdeaths_map(ISO),PrettyPrint=true);
        allcountry_hccdeaths = strcat(allcountry_hccdeaths," COUNTRY: ",ISO,thiscountry_hccdeaths);
    end
    writelines(allcountry_hccdeaths,"raw_params/country_s_e_HCCdeaths.txt");

    

    writematrix(dwvec, "raw_params/dwvec.csv");
    allcountry_params = "";
    for country_num = 1:num_countries
        ISO = ListOfISOs{country_num};
        thiscountry_params = jsonencode(params_map(ISO),PrettyPrint=true);
        allcountry_params = allcountry_params + " COUNTRY: " + ISO + thiscountry_params + newline;
    end
    writelines(allcountry_params,"raw_params/params.txt");



    %% stochas_params_mat is a 200x883 matrix. I think the 200 refers to the number of stochastic runs.
    %% MP question: what is the 883? 883 is prime.
    writematrix(stochas_params_mat, "raw_params/stochas_params_mat.csv");



    allcountry_num_in_treatment_2016 = "";
    for country_num = 1:110
        ISO = ListOfISOs{country_num};
        thiscountry_num_in_treatment_2016 = jsonencode(num_in_treatment_2016_map(ISO),PrettyPrint=true);
        thiscountry_num_in_treatment_2016 = ISO + ":" + thiscountry_num_in_treatment_2016 + newline
        allcountry_num_in_treatment_2016 = strcat(allcountry_num_in_treatment_2016,thiscountry_num_in_treatment_2016);
    end
    writelines(allcountry_num_in_treatment_2016,"raw_params/num_in_treatment_2016.txt");




    allcountry_treatment_rates = "";
    for country_num = 1:110
        ISO = ListOfISOs{country_num};
        thiscountry_treatment_rates = jsonencode(treatment_rates_map(ISO),PrettyPrint=true);
        thiscountry_treatment_rates = ISO + ":" + thiscountry_treatment_rates + newline
        allcountry_treatment_rates = strcat(allcountry_treatment_rates,thiscountry_treatment_rates);
    end
    writelines(allcountry_treatment_rates,"raw_params/treatment_rates.txt");

end
