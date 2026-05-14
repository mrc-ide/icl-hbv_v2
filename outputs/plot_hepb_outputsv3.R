rm(list=ls())
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library("gridExtra")

## Output graphs saved here:
graph.dir <- "C:/Users/mpickles/OneDrive - Imperial College London/Dropbox_copy/Hepatits B/Presentations/29Apr2026/"


default.barplot.cols <- c("#f1595f","#599ad3","#f9a65a","#9e66ab","#cd7058","#d77fb3")
##comparison.cols <- c('Baseline'=default.barplot.cols[1], 'Intervention'=default.barplot.cols[2])
##barplot.cols.grey <- brewer.pal(4,"Greys")
plot.alpha.baseline <- 0.3

WHOregion.cols <- brewer.pal(7,"Spectral")

CountryClassification <- read.csv("C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/resources/CountryClassifier.csv",header=TRUE)

## GBD aggregated regions:
# Central Europe, eastern Europe, and central Asia
# High income
# Latin America and Caribbean
# North Africa and Middle East
# South Asia
# Southeast Asia, east Asia, and Oceania
# Sub-Saharan Africa


## WB classification:
##Low income
##Upper middle income
##Upper middle income
##High income


##ISO <- "ETH"
##ISO <- "GMB"

### Code from https://statisticsglobe.com/add-common-legend-to-combined-ggplot2-plots-in-r/
# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}


read_country_df <- function(ISO,years.to.get)
{
  df <- read.csv(paste0("C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/HPC/v1_runs/summary_outputs_",ISO,".csv"),sep=",", header=T, stringsAsFactors = FALSE)
  ##df <- df[df$Year %in% years.to.get,!(colnames(df) %in% c("country_ISO"))]
  df <- df[df$Year %in% years.to.get,]
  
  ## For now we remove the Wuenic 2019 scenario, in-facility BD, MAP and CPAD scenarios:
  df <- df[df$scenario %in% c(2,4,5,6,7,8),]
  
  df$scenario_name <- "X"
  ##df$scenario_name[df$scenario %in% 1] <- "WUENIC2019"
  df$scenario_name[df$scenario %in% 2] <- "Status quo" ## WUENIC 2025 + 2016 treat continued
  df$scenario_name[df$scenario %in% 4] <- "BD WHO target"
  df$scenario_name[df$scenario %in% 5] <- "HepB3 WHO target"
  df$scenario_name[df$scenario %in% 6] <- "Treat 40% coverage"
  df$scenario_name[df$scenario %in% 7] <- "Treat WHO target 80%"
  df$scenario_name[df$scenario %in% 8] <- "WHO targets"
  ##df$scenario_name[df$scenario %in% 9] <- "MAP"
  ##df$scenario_name[df$scenario %in% 10] <- "CPAD"
  
  
  ## Don't need the scenario variable any more:
  df <- subset(df, select = -scenario)
  
  ## Costs in US$ - placeholders
  unitcost_BD <- 0.38          
  unitcost_MAP <- 3
  unitcost_CPAD <- 5
  unitcost_HB3 <- 2.35
  cost_PAP <- 1
  cost_PAP_screen <- 1.23
  cost_treatelig_screen <- 1.23
  annualcost_treat <- 43.80
  annualcost_treat3TC <- 30
  
  ##df$cost <- unitcost_BD*df$NBirthDose + unitcost_MAP*df$NBD_MAP + unitcost_CPAD*df$NBD_CPAD + unitcost_HB3*df$N_InfantVacc + cost_PAP*(df$N_PAP_EAgHVL+df$N_PAP_EAgLVL+df$N_PAP_SAgHVL+df$N_PAP_SAgLVL) + cost_PAP_screen*df$N_screen_PAP + cost_treatelig_screen*df$N_starting_treatment + annualcost_treat*df$n_cases_TDFtreat + annualcost_treat3TC*(df$n_cases_3TCtreat+df$n_cases_3TCfailed)
  
  df$total_chronic <- df$n_cases_ImmTol + df$n_cases_ImmReact + df$n_cases_AsymptCarr + df$n_cases_Chr + df$n_cases_CompCirr + df$n_cases_DecompCirr + df$n_cases_HCC + df$n_cases_TDFtreat + df$n_cases_3TCtreat + df$n_cases_3TCfailed
  
  ## Proportion of eligible people who are on treatment:
  df$prop_treat <- (df$n_cases_TDFtreat)/(df$n_cases_TDFtreat + df$n_cases_ImmReact + df$n_cases_Chr + df$n_cases_CompCirr + df$n_cases_DecompCirr)
  df$Npop <- df$n_Susc + df$n_cases_ImmTol + df$n_cases_ImmReact + df$n_cases_AsymptCarr + df$n_cases_Chr + df$n_cases_CompCirr + df$n_cases_DecompCirr + df$n_cases_HCC + df$n_cases_Immune + df$n_cases_TDFtreat + df$n_cases_3TCtreat + df$n_cases_3TCfailed + df$n_cases_Acute
  df$prev.chronic <- df$total_chronic/df$Npop
  df$incidence_total <- df$incidence_neonatal + df$incidence_horiz0to4 + df$incidence_horiz5plus
  
  return(df)
}


make_stats_by_particle <- function(df)
{
  ## Boundaries for credible interval (95%)
  ll <- 0.025
  ul <- 0.975
  
  df_stats_by_run <- df %>% group_by(Year,scenario_name) %>% summarise(median_total_chronic=median(total_chronic),
                                                                       ll_total_chronic=percentile(total_chronic,ll),
                                                                       ul_total_chronic=percentile(total_chronic,ul),
                                                                       median_prop_treat=median(prop_treat),
                                                                       ll_prop_treat=percentile(prop_treat,ll),
                                                                       ul_prop_treat=percentile(prop_treat,ul),
                                                                       median_N=median(N),
                                                                       median_cost=median(cost),
                                                                       median_incidence_neonatal=median(incidence_neonatal),
                                                                       ll_incidence_neonatal=percentile(incidence_neonatal,ll),
                                                                       ul_incidence_neonatal=percentile(incidence_neonatal,ul),
                                                                       median_incidence_horiz0to4=median(incidence_horiz0to4),
                                                                       median_incidence_horiz0to4=median(incidence_horiz0to4),
                                                                       median_incidence_horiz5plus=median(incidence_horiz5plus),median_deaths=median(newdeath))
  
  return(df_stats_by_run)
}


## Regional stats only has some variables:
make_regional_stats_by_particle <- function(df)
{
  this.quantile <- 0.025
  df_stats_by_run <- df %>% group_by(Year,scenario_name) %>% summarise(mean_total_chronic=mean(total_chronic),
                                                                       sd_total_chronic=sd(total_chronic),
                                                                       median_total_chronic=median(total_chronic),
                                                                       ll_total_chronic=quantile(total_chronic,this.quantile),
                                                                       ul_total_chronic=quantile(total_chronic,(1-this.quantile)),
                                                                       mean_incidence_neonatal=mean(incidence_neonatal),
                                                                       sd_incidence_neonatal=sd(incidence_neonatal),
                                                                       median_incidence_neonatal=median(incidence_neonatal),
                                                                       ll_incidence_neonatal=quantile(incidence_neonatal,this.quantile),
                                                                       ul_incidence_neonatal=quantile(incidence_neonatal,(1-this.quantile)),
                                                                       mean_incidence_horiz0to4=mean(incidence_horiz0to4),
                                                                       sd_incidence_horiz0to4=sd(incidence_horiz0to4),
                                                                       median_incidence_horiz0to4=median(incidence_horiz0to4),
                                                                       ll_incidence_horiz0to4=quantile(incidence_horiz0to4,this.quantile),
                                                                       ul_incidence_horiz0to4=quantile(incidence_horiz0to4,(1-this.quantile)),
                                                                       mean_incidence_horiz5plus=mean(incidence_horiz5plus),
                                                                       sd_incidence_horiz5plus=sd(incidence_horiz5plus),
                                                                       median_incidence_horiz5plus=median(incidence_horiz5plus),
                                                                       ll_incidence_horiz5plus=quantile(incidence_horiz5plus,this.quantile),
                                                                       ul_incidence_horiz5plus=quantile(incidence_horiz5plus,(1-this.quantile)),
                                                                       mean_incidence_total=mean(incidence_total),
                                                                       sd_incidence_total=sd(incidence_total),
                                                                       median_incidence_total=median(incidence_total),
                                                                       ll_incidence_total=quantile(incidence_total,this.quantile),
                                                                       ul_incidence_total=quantile(incidence_total,(1-this.quantile)),
                                                                       mean_deaths=mean(newdeath),
                                                                       sd_deaths=sd(newdeath),
                                                                       median_deaths=median(newdeath),
                                                                       ll_deaths=quantile(newdeath,this.quantile),
                                                                       ul_deaths=quantile(newdeath,(1-this.quantile)),
                                                                       mean_prop_treat=mean(prop_treat),
                                                                       sd_prop_treat=sd(prop_treat),
                                                                       median_prop_treat=median(prop_treat),
                                                                       ll_prop_treat=quantile(prop_treat,this.quantile),
                                                                       ul_prop_treat=quantile(prop_treat,(1-this.quantile))
  )
  
  return(df_stats_by_run)
}


##make_WHO_region_df_by_run()
##df_WHO_Africa_byrun <- data.frame(matrix(ncol = 33, nrow = 0))
## Check if each df exists yet:
df.whoafrica.flag <- FALSE
df.whoamericas.flag <- FALSE
df.whoemed.flag <- FALSE
df.whoeurope.flag <- FALSE
df.whoseasia.flag <- FALSE
df.whowpacific.flag <- FALSE

## GBD aggregated regions:
df.gbdEurAsia.flag <- FALSE      # Central Europe, eastern Europe, and central Asia
df.gbdHighIncome.flag <- FALSE   # High income
df.gbdLatinAmCarib.flag <- FALSE # Latin America and Caribbean
df.gbdMENA.flag <- FALSE         # North Africa and Middle East
df.gbdSAsia.flag <- FALSE        # South Asia
df.gbdESEAsiaOce.flag <- FALSE   # Southeast Asia, east Asia, and Oceania
df.gbdSSA.flag <- FALSE          # Sub-Saharan Africa

## World Bank income classifications:
df.wbLow.flag <- FALSE        # Low income
df.wbLowerMid.flag <- FALSE   # Lower middle income
df.wbUpperMid.flag <- FALSE   # Upper middle income
df.wbHigh.flag <- FALSE       # High income



years.to.get <- seq(2000,2100)
nyears <- length(years.to.get)
for(i in 1:nrow(CountryClassification))
{
  country_ISO = CountryClassification$ISO[i]
  print(country_ISO)
  this.df <- read_country_df(country_ISO,years.to.get)
  
  if(i==1){
    ## Use the first country to set up variables:
    nruns <- max(this.df$run_number)
    nscenarios <- length(unique(this.df$scenario_name))
    
    ## We use this later on to add scenarios back to each regional DF:
    list.of.scenarios <- this.df$scenario_name
    list.of.years <- this.df$Year
    
  }
  
  shuffle.df <- data.frame(run_number=1:nruns, shuffled_run_number=sample(1:nruns))
  
  this.df <- merge(this.df,shuffle.df,by="run_number")
  ## Remove shuffle for now.
  this.df.shuffled <- this.df[order(this.df$scenario_name, this.df$shuffled_run_number, this.df$Year),]
  ##this.df.shuffled <- this.df[order(this.df$scenario_name, this.df$Year),]

  ## Shuffle the particles (stochastic runs):
  ##this.df.shuffled <- this.df[shuffle.this.country, colnames(this.df) %in% c("incidence_neonatal","incidence_horiz0to4","incidence_horiz5plus","newdeath","total_chronic")]
  
  this.df.shuffled <- this.df.shuffled[, colnames(this.df.shuffled) %in% c("incidence_neonatal","incidence_horiz0to4","incidence_horiz5plus","incidence_total","newdeath","total_chronic","prop_treat")]
  
  ## Group by WHO region:
  WHO_region = CountryClassification$WHO_region[i]
  if(WHO_region=="Africa"){
    if(df.whoafrica.flag==TRUE){
      df_WHO_Africa_byrun <- df_WHO_Africa_byrun + this.df.shuffled
    }else{
      df_WHO_Africa_byrun <- this.df.shuffled
      df.whoafrica.flag <- TRUE
    }
  }else if(WHO_region=="Americas"){
    if(df.whoamericas.flag==TRUE){
      df_WHO_Americas_byrun <- df_WHO_Americas_byrun + this.df.shuffled
    }else{
      df_WHO_Americas_byrun <- this.df.shuffled
      df.whoamericas.flag <- TRUE
    }
  }else if(WHO_region=="Eastern Mediterranean"){
    if(df.whoemed.flag==TRUE){
      df_WHO_EMed_byrun <- df_WHO_EMed_byrun + this.df.shuffled
    }else{
      df_WHO_EMed_byrun <- this.df.shuffled
      df.whoemed.flag <- TRUE
    }
  }else if(WHO_region=="Europe"){
    if(df.whoeurope.flag==TRUE){
      df_WHO_Europe_byrun <- df_WHO_Europe_byrun + this.df.shuffled
    }else{
      df_WHO_Europe_byrun <- this.df.shuffled
      df.whoeurope.flag <- TRUE
    }
  }else if(WHO_region=="South-East Asia"){
    if(df.whoseasia.flag==TRUE){
      df_WHO_SEAsia_byrun <- df_WHO_SEAsia_byrun + this.df.shuffled
    }else{
      df_WHO_SEAsia_byrun <- this.df.shuffled
      df.whoseasia.flag <- TRUE
    }
  }else if(WHO_region=="Western Pacific"){
    if(df.whowpacific.flag==TRUE){
      df_WHO_WPacific_byrun <- df_WHO_WPacific_byrun + this.df.shuffled
    }else{
      df_WHO_WPacific_byrun <- this.df.shuffled
      df.whowpacific.flag <- TRUE
    }
  }else
    print("Unknown region")
  
  ## Group by GBD aggregated region:
  ## Note - no high income countries currently modelled.
  GBD_region = CountryClassification$GBD_region_aggregated[i]
  if(GBD_region=="Sub-Saharan Africa"){
    if(df.gbdSSA.flag==TRUE){
      df_GBD_SSA_byrun <- df_GBD_SSA_byrun + this.df.shuffled
    }else{
      df_GBD_SSA_byrun <- this.df.shuffled
      df.gbdSSA.flag <- TRUE
    }
  }else if(GBD_region=="Southeast Asia, east Asia, and Oceania"){
    if(df.gbdESEAsiaOce.flag==TRUE){
      df_GBD_ESEAOc_byrun <- df_GBD_ESEAOc_byrun + this.df.shuffled
    }else{
      df_GBD_ESEAOc_byrun <- this.df.shuffled
      df.gbdESEAsiaOce.flag <- TRUE
    }
  }else if(GBD_region=="Latin America and Caribbean"){
    if(df.gbdLatinAmCarib.flag==TRUE){
      df_GBD_LAmCar_byrun <- df_GBD_LAmCar_byrun + this.df.shuffled
    }else{
      df_GBD_LAmCar_byrun <- this.df.shuffled
      df.gbdLatinAmCarib.flag <- TRUE
    }
  }else if(GBD_region=="Central Europe, eastern Europe, and central Asia"){
    if(df.gbdEurAsia.flag==TRUE){
      df_GBD_EurAs_byrun <- df_GBD_EurAs_byrun + this.df.shuffled
    }else{
      df_GBD_EurAs_byrun <- this.df.shuffled
      df.gbdEurAsia.flag <- TRUE
    }
  }else if(GBD_region=="North Africa and Middle East"){
    if(df.gbdMENA.flag==TRUE){
      df_GBD_MENA_byrun <- df_GBD_MENA_byrun + this.df.shuffled
    }else{
      df_GBD_MENA_byrun <- this.df.shuffled
      df.gbdMENA.flag <- TRUE
    }
  }else if(GBD_region=="South Asia"){
    if(df.gbdSAsia.flag==TRUE){
      df_GBD_SAsia_byrun <- df_GBD_SAsia_byrun + this.df.shuffled
    }else{
      df_GBD_SAsia_byrun <- this.df.shuffled
      df.gbdSAsia.flag <- TRUE
    }
  }else
    print("Unknown region")

  ### Group by World Bank income classification:
  ## Use income group, not WB regions:
  WB_Region = CountryClassification$WB_income_group[i]
  if(WB_Region=="High income"){
    if(df.wbHigh.flag==TRUE){
      df_WB_HighInc_byrun <- df_WB_HighInc_byrun + this.df.shuffled
    }else{
      df_WB_HighInc_byrun <- this.df.shuffled
      df.wbHigh.flag <- TRUE
    }
  }else if(WB_Region=="Upper middle income"){
    if(df.wbUpperMid.flag==TRUE){
      df_WB_UpperMid_byrun <- df_WB_UpperMid_byrun + this.df.shuffled
    }else{
      df_WB_UpperMid_byrun <- this.df.shuffled
      df.wbUpperMid.flag <- TRUE
    }
  }else if(WB_Region=="Lower middle income"){
    if(df.wbLowerMid.flag==TRUE){
      df_WB_LowerMid_byrun <- df_WB_LowerMid_byrun + this.df.shuffled
    }else{
      df_WB_LowerMid_byrun <- this.df.shuffled
      df.wbLowerMid.flag <- TRUE
    }
  }else if(WB_Region=="Low income"){
    if(df.wbLow.flag==TRUE){
      df_WB_Low_byrun <- df_WB_Low_byrun + this.df.shuffled
    }else{
      df_WB_Low_byrun <- this.df.shuffled
      df.wbLow.flag <- TRUE
    }
  }else
    print("Unknown region")
  
  
  

  
   
}
## Now make world total:
df_world_byrun <- df_WHO_Africa_byrun + df_WHO_Americas_byrun + df_WHO_EMed_byrun + df_WHO_Europe_byrun + df_WHO_SEAsia_byrun + df_WHO_WPacific_byrun


df_WB_world_byrun <- df_WB_HighInc_byrun + df_WB_UpperMid_byrun + df_WB_LowerMid_byrun + df_WB_Low_byrun

df_GBD_world_byrun <- df_GBD_SSA_byrun + df_GBD_ESEAOc_byrun + df_GBD_LAmCar_byrun + df_GBD_EurAs_byrun + df_GBD_MENA_byrun + df_GBD_SAsia_byrun

## Not quite identical, but within numerical accuracy:
max(abs(df_world_byrun - df_GBD_world_byrun))
## 8.940697e-08
max(abs(df_world_byrun - df_WB_world_byrun))
## 1.192093e-07

##df_world_byrun <- cbind(df_world_byrun, list.of.scenarios)

add_back_scenarios <- function(df, list.of.scenarios,list.of.years){
  df <- cbind(df, list.of.scenarios)
  df <- cbind(df, list.of.years)
  ## Add names if needed:
  if(!('scenario_name' %in% names(df)))
    names(df)[names(df) == 'list.of.scenarios'] <- 'scenario_name'
  else
    print("Error - duplicated scenario_name")
  if(!('Year' %in% names(df)))
    names(df)[names(df) == 'list.of.years'] <- 'Year'
  else
    print("Error - duplicated Year")
  return(df)
}



df_WHO_Africa_byrun <- add_back_scenarios(df_WHO_Africa_byrun, list.of.scenarios, list.of.years)
df_WHO_Americas_byrun <- add_back_scenarios(df_WHO_Americas_byrun, list.of.scenarios, list.of.years)
df_WHO_EMed_byrun <- add_back_scenarios(df_WHO_EMed_byrun, list.of.scenarios, list.of.years)
df_WHO_Europe_byrun <- add_back_scenarios(df_WHO_Europe_byrun, list.of.scenarios, list.of.years)
df_WHO_SEAsia_byrun <- add_back_scenarios(df_WHO_SEAsia_byrun, list.of.scenarios, list.of.years)
df_WHO_WPacific_byrun <- add_back_scenarios(df_WHO_WPacific_byrun, list.of.scenarios, list.of.years)
df_world_byrun <- add_back_scenarios(df_world_byrun, list.of.scenarios, list.of.years)


df_WB_HighInc_byrun <- add_back_scenarios(df_WB_HighInc_byrun, list.of.scenarios, list.of.years)
df_WB_UpperMid_byrun <- add_back_scenarios(df_WB_UpperMid_byrun, list.of.scenarios, list.of.years)
df_WB_LowerMid_byrun <- add_back_scenarios(df_WB_LowerMid_byrun, list.of.scenarios, list.of.years)
df_WB_Low_byrun <- add_back_scenarios(df_WB_Low_byrun, list.of.scenarios, list.of.years)
df_WB_world_byrun <- add_back_scenarios(df_WB_world_byrun, list.of.scenarios, list.of.years)


df_GBD_SSA_byrun <- add_back_scenarios(df_GBD_SSA_byrun, list.of.scenarios, list.of.years)
df_GBD_ESEAOc_byrun <- add_back_scenarios(df_GBD_ESEAOc_byrun, list.of.scenarios, list.of.years)
df_GBD_LAmCar_byrun <- add_back_scenarios(df_GBD_LAmCar_byrun, list.of.scenarios, list.of.years)
df_GBD_EurAs_byrun <- add_back_scenarios(df_GBD_EurAs_byrun, list.of.scenarios, list.of.years)
df_GBD_MENA_byrun <- add_back_scenarios(df_GBD_MENA_byrun, list.of.scenarios, list.of.years)
df_GBD_SAsia_byrun <- add_back_scenarios(df_GBD_SAsia_byrun, list.of.scenarios, list.of.years)
df_GBD_world_byrun <- add_back_scenarios(df_GBD_world_byrun, list.of.scenarios, list.of.years)

## We want to make stacked plots: (not needed?)
# cum_df_WHO_Africa_byrun <- df_WHO_Africa_byrun
# cum_df_WHO_Americas_byrun <- df_WHO_Africa_byrun + df_WHO_Americas_byrun
# cum_df_WHO_EMed_byrun <- df_WHO_Africa_byrun + df_WHO_Americas_byrun + df_WHO_EMed_byrun
# cum_df_WHO_Europe_byrun <- df_WHO_Africa_byrun + df_WHO_Americas_byrun + df_WHO_EMed_byrun + df_WHO_Europe_byrun
# cum_df_WHO_SEAsia_byrun <- df_WHO_Africa_byrun + df_WHO_Americas_byrun + df_WHO_EMed_byrun + df_WHO_Europe_byrun + df_WHO_SEAsia_byrun + df_WHO_WPacific_byrun
# cum_df_WHO_WPacific_byrun <- df_world_byrun
# 
# cum_df_WHO_Africa_byrun <- add_back_scenarios(cum_df_WHO_Africa_byrun, list.of.scenarios, list.of.years)
# cum_df_WHO_Americas_byrun <- add_back_scenarios(cum_df_WHO_Americas_byrun, list.of.scenarios, list.of.years)
# cum_df_WHO_EMed_byrun <- add_back_scenarios(cum_df_WHO_EMed_byrun, list.of.scenarios, list.of.years)
# cum_df_WHO_Europe_byrun <- add_back_scenarios(cum_df_WHO_Europe_byrun, list.of.scenarios, list.of.years)
# cum_df_WHO_SEAsia_byrun <- add_back_scenarios(cum_df_WHO_SEAsia_byrun, list.of.scenarios, list.of.years)
# cum_df_WHO_WPacific_byrun <- add_back_scenarios(cum_df_WHO_WPacific_byrun, list.of.scenarios, list.of.years)


df_WHO_Africa_stats <- make_regional_stats_by_particle(df_WHO_Africa_byrun)
df_WHO_Americas_stats <- make_regional_stats_by_particle(df_WHO_Americas_byrun)
df_WHO_EMed_stats <- make_regional_stats_by_particle(df_WHO_EMed_byrun)
df_WHO_Europe_stats <- make_regional_stats_by_particle(df_WHO_Europe_byrun)
df_WHO_SEAsia_stats <- make_regional_stats_by_particle(df_WHO_SEAsia_byrun)
df_WHO_WPacific_stats <- make_regional_stats_by_particle(df_WHO_WPacific_byrun)
df_world_stats <- make_regional_stats_by_particle(df_world_byrun)


df_WB_HighInc_stats <- make_regional_stats_by_particle(df_WB_HighInc_byrun)
df_WB_UpperMid_stats <- make_regional_stats_by_particle(df_WB_UpperMid_byrun)
df_WB_LowerMid_stats <- make_regional_stats_by_particle(df_WB_LowerMid_byrun)
df_WB_Low_stats <- make_regional_stats_by_particle(df_WB_Low_byrun)
df_WB_world_stats <- make_regional_stats_by_particle(df_WB_world_byrun)


df_GBD_SSA_stats <- make_regional_stats_by_particle(df_GBD_SSA_byrun)
df_GBD_ESEAOc_stats <- make_regional_stats_by_particle(df_GBD_ESEAOc_byrun)
df_GBD_LAmCar_stats <- make_regional_stats_by_particle(df_GBD_LAmCar_byrun)
df_GBD_EurAs_stats <- make_regional_stats_by_particle(df_GBD_EurAs_byrun)
df_GBD_MENA_stats <- make_regional_stats_by_particle(df_GBD_MENA_byrun)
df_GBD_SAsia_stats <- make_regional_stats_by_particle(df_GBD_SAsia_byrun)
df_GBD_world_stats <- make_regional_stats_by_particle(df_GBD_world_byrun)


# cum_df_WHO_Africa_stats <- make_regional_stats_by_particle(cum_df_WHO_Africa_byrun)
# cum_df_WHO_Americas_stats <- make_regional_stats_by_particle(cum_df_WHO_Americas_byrun)
# cum_df_WHO_EMed_stats <- make_regional_stats_by_particle(cum_df_WHO_EMed_byrun)
# cum_df_WHO_Europe_stats <- make_regional_stats_by_particle(cum_df_WHO_Europe_byrun)
# cum_df_WHO_SEAsia_stats <- make_regional_stats_by_particle(cum_df_WHO_SEAsia_byrun)
# cum_df_WHO_WPacific_stats <- make_regional_stats_by_particle(cum_df_WHO_WPacific_byrun)
df_WHO_Africa_stats$Region <- "Africa"
df_WHO_Americas_stats$Region <- "Americas"
df_WHO_EMed_stats$Region <- "Eastern Mediterranean"
df_WHO_Europe_stats$Region <- "Europe"
df_WHO_SEAsia_stats$Region <- "South-East Asia"
df_WHO_WPacific_stats$Region <- "Western Pacific"

df_WB_HighInc_stats$Region <- "High income"
df_WB_UpperMid_stats$Region <- "Upper middle income"
df_WB_LowerMid_stats$Region <- "Lower middle income"
df_WB_Low_stats$Region <- "Low income"

df_GBD_SSA_stats$Region <- "Sub-Saharan Africa"
df_GBD_ESEAOc_stats$Region <- "Southeast Asia, east Asia, and Oceania"
df_GBD_LAmCar_stats$Region <- "Latin America and Caribbean"
df_GBD_EurAs_stats$Region <- "Central Europe, eastern Europe, and central Asia"
df_GBD_MENA_stats$Region <- "North Africa and Middle East"
df_GBD_SAsia_stats$Region <- "South Asia"
##df_GBD_HighIncome_stats$Region <- "High income"


df_WHO_stats <- rbind(df_WHO_Africa_stats, df_WHO_Americas_stats)
df_WHO_stats <- rbind(df_WHO_stats, df_WHO_EMed_stats)
df_WHO_stats <- rbind(df_WHO_stats, df_WHO_Europe_stats)
df_WHO_stats <- rbind(df_WHO_stats, df_WHO_SEAsia_stats)
df_WHO_stats <- rbind(df_WHO_stats, df_WHO_WPacific_stats)
WHO_region_order <- c("Europe", "Eastern Mediterranean", "Americas", "South-East Asia", "Western Pacific", "Africa")
df_WHO_stats$Region <- factor(df_WHO_stats$Region, levels=WHO_region_order)

## WB:
df_WB_stats <- rbind(df_WB_HighInc_stats, df_WB_UpperMid_stats)
df_WB_stats <- rbind(df_WB_stats, df_WB_LowerMid_stats)
df_WB_stats <- rbind(df_WB_stats, df_WB_Low_stats)
WB_region_order <- c("High income","Upper middle income","Lower middle income","Low income")
df_WB_stats$Region <- factor(df_WB_stats$Region, levels=WB_region_order)

## GBD:
df_GBD_stats <- rbind(df_GBD_SSA_stats, df_GBD_ESEAOc_stats)
df_GBD_stats <- rbind(df_GBD_stats, df_GBD_LAmCar_stats)
df_GBD_stats <- rbind(df_GBD_stats, df_GBD_EurAs_stats)
df_GBD_stats <- rbind(df_GBD_stats, df_GBD_MENA_stats)
df_GBD_stats <- rbind(df_GBD_stats, df_GBD_SAsia_stats)
##df_GBD_stats <- rbind(df_GBD_stats, df_GBD_HighIncome_stats)
GBD_region_order <- c("Sub-Saharan Africa", "Southeast Asia, east Asia, and Oceania", "Latin America and Caribbean", "Central Europe, eastern Europe, and central Asia", "North Africa and Middle East", "South Asia")
##"High income"
df_GBD_stats$Region <- factor(df_GBD_stats$Region, levels=GBD_region_order)





## A simple plot: 
# this.scenario <- "Status quo"
# plot.df <- df_WHO_stats[df_WHO_stats$scenario_name %in% this.scenario,]
# p <- ggplot(plot.df, aes(x=Year, y=mean_total_chronic/1e6, fill=Region)) + 
#   geom_area() +
#   xlim(2020,2100) +
#   ylab("Mean total chronic (millions)")




## Simple area plot with overall CI:
this.scenario <- "Status quo"
plot.df <- df_WHO_stats[df_WHO_stats$scenario_name %in% this.scenario,]
total.ci.df <- df_world_stats[df_world_stats$scenario_name %in% this.scenario,]
total.ci.df$Region <- ""
ggplot(plot.df, aes(x=Year, y=mean_total_chronic, fill=Region)) + 
  geom_area() +
  geom_line(data=total.ci.df, aes(x=Year, y=mean_total_chronic +1.96*sd_total_chronic), color="black", size=0.5, linetype="dashed") +
  geom_line(data=total.ci.df, aes(x=Year, y=mean_total_chronic -1.96*sd_total_chronic), color="black", size=0.5, linetype="dashed") +
  xlim(2020,2100) +
  theme_bw() +
  labs(fill="WHO region") +
  theme(legend.position = "bottom") 




###################################################################################
## Testing - trying to make a plotting function
###################################################################################


## Generates a stacked area plot by region
## Plot "yvar" against Year (dividing yvar by y.factor e.g. to scale by millions)
make_area_plot <- function(plot.df, total.ci.df, y.lab, x.min=2020, x.max=2100, legend.lab)
{
  p <- ggplot(plot.df, aes(x=Year, y=y.var, fill=Region)) + 
      geom_area() +
      geom_line(data=total.ci.df, aes(x=Year, y=y.var +1.96*sd.var), color="black", size=0.5, linetype="dashed") +
      geom_line(data=total.ci.df, aes(x=Year, y=y.var -1.96*sd.var), color="black", size=0.5, linetype="dashed") +
      xlim(x.min,x.max) +
      ylim(0,NA) +
      theme_bw() +
      ylab(y.lab) + 
      labs(fill=legend.lab) +
      theme(legend.position = "bottom",axis.text=element_text(size=16),
            axis.title=element_text(size=20,face="bold"), legend.text = element_text(size=18), 
            legend.title=element_text(size=16, face="bold")) 
  return(p)  
}


## Plot chronic HIV area plot with CI by WHO region for status quo.
this.scenario <- "Status quo"
plot.df <- df_WHO_stats[df_WHO_stats$scenario_name %in% this.scenario,]
plot.df$y.var <- plot.df$mean_total_chronic/1e6
total.ci.df <- df_world_stats[df_world_stats$scenario_name %in% this.scenario,]
##total.ci.df$Region <- "World 95% CI"
total.ci.df$Region <- ""
total.ci.df$y.var <- total.ci.df$mean_total_chronic/1e6
total.ci.df$sd.var <- total.ci.df$sd_total_chronic/1e6
pdf(paste0(graph.dir,"Chronic_areaplot_byWHO.pdf"),width=16,height=12)
make_area_plot(plot.df, total.ci.df, y.lab="Mean total chronic (millions)", x.min=2020, x.max=2100, legend.lab="WHO region")
dev.off()


## Plot chronic HIV area plot with CI by GBD region for status quo.
this.scenario <- "Status quo"
plot.df <- df_GBD_stats[df_GBD_stats$scenario_name %in% this.scenario,]
plot.df$y.var <- plot.df$mean_total_chronic/1e6
total.ci.df <- df_world_stats[df_world_stats$scenario_name %in% this.scenario,]
##total.ci.df$Region <- "World 95% CI"
total.ci.df$Region <- ""
total.ci.df$y.var <- total.ci.df$mean_total_chronic/1e6
total.ci.df$sd.var <- total.ci.df$sd_total_chronic/1e6
pdf(paste0(graph.dir,"Chronic_areaplot_byGBD.pdf"),width=16,height=12)
make_area_plot(plot.df, total.ci.df, y.lab="Mean total chronic (millions)", x.min=2020, x.max=2100, legend.lab="GBD region")
dev.off()


## Plot chronic HIV area plot with CI by WB region for status quo.
this.scenario <- "Status quo"
plot.df <- df_WB_stats[df_WB_stats$scenario_name %in% this.scenario,]
plot.df$y.var <- plot.df$mean_total_chronic/1e6
total.ci.df <- df_world_stats[df_world_stats$scenario_name %in% this.scenario,]
##total.ci.df$Region <- "World 95% CI"
total.ci.df$Region <- ""
total.ci.df$y.var <- total.ci.df$mean_total_chronic/1e6
total.ci.df$sd.var <- total.ci.df$sd_total_chronic/1e6
pdf(paste0(graph.dir,"Chronic_areaplot_byWB.pdf"),width=16,height=12)
make_area_plot(plot.df, total.ci.df, y.lab="Mean total chronic (millions)", x.min=2020, x.max=2100, legend.lab="World bank region")
dev.off()

############################################################################
## Plot deaths

## Plot chronic HIV area plot with CI by WHO region for status quo.
this.scenario <- "Status quo"
plot.df <- df_WHO_stats[df_WHO_stats$scenario_name %in% this.scenario,]
plot.df$y.var <- plot.df$mean_deaths/1e3
total.ci.df <- df_world_stats[df_world_stats$scenario_name %in% this.scenario,]
##total.ci.df$Region <- "World 95% CI"
total.ci.df$Region <- ""
total.ci.df$y.var <- total.ci.df$mean_deaths/1e3
total.ci.df$sd.var <- total.ci.df$sd_deaths/1e3
pdf(paste0(graph.dir,"Deaths_areaplot_byWHO.pdf"),width=16,height=12)
p <- make_area_plot(plot.df, total.ci.df, y.lab="Mean annual deaths (thousands)", x.min=2020, x.max=2100, legend.lab="WHO region")
baseline.2015 <- total.ci.df$y.var[total.ci.df$Year %in% 2015]
elim.target <- baseline.2015*(1-0.65) ## 65% reduction in deaths
p <- p + geom_hline(yintercept = elim.target, color="black", size=1, linetype="dotted") +
  geom_vline(xintercept = 2030, color="blue", size=0.5, linetype="dashed")
p
dev.off()

## Plot chronic HIV area plot with CI by WB region for status quo.
this.scenario <- "Status quo"
plot.df <- df_WB_stats[df_WB_stats$scenario_name %in% this.scenario,]
plot.df$y.var <- plot.df$mean_deaths/1e3
total.ci.df <- df_world_stats[df_world_stats$scenario_name %in% this.scenario,]
##total.ci.df$Region <- "World 95% CI"
total.ci.df$Region <- ""
total.ci.df$y.var <- total.ci.df$mean_deaths/1e3
total.ci.df$sd.var <- total.ci.df$sd_deaths/1e3
pdf(paste0(graph.dir,"Deaths_areaplot_byWB.pdf"),width=16,height=12)
p <- make_area_plot(plot.df, total.ci.df, y.lab="Mean annual deaths (thousands)", x.min=2020, x.max=2100, legend.lab="WB region")
baseline.2015 <- total.ci.df$y.var[total.ci.df$Year %in% 2015]
elim.target <- baseline.2015*(1-0.65) ## 65% reduction in deaths
p <- p + geom_hline(yintercept = elim.target, color="black", size=1, linetype="dotted") +
  geom_vline(xintercept = 2030, color="blue", size=0.5, linetype="dashed")
p
dev.off()

## Plot chronic HIV area plot with CI by GBD region for status quo.
this.scenario <- "Status quo"
plot.df <- df_GBD_stats[df_GBD_stats$scenario_name %in% this.scenario,]
plot.df$y.var <- plot.df$mean_deaths/1e3
total.ci.df <- df_world_stats[df_world_stats$scenario_name %in% this.scenario,]
##total.ci.df$Region <- "World 95% CI"
total.ci.df$Region <- ""
total.ci.df$y.var <- total.ci.df$mean_deaths/1e3
total.ci.df$sd.var <- total.ci.df$sd_deaths/1e3
pdf(paste0(graph.dir,"Deaths_areaplot_byGBD.pdf"),width=16,height=12)
p <- make_area_plot(plot.df, total.ci.df, y.lab="Mean annual deaths (thousands)", x.min=2020, x.max=2100, legend.lab="GBD region")
baseline.2015 <- total.ci.df$y.var[total.ci.df$Year %in% 2015]
elim.target <- baseline.2015*(1-0.65) ## 65% reduction in deaths
p <- p + geom_hline(yintercept = elim.target, color="black", size=1, linetype="dotted") +
  geom_vline(xintercept = 2030, color="blue", size=0.5, linetype="dashed")
p
dev.off()



############################################################################
## Plot new cases:

scenarios.to.keep <- c("Status quo","BD WHO target","HepB3 WHO target","Treat 40% coverage","Treat WHO target 80%","WHO targets")
df_WHO_stats <- df_WHO_stats[df_WHO_stats$scenario_name %in% scenarios.to.keep, ]
df_WHO_stats$scenario_name <- factor(df_WHO_stats$scenario_name, levels=scenarios.to.keep)
## Plot chronic HIV area plot with CI by WHO region for status quo.
##this.scenario <- "Status quo"
make_line_plot_by_scenario <- function(plot.df, region.name, y.scale, y.lab, y.min, y.max, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=NA)
{

  plot.df <- plot.df[plot.df$Region %in% region.name,]
  
    p <- ggplot(plot.df, aes(x=Year, y=yvar, color=scenario_name)) + 
    geom_line() +
    theme_bw() +
    ylab(y.lab) + 
    xlim(x.min,x.max) +
    ylim(y.min,y.max) +
      ggtitle(region.name) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=16),
          plot.title = element_text(size = 20, face = "bold")) 

  if(show.legend==TRUE)
    p <- p + 
        labs(color=legend.lab) +
        theme(legend.position = "bottom", legend.text = element_text(size=18), 
          legend.title=element_text(size=16, face="bold")) +
        guides(color = guide_legend(override.aes = list(linewidth = 2))) ## wider line in legend
  else
    p <- p + theme(legend.position = "none")
  if(!is.na(reduction.target)){
    baseline.year <- 2015 ## WHO targets relative to 2015
    baseline.2015 <- plot.df$yvar[plot.df$Year %in% baseline.year & plot.df$scenario_name %in% "Status quo"]
    elim.target <- baseline.2015*(1-reduction.target)
    p <- p + geom_hline(yintercept = elim.target, color="black", size=0.5, linetype="dotted") +
      geom_vline(xintercept = 2030, color="blue", size=0.5, linetype="dashed")
  }
    
  return(p)  
}

##plot.df <- df_WHO_stats[df_WHO_stats$Region %in% "Africa",]
## We want to plot incidence in 0-4 y.o.s (in thousands):
df_WHO_stats$yvar <- df_WHO_stats$median_incidence_horiz0to4/1e3

p.dummy <- make_line_plot_by_scenario(df_WHO_stats, region.name="Africa", y.lab="Incident cases infants 0-4 (thousands)", y.min=0, y.max=250, x.min=2020, x.max=2100, legend.lab="Model scenario", show.legend=TRUE, reduction.target=0.95)
plot.legend <- extract_legend(p.dummy)

## Plot incident cases in Africa, compared to 95% reduction in target (compared to 2015):
p.africa <- make_line_plot_by_scenario(df_WHO_stats, region.name="Africa", y.lab="Incident cases infants 0-4 (thousands)", y.min=0, y.max=250, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)

#####
##plot.df <- df_WHO_stats[df_WHO_stats$Region %in% "Europe",]
p.europe <- make_line_plot_by_scenario(df_WHO_stats, region.name="Europe", y.lab="Incident cases infants 0-4 (thousands)", y.min=0, y.max=1.5, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.emed <- make_line_plot_by_scenario(df_WHO_stats, region.name="Eastern Mediterranean", y.lab="Incident cases infants 0-4 (thousands)", y.min=0, y.max=20, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.americas <- make_line_plot_by_scenario(df_WHO_stats, region.name="Americas", y.lab="Incident cases infants 0-4 (thousands)", y.min=0, y.max=1.5, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.seasia <- make_line_plot_by_scenario(df_WHO_stats, region.name="South-East Asia", y.lab="Incident cases infants 0-4 (thousands)", y.min=0, y.max=20, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.wpac <- make_line_plot_by_scenario(df_WHO_stats, region.name="Western Pacific", y.lab="Incident cases infants 0-4 (thousands)", y.min=0, y.max=25, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)


pdf(paste0(graph.dir,"Incidentcases0_4panel_byWHOregion.pdf"),width=16,height=12)
grid.arrange(arrangeGrob(p.africa, p.europe, p.emed, p.americas, p.seasia,p.wpac, ncol = 3)
             , plot.legend, nrow = 2, heights = c(12, 2))
dev.off()


##
df_WHO_stats$yvar <- df_WHO_stats$median_incidence_neonatal/1e3
p.africa <- make_line_plot_by_scenario(df_WHO_stats, region.name="Africa", y.lab="Incident neonatal cases (thousands)", y.min=0, y.max=700, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.europe <- make_line_plot_by_scenario(df_WHO_stats, region.name="Europe", y.lab="Incident neonatal cases (thousands)", y.min=0, y.max=5, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.emed <- make_line_plot_by_scenario(df_WHO_stats, region.name="Eastern Mediterranean", y.lab="Incident neonatal cases (thousands)", y.min=0, y.max=100, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.americas <- make_line_plot_by_scenario(df_WHO_stats, region.name="Americas", y.lab="Incident neonatal cases (thousands)", y.min=0, y.max=6, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.seasia <- make_line_plot_by_scenario(df_WHO_stats, region.name="South-East Asia", y.lab="Incident neonatal cases (thousands)", y.min=0, y.max=120, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.wpac <- make_line_plot_by_scenario(df_WHO_stats, region.name="Western Pacific", y.lab="Incident neonatal cases (thousands)", y.min=0, y.max=100, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)


pdf(paste0(graph.dir,"Incidentcases_neonatal_panel_byWHOregion.pdf"),width=16,height=12)
grid.arrange(arrangeGrob(p.africa, p.europe, p.emed, p.americas, p.seasia,p.wpac, ncol = 3)
             , plot.legend, nrow = 2, heights = c(12, 2))
dev.off()


p.africa <- make_line_plot_by_scenario(df_WHO_stats, region.name="Africa", y.lab="Incident cases 5+ (thousands)", y.min=0, y.max=700, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.europe <- make_line_plot_by_scenario(df_WHO_stats, region.name="Europe", y.lab="Incident cases 5+ (thousands)", y.min=0, y.max=6, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.emed <- make_line_plot_by_scenario(df_WHO_stats, region.name="Eastern Mediterranean", y.lab="Incident cases 5+ (thousands)", y.min=0, y.max=100, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.americas <- make_line_plot_by_scenario(df_WHO_stats, region.name="Americas", y.lab="Incident cases 5+ (thousands)", y.min=0, y.max=6, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.seasia <- make_line_plot_by_scenario(df_WHO_stats, region.name="South-East Asia", y.lab="Incident cases 5+ (thousands)", y.min=0, y.max=150, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.wpac <- make_line_plot_by_scenario(df_WHO_stats, region.name="Western Pacific", y.lab="Incident cases 5+ (thousands)", y.min=0, y.max=100, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)


pdf(paste0(graph.dir,"Incidentcases5plus_panel_byWHOregion.pdf"),width=16,height=12)
grid.arrange(arrangeGrob(p.africa, p.europe, p.emed, p.americas, p.seasia,p.wpac, ncol = 3)
             , plot.legend, nrow = 2, heights = c(12, 2))
dev.off()



ggplot(df_WHO_stats, aes(x=Year, y=mean_incidence_horiz0to4, color=scenario_name)) + 
  geom_line() +
  facet_wrap(~Region) +
  geom_line(data=total.ci.df, aes(x=Year, y=mean_total_chronic +1.96*sd_total_chronic), color="black", size=0.5, linetype="dashed") +
  geom_line(data=total.ci.df, aes(x=Year, y=mean_total_chronic -1.96*sd_total_chronic), color="black", size=0.5, linetype="dashed") +
  xlim(x.min,2100) +
  theme_bw() +
  labs(fill="WHO region") +
  theme(legend.position = "bottom") 

plot.df <- df_WHO_stats[df_WHO_stats$scenario_name %in% this.scenario,]
plot.df$y.var <- plot.df$mean_total_chronic/1e6
total.ci.df <- df_world_stats[df_world_stats$scenario_name %in% this.scenario,]
##total.ci.df$Region <- "World 95% CI"
total.ci.df$Region <- ""
total.ci.df$y.var <- total.ci.df$mean_total_chronic/1e6
total.ci.df$sd.var <- total.ci.df$sd_total_chronic/1e6
pdf(paste0(graph.dir,"Chronic_areaplot_byWHO.pdf"),width=16,height=12)
make_area_plot(plot.df, total.ci.df, y.lab="Mean total chronic (millions)", x.min=2020, x.max=2100, legend.lab="WHO region")
dev.off()




################################
## Treatment:

df_WHO_stats$yvar <- df_WHO_stats$median_prop_treat
p.africa <- make_line_plot_by_scenario(df_WHO_stats, region.name="Africa", y.lab="Treatment coverage (%)", y.min=0, y.max=50, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.europe <- make_line_plot_by_scenario(df_WHO_stats, region.name="Europe", y.lab="Treatment coverage (%)", y.min=0, y.max=50, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.emed <- make_line_plot_by_scenario(df_WHO_stats, region.name="Eastern Mediterranean", y.lab="Treatment coverage (%)", y.min=0, y.max=50, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.americas <- make_line_plot_by_scenario(df_WHO_stats, region.name="Americas", y.lab="Treatment coverage (%)", y.min=0, y.max=50, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.seasia <- make_line_plot_by_scenario(df_WHO_stats, region.name="South-East Asia", y.lab="Treatment coverage (%)", y.min=0, y.max=50, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
p.wpac <- make_line_plot_by_scenario(df_WHO_stats, region.name="Western Pacific", y.lab="Treatment coverage (%)", y.min=0, y.max=50, x.min=2020, x.max=2100, legend.lab=NA, show.legend=FALSE, reduction.target=0.95)
pdf(paste0(graph.dir,"Treatment_panel_byWHOregion.pdf"),width=16,height=12)
grid.arrange(arrangeGrob(p.africa, p.europe, p.emed, p.americas, p.seasia,p.wpac, ncol = 3)
             , plot.legend, nrow = 2, heights = c(12, 2))
dev.off()













##################################################################################


plot.df <- df_WHO_stats[df_WHO_stats$scenario_name %in% this.scenario,]
total.ci.df <- df_world_stats[df_world_stats$scenario_name %in% this.scenario,]


##total.ci.df$Region <- "World 95% CI"
total.ci.df$Region <- ""
ggplot(plot.df, aes(x=Year, y=mean_total_chronic, fill=Region)) + 
  geom_area() +
  geom_line(data=total.ci.df, aes(x=Year, y=mean_total_chronic +1.96*sd_total_chronic), color="black", size=0.5, linetype="dashed") +
  geom_line(data=total.ci.df, aes(x=Year, y=mean_total_chronic -1.96*sd_total_chronic), color="black", size=0.5, linetype="dashed") +
  xlim(2020,2100) +
  theme_bw() +
  labs(fill="WHO region") +
  theme(legend.position = "bottom") 

ggplot() +
  geom_line(data=, aes(x=Year, y=mean_total_chronic, color=WHOregion.cols[1]), linewidth=1, show.legend=FALSE) +
  geom_line(data=cum_df_WHO_Americas_stats[cum_df_WHO_Americas_stats$scenario_name %in% this.scenario,], aes(x=Year, y=mean_total_chronic, color=WHOregion.cols[2]), linewidth=1, show.legend=FALSE) +
  geom_line(data=cum_df_WHO_EMed_stats[cum_df_WHO_EMed_stats$scenario_name %in% this.scenario,], aes(x=Year, y=mean_total_chronic, color=WHOregion.cols[3]), linewidth=1, show.legend=FALSE) +
  geom_line(data=cum_df_WHO_SEAsia_stats[cum_df_WHO_SEAsia_stats$scenario_name %in% this.scenario,], aes(x=Year, y=mean_total_chronic, color=WHOregion.cols[4]), linewidth=1, show.legend=FALSE) +
  geom_line(data=cum_df_WHO_WPacific_stats[cum_df_WHO_WPacific_stats$scenario_name %in% this.scenario,], aes(x=Year, y=mean_total_chronic, color=WHOregion.cols[5]), linewidth=1, show.legend=FALSE) +
  geom_ribbon(data=plot.df, 
              aes(x=Year, ymin=100*ll,ymax=100*ul, fill="Baseline"), alpha=plot.alpha.baseline) +
  ##facet_wrap(~pop) +
  ylab(y.label) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  



df_long <- melt(df, id.vars=c("Year","run_number","scenario_name"))
##df_long <- df_long[,!(colnames(df_long) %in% "scenario")]

df_medians <- df %>% group_by(Year,scenario_name) %>% summarise(median_total_chronic=median(total_chronic),median_cost=median(cost),median_incidence_neonatal=median(incidence_neonatal),median_incidence_horiz0to4=median(incidence_horiz0to4),median_incidence_horiz5plus=median(incidence_horiz5plus),median_deaths=median(newdeath))

##df_medians <- df_medians[df_medians$Year %in% seq(1980,2100),]


df_medians_long <- melt(df_medians[,colnames(df_medians) %in% c("Year","scenario_name","median_incidence_neonatal","median_incidence_horiz0to4","median_incidence_horiz5plus")] , id.vars=c("Year","scenario_name"))


pdf(paste0("Areaplot_incidence_statusquo_",ISO,".pdf"))
ggplot(df_medians_long[df_medians_long$scenario_name %in% "Status quo",], aes(x=Year, y=value, fill=variable)) + 
    geom_area()
dev.off()

pdf(paste0("Areaplot_incidence_WHO_targets_",ISO,".pdf"))
ggplot(df_medians_long[df_medians_long$scenario_name %in% "WHO targets",], aes(x=Year, y=value, fill=variable)) + 
    geom_area()
dev.off()

                                        #
pdf(paste0("N_chronicHepB_",ISO,".pdf"))
ggplot(df_medians, aes(x=Year, y=median_total_chronic/1e6, group=scenario_name, color=scenario_name)) +
    theme_bw() +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) + 
    geom_line() +
    xlab("") +
    ylab("Median number of people living with chronic hepatitis B (millions)")
dev.off()




pdf(paste0("Median_newcaseschronic_",ISO,".pdf"))
ggplot(df_medians, aes(x=Year, y=(median_incidence_neonatal + median_incidence_horiz0to4+median_incidence_horiz5plus)/1e3, group=scenario_name, color=scenario_name)) +
    theme_bw() +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) + 
    geom_line() +
    xlab("") +
    ylab("Median new cases of chronic hepatitis B carriage (thousands)")
dev.off()


pdf(paste0("Median_deaths_",ISO,".pdf"))
ggplot(df_medians, aes(x=Year, y=median_deaths/1e3, group=scenario_name, color=scenario_name)) +
    theme_bw() +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) + 
    geom_line() +
    xlab("") +
    ylab("Median number of HBV deaths (thousands)")
dev.off()


pdf(paste0("MedianCosts_",ISO,".pdf"))
ggplot(df_medians, aes(x=Year, y=median_cost, group=scenario_name, color=scenario_name)) +
    theme_bw() +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) + 
    geom_line() +
    xlab("") +
    ylab("Median cost (US$)")
dev.off()



## plotit <- function(plot.df)
## {
##     p <- ggplot(plot.df, aes(x=Year, y=output.to.plot, group=scenario_name, color=scenario_name)) +
##         theme_bw() +
##         theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) + 
##         geom_line() + xlab("")
##     if(!is.na(y.lab))
##         p <- p + ylab(y.lab)
##     if(show.legend==FALSE){
##         p <- p + theme(legend.position="none")
##     }else{
##         p <- p + theme(legend.position="bottom")
##     }
##     if(!is.na(ymin))
##         p <- p + coord_cartesian(ylim=c(ymin,ymax))




ggplot() +
  geom_line(data=plot.df, aes(x=Year, y=100*Median, color="Baseline"), size=1, show.legend=FALSE) +
  geom_line(data=plot.df, aes(x=Year, y=100*ll, color="Baseline"), size=1, linetype="dashed", show.legend=FALSE) +
  geom_line(data=plot.df, aes(x=Year, y=100*ul, color="Baseline"), size=1, linetype="dashed", show.legend=FALSE) +
  geom_ribbon(data=plot.df, 
              aes(x=Year, ymin=100*ll,ymax=100*ul, fill="Baseline"), alpha=plot.alpha.baseline) +
  ##facet_wrap(~pop) +
  ylab(y.label) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  ## geom_errorbar(data=validation.data.df,aes(x=t,ymin=100*ll,ymax=100*ul),color="blue") +
  ## geom_point(data=validation.data.df,aes(x=t,y=100*prevalence),color="blue") +
  ## geom_errorbar(data=calibration.data.df,aes(x=t,ymin=100*ll,ymax=100*ul),color="black",width=2, size=1) +
  ## geom_point(data=calibration.data.df,aes(x=t,y=100*prevalence),color="black", size=2)+
  ## Add in the intervention:
  geom_line(data=plot.df.comparator, aes(x=Year, y=100*Median, color="Intervention"), size=1, show.legend=FALSE) +
  geom_line(data=plot.df.comparator, aes(x=Year, y=100*ll, color="Intervention"), size=1, linetype="dashed", show.legend=FALSE) +
  geom_line(data=plot.df.comparator, aes(x=Year, y=100*ul, color="Intervention"), size=1, linetype="dashed", show.legend=FALSE) +
  geom_ribbon(data=plot.df.comparator, 
              aes(x=Year, ymin=100*ll,ymax=100*ul, fill="Intervention"), alpha=plot.alpha.intervention)
