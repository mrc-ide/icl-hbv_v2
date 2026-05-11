
rm(list=ls())
GETDATA = 0

library(dplyr)
library(ggplot2)
##library(gifski)
##library(gganimate)
library(scales)
library(animation)
library(gganimate)
library(tidyr)

n.scenarios <- 10
n.runs <- 50

years.to.get <- seq(2000,2100)
##ISO <- "ETH"
ISO <- "GMB"

if(ISO=="GMB"){
  this.country = "Gambia"
}else if(ISO=="ETH"){
  this.country = "Ethiopia"
}

## Output graphs saved here:
graph.dir <- "C:/Users/mpickles/OneDrive - Imperial College London/Dropbox_copy/Hepatits B/Presentations/5May2026/"


header <- read.csv("C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/HPC/individual_runs/30/header.txt",header=T,sep=",")
header <- colnames(header)
##CountryClassification <- read.csv("C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/resources/CountryClassifier.csv",header=TRUE)

### Code from https://statisticsglobe.com/add-common-legend-to-combined-ggplot2-plots-in-r/
# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}



read_country_summary_df <- function(ISO,years.to.get)
{
  df <- read.csv(paste0("C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/HPC/summary_outputs_",ISO,".csv"),sep=",", header=T, stringsAsFactors = FALSE)
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



## Read a single run output for a single scenario in a single country:
read_country_df <- function(ISO,i_run,scenario_num,years.to.get,header)
{
  if(ISO=="GMB"){
    ISO_num=36
  }else if(ISO=="ETH"){
    ISO_num=30
  }else{
    print("Unknown ISO")
  }
  
  df <- read.csv(paste0("C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/HPC/individual_runs/",as.character(ISO_num),"/results_",ISO,"_scenario",as.character(scenario_num),"_default_run_",as.character(i_run),".csv"), sep=",", header=F, stringsAsFactors = FALSE)
  ##df <- df[df$Year %in% years.to.get,!(colnames(df) %in% c("country_ISO"))]
  colnames(df) <- header
  df <- df[df$Year %in% years.to.get,]
  nagegps <- 20
  ## Women, no treatment:
  i_Susc_F_Notreat =       2:(nagegps+1)
  i_ImmTol_F_Notreat =     (nagegps+2):(2*nagegps+1)
  i_ImmReact_F_Notreat =   (2*nagegps+2):(3*nagegps+1)
  i_AsymptCarr_F_Notreat = (3*nagegps+2):(4*nagegps+1)
  i_Chronic_F_Notreat =    (4*nagegps+2):(5*nagegps+1)
  i_CompCirr_F_Notreat =   (5*nagegps+2):(6*nagegps+1)
  i_DecompCirr_F_Notreat = (6*nagegps+2):(7*nagegps+1)
  i_HCC_F_Notreat =        (7*nagegps+2):(8*nagegps+1)
  i_Immune_F_Notreat =     (8*nagegps+2):(9*nagegps+1)
  i_TDFtreat_F_Notreat =   (9*nagegps+2):(10*nagegps+1)
  i_HBVdeath_F_Notreat =   (10*nagegps+2):(11*nagegps+1)
  i_ever3TC_F_Notreat =    (11*nagegps+2):(13*nagegps+1) ## 2 compartments
  i_NonSevAcute_F_Notreat = (13*nagegps+2):(14*nagegps+1)
  i_SevereAcute_F_Notreat = (14*nagegps+2):(15*nagegps+1)
  
  # header[i_Susc_F_Notreat]
  # header[i_ImmTol_F_Notreat]
  # header[i_ImmReact_F_Notreat]
  # header[i_AsymptCarr_F_Notreat]
  # header[i_Chronic_F_Notreat]
  # header[i_CompCirr_F_Notreat]
  # header[i_DecompCirr_F_Notreat]
  # header[i_HCC_F_Notreat]
  # header[i_Immune_F_Notreat]
  # header[i_TDFtreat_F_Notreat]
  # header[i_HBVdeath_F_Notreat]
  # header[i_ever3TC_F_Notreat] ## 2 groups.
  # header[i_NonSevAcute_F_Notreat]
  # header[i_SevereAcute_F_Notreat]

  ## Men
  i_Susc_M_Notreat =       (15*nagegps+2):(16*nagegps+1)
  i_ImmTol_M_Notreat =     (16*nagegps+2):(17*nagegps+1)
  i_ImmReact_M_Notreat =   (17*nagegps+2):(18*nagegps+1)
  i_AsymptCarr_M_Notreat = (18*nagegps+2):(19*nagegps+1)
  i_Chronic_M_Notreat =    (19*nagegps+2):(20*nagegps+1)
  i_CompCirr_M_Notreat =   (20*nagegps+2):(21*nagegps+1)
  i_DecompCirr_M_Notreat = (21*nagegps+2):(22*nagegps+1)
  i_HCC_M_Notreat =        (22*nagegps+2):(23*nagegps+1)
  i_Immune_M_Notreat =     (23*nagegps+2):(24*nagegps+1)
  i_TDFtreat_M_Notreat =   (24*nagegps+2):(25*nagegps+1)
  i_HBVdeath_M_Notreat =   (25*nagegps+2):(26*nagegps+1)
  i_ever3TC_M_Notreat =    (26*nagegps+2):(28*nagegps+1) ## 2 compartments
  i_NonSevAcute_M_Notreat = (28*nagegps+2):(29*nagegps+1)
  i_SevereAcute_M_Notreat = (29*nagegps+2):(30*nagegps+1)
  
  # header[i_Susc_M_Notreat]
  # header[i_ImmTol_M_Notreat]
  # header[i_ImmReact_M_Notreat]
  # header[i_AsymptCarr_M_Notreat]
  # header[i_Chronic_M_Notreat]
  # header[i_CompCirr_M_Notreat]
  # header[i_DecompCirr_M_Notreat]
  # header[i_HCC_M_Notreat]
  # header[i_Immune_M_Notreat]
  # header[i_TDFtreat_M_Notreat]
  # header[i_HBVdeath_M_Notreat]
  # header[i_ever3TC_M_Notreat]
  # header[i_NonSevAcute_M_Notreat]
  # header[i_SevereAcute_M_Notreat]
  #############################################################
  ### Now the Treat groups:
  ## Women, no treatment:
  i_Susc_F_Treat =       (30*nagegps+2):(31*nagegps+1)
  i_ImmTol_F_Treat =     (31*nagegps+2):(32*nagegps+1)
  i_ImmReact_F_Treat =   (32*nagegps+2):(33*nagegps+1)
  i_AsymptCarr_F_Treat = (33*nagegps+2):(34*nagegps+1)
  i_Chronic_F_Treat =    (34*nagegps+2):(35*nagegps+1)
  i_CompCirr_F_Treat =   (35*nagegps+2):(36*nagegps+1)
  i_DecompCirr_F_Treat = (36*nagegps+2):(37*nagegps+1)
  i_HCC_F_Treat =        (37*nagegps+2):(38*nagegps+1)
  i_Immune_F_Treat =     (38*nagegps+2):(39*nagegps+1)
  i_TDFtreat_F_Treat =   (39*nagegps+2):(40*nagegps+1)
  i_HBVdeath_F_Treat =   (40*nagegps+2):(41*nagegps+1)
  i_ever3TC_F_Treat =    (41*nagegps+2):(43*nagegps+1) ## 2 compartments
  i_NonSevAcute_F_Treat = (43*nagegps+2):(44*nagegps+1)
  i_SevereAcute_F_Treat = (44*nagegps+2):(45*nagegps+1)
  
  # header[i_Susc_F_Treat]
  # header[i_ImmTol_F_Treat]
  # header[i_ImmReact_F_Treat]
  # header[i_AsymptCarr_F_Treat]
  # header[i_Chronic_F_Treat]
  # header[i_CompCirr_F_Treat]
  # header[i_DecompCirr_F_Treat]
  # header[i_HCC_F_Treat]
  # header[i_Immune_F_Treat]
  # header[i_TDFtreat_F_Treat]
  # header[i_HBVdeath_F_Treat]
  # header[i_ever3TC_F_Treat] ## 2 groups.
  # header[i_NonSevAcute_F_Treat]
  # header[i_SevereAcute_F_Treat]
  
  ## Men
  i_Susc_M_Treat =       (45*nagegps+2):(46*nagegps+1)
  i_ImmTol_M_Treat =     (46*nagegps+2):(47*nagegps+1)
  i_ImmReact_M_Treat =   (47*nagegps+2):(48*nagegps+1)
  i_AsymptCarr_M_Treat = (48*nagegps+2):(49*nagegps+1)
  i_Chronic_M_Treat =    (49*nagegps+2):(50*nagegps+1)
  i_CompCirr_M_Treat =   (50*nagegps+2):(51*nagegps+1)
  i_DecompCirr_M_Treat = (51*nagegps+2):(52*nagegps+1)
  i_HCC_M_Treat =        (52*nagegps+2):(53*nagegps+1)
  i_Immune_M_Treat =     (53*nagegps+2):(54*nagegps+1)
  i_TDFtreat_M_Treat =   (54*nagegps+2):(55*nagegps+1)
  i_HBVdeath_M_Treat =   (55*nagegps+2):(56*nagegps+1)
  i_ever3TC_M_Treat =    (56*nagegps+2):(58*nagegps+1) ## 2 compartments
  i_NonSevAcute_M_Treat = (58*nagegps+2):(59*nagegps+1)
  i_SevereAcute_M_Treat = (59*nagegps+2):(60*nagegps+1)
  chronic_states = c(2,3,4,5,6,7,8,10,12,13)-1
  all_states_exc_death = c(1,2,3,4,5,6,7,8,9,10,12,13,14,15)-1
  i_U5_chronic = c(2+chronic_states*20,2+chronic_states*20+300,2+chronic_states*20+600,2+chronic_states*20+900)
  ##header[i_U5_chronic]

  i_U5_all = c(2+all_states_exc_death*20,2+all_states_exc_death*20+300,2+all_states_exc_death*20+600,2+all_states_exc_death*20+900)
  ##header[i_U5_all]
  
  ##header[i_Susc_M_Treat]
  # header[i_ImmTol_M_Treat]
  # header[i_ImmReact_M_Treat]
  # header[i_AsymptCarr_M_Treat]
  # header[i_Chronic_M_Treat]
  # header[i_CompCirr_M_Treat]
  # header[i_DecompCirr_M_Treat]
  # header[i_HCC_M_Treat]
  # header[i_Immune_M_Treat]
  # header[i_TDFtreat_M_Treat]
  # header[i_HBVdeath_M_Treat]
  # header[i_ever3TC_M_Treat]
  # header[i_NonSevAcute_M_Treat]
  # header[i_SevereAcute_M_Notreat]
  
  i_incidence_neonatal = which(header=="Incidence_neonatal")
  i_incidence_horiz0to4 = which(header=="IncidAge0_4")
  i_incidence_horiz5plus = seq(which(header=="IncidAge5_9"),which(header=="IncidAge95_99"))
  i_newdeath = seq(which(header=="DeathAge0_4"),which(header=="DeathAge95_99"))
  i_NBirthDose = which(header=="NBirthDose")
  i_NBD_MAP = which(header=="NBD_MAP")
  i_NBD_CPAD = which(header=="NBD_CPAD")
  i_N_InfantVacc = which(header=="N_InfantVacc")
  i_N_PAP_EAgHVL = which(header=="N_PAP_EAgHVL")
  i_N_PAP_EAgLVL = which(header=="N_PAP_EAgLVL")
  i_N_PAP_SAgHVL = which(header=="N_PAP_SAgHVL")
  i_N_PAP_SAgLVL = which(header=="N_PAP_SAgLVL")
  i_N_screen_PAP = which(header=="N_screen_PAP")
  i_N_starting_treatment = which(header=="N_starting_treatment")
  i_DALY = which(header=="DALYs")
  
  ### 
  df_processed = data.frame(Year=df$Year,i_run=i_run,scenario=scenario_num)
  df_processed$n_Susc = rowSums(df[,c(i_Susc_F_Notreat, i_Susc_M_Notreat, i_Susc_F_Treat, i_Susc_M_Treat)])
  df_processed$n_ImmTol = rowSums(df[,c(i_ImmTol_F_Notreat, i_ImmTol_M_Notreat, i_ImmTol_F_Treat, i_ImmTol_M_Treat)])
  df_processed$n_ImmReact = rowSums(df[,c(i_ImmReact_F_Notreat, i_ImmReact_M_Notreat, i_ImmReact_F_Treat, i_ImmReact_M_Treat)])
  df_processed$n_AsymptCarr = rowSums(df[,c(i_AsymptCarr_F_Notreat, i_AsymptCarr_M_Notreat, i_AsymptCarr_F_Treat, i_AsymptCarr_M_Treat)])
  df_processed$n_ChronCompartment = rowSums(df[,c(i_Chronic_F_Notreat, i_Chronic_M_Notreat, i_Chronic_F_Treat, i_Chronic_M_Treat)])
  df_processed$n_CompCirr = rowSums(df[,c(i_CompCirr_F_Notreat, i_CompCirr_M_Notreat, i_CompCirr_F_Treat, i_CompCirr_M_Treat)])
  df_processed$n_DecompCirr = rowSums(df[,c(i_DecompCirr_F_Notreat, i_DecompCirr_M_Notreat, i_DecompCirr_F_Treat, i_DecompCirr_M_Treat)])
  df_processed$n_HCC = rowSums(df[,c(i_HCC_F_Notreat, i_HCC_M_Notreat, i_HCC_F_Treat, i_HCC_M_Treat)])
  df_processed$n_TDFtreat = rowSums(df[,c(i_TDFtreat_F_Notreat, i_TDFtreat_M_Notreat, i_TDFtreat_F_Treat, i_TDFtreat_M_Treat)])
  df_processed$n_Immune = rowSums(df[,c(i_Immune_F_Notreat, i_Immune_M_Notreat, i_Immune_F_Treat, i_Immune_M_Treat)])
  ##df_processed$n_HBVdeath = rowSums(df[,c(i_HBVdeath_F_Notreat, i_HBVdeath_M_Notreat, i_HBVdeath_F_Treat, i_HBVdeath_M_Treat)])
  ##df_processed$n_ever3TC = rowSums(df[,c(i_ever3TC_F_Notreat, i_ever3TC_M_Notreat, i_ever3TC_F_Treat, i_ever3TC_M_Treat)])
  ##df_processed$n_NonSevAcute = rowSums(df[,c(i_NonSevAcute_F_Notreat, i_NonSevAcute_M_Notreat, i_NonSevAcute_F_Treat, i_NonSevAcute_M_Treat)])
  ##df_processed$n_SevereAcute = rowSums(df[,c(i_SevereAcute_F_Notreat, i_SevereAcute_M_Notreat, i_SevereAcute_F_Treat, i_SevereAcute_M_Treat)])
  df_processed$n_Acute = rowSums(df[,c(i_NonSevAcute_F_Notreat, i_NonSevAcute_M_Notreat, i_NonSevAcute_F_Treat, i_NonSevAcute_M_Treat,i_SevereAcute_F_Notreat, i_SevereAcute_M_Notreat, i_SevereAcute_F_Treat, i_SevereAcute_M_Treat)])
  
  df_processed$n_Chronic_sAgPos = df_processed$n_ImmTol + df_processed$n_ImmReact + df_processed$n_AsymptCarr + df_processed$n_ChronCompartment + df_processed$n_CompCirr + df_processed$n_DecompCirr + df_processed$n_HCC + df_processed$n_TDFtreat
  df_processed$Npop <- df_processed$n_Susc +df_processed$n_Chronic_sAgPos  + df_processed$n_Immune + df_processed$n_Acute
  
  df_processed$n_U5_chronic = rowSums(df[,i_U5_chronic])
  df_processed$n_U5_all = rowSums(df[,i_U5_all])

  ##df_processed$NoTreatStratification <- rowSums(df[,c(i_Susc_F_Notreat, i_Susc_M_Notreat, i_ImmTol_F_Notreat, i_ImmTol_M_Notreat, i_ImmReact_F_Notreat, i_ImmReact_M_Notreat, i_AsymptCarr_F_Notreat, i_AsymptCarr_M_Notreat, i_Chronic_F_Notreat, i_Chronic_M_Notreat, i_CompCirr_F_Notreat, i_CompCirr_M_Notreat, i_DecompCirr_F_Notreat, i_DecompCirr_M_Notreat, i_HCC_F_Notreat, i_HCC_M_Notreat, i_TDFtreat_F_Notreat, i_TDFtreat_M_Notreat, i_Immune_F_Notreat, i_Immune_M_Notreat, i_NonSevAcute_F_Notreat, i_NonSevAcute_M_Notreat,i_SevereAcute_F_Notreat, i_SevereAcute_M_Notreat)])
  ##df_processed$TreatStratification <- rowSums(df[,c(i_Susc_F_Treat, i_Susc_M_Treat, i_ImmTol_F_Treat, i_ImmTol_M_Treat, i_ImmReact_F_Treat, i_ImmReact_M_Treat, i_AsymptCarr_F_Treat, i_AsymptCarr_M_Treat, i_Chronic_F_Treat, i_Chronic_M_Treat, i_CompCirr_F_Treat, i_CompCirr_M_Treat, i_DecompCirr_F_Treat, i_DecompCirr_M_Treat, i_HCC_F_Treat, i_HCC_M_Treat, i_TDFtreat_F_Treat, i_TDFtreat_M_Treat, i_Immune_F_Treat, i_Immune_M_Treat, i_NonSevAcute_F_Treat, i_NonSevAcute_M_Treat,i_SevereAcute_F_Treat, i_SevereAcute_M_Treat)])
  
  ##max(df_processedf$TreatStratification)
  ## 0
  
  df_processed$prop_treat = df_processed$n_TDFtreat / (df_processed$n_TDFtreat + df_processed$n_ImmReact + df_processed$n_ChronCompartment + df_processed$n_CompCirr + df_processed$n_DecompCirr)
  df_processed$prop_chronic_sAgpos = df_processed$n_Chronic_sAgPos / df_processed$Npop
  
  
  
  # i_incidence_neonatal = which(header=="Incidence_neonatal")
  # i_incidence_horiz0to4 = which(header=="IncidAge0_4")
  # i_incidence_horiz5plus = seq(which(header=="IncidAge5_9"),which(header=="IncidAge95_99"))
  # i_newdeath = seq(which(header=="DeathAge0_4"),which(header=="DeathAge95_99"))
  # i_NBirthDose = which(header=="NBirthDose")
  # i_NBD_MAP = which(header=="NBD_MAP")
  # i_NBD_CPAD = which(header=="NBD_CPAD")
  # i_N_InfantVacc = which(header=="N_InfantVacc")
  # i_N_PAP_EAgHVL = which(header=="N_PAP_EAgHVL")
  # i_N_PAP_EAgLVL = which(header=="N_PAP_EAgLVL")
  # i_N_PAP_SAgHVL = which(header=="N_PAP_SAgHVL")
  # i_N_PAP_SAgLVL = which(header=="N_PAP_SAgLVL")
  # i_N_screen_PAP = which(header=="N_screen_PAP")
  # i_N_starting_treatment = which(header=="N_starting_treatment")
  # i_DALY = which(header=="DALYs")

  ##rowSums(df[,c(i_Susc_F_Notreat, i_Susc_M_Notreat, i_Susc_F_Treat, i_Susc_M_Treat)])
  df_processed$incidence_neonatal = df[,i_incidence_neonatal]
  df_processed$incidence_horiz0to4 = df[,i_incidence_horiz0to4]
  df_processed$incidence_horiz5plus = rowSums(df[,i_incidence_horiz5plus])
  df_processed$newdeath = rowSums(df[,i_newdeath])
  df_processed$NBirthDose = df[,i_NBirthDose]
  df_processed$NBD_MAP = df[,i_NBD_MAP]
  df_processed$NBD_CPAD = df[,i_NBD_CPAD]
  df_processed$N_InfantVacc = df[,i_N_InfantVacc]
  df_processed$N_PAP_EAgHVL = df[,i_N_PAP_EAgHVL]
  df_processed$N_PAP_EAgLVL = df[,i_N_PAP_EAgLVL]
  df_processed$N_PAP_SAgHVL = df[,i_N_PAP_SAgHVL]
  df_processed$N_PAP_SAgLVL = df[,i_N_PAP_SAgLVL]
  df_processed$N_screen_PAP = df[,i_N_screen_PAP]
  df_processed$N_starting_treatment = df[,i_N_starting_treatment]
  df_processed$N_DALYs = df[,i_DALY]
  
  df_processed$overall_incidence <- df_processed$incidence_neonatal + df_processed$incidence_horiz0to4 + df_processed$incidence_horiz5plus
  i_2015 <- which(df_processed$Year %in% 2015)
  df_processed$decrease_incidentcases = df_processed$overall_incidence / df_processed$overall_incidence[i_2015]
  df_processed$decrease_mortality = df_processed$newdeath / df_processed$newdeath[i_2015]
  
  
  return(df_processed)
  
}



# ## Regional stats only has some variables:
# make_regional_stats_by_particle <- function(df)
# {
#   this.quantile <- 0.025
#   df_stats_by_run <- df %>% group_by(Year,scenario_name) %>% summarise(mean_total_chronic=mean(total_chronic),
#                                                                        sd_total_chronic=sd(total_chronic),
#                                                                        median_total_chronic=median(total_chronic),
#                                                                        ll_total_chronic=quantile(total_chronic,this.quantile),
#                                                                        ul_total_chronic=quantile(total_chronic,(1-this.quantile)),
#                                                                        mean_incidence_neonatal=mean(incidence_neonatal),
#                                                                        sd_incidence_neonatal=sd(incidence_neonatal),
#                                                                        median_incidence_neonatal=median(incidence_neonatal),
#                                                                        ll_incidence_neonatal=quantile(incidence_neonatal,this.quantile),
#                                                                        ul_incidence_neonatal=quantile(incidence_neonatal,(1-this.quantile)),
#                                                                        mean_incidence_horiz0to4=mean(incidence_horiz0to4),
#                                                                        sd_incidence_horiz0to4=sd(incidence_horiz0to4),
#                                                                        median_incidence_horiz0to4=median(incidence_horiz0to4),
#                                                                        ll_incidence_horiz0to4=quantile(incidence_horiz0to4,this.quantile),
#                                                                        ul_incidence_horiz0to4=quantile(incidence_horiz0to4,(1-this.quantile)),
#                                                                        mean_incidence_horiz5plus=mean(incidence_horiz5plus),
#                                                                        sd_incidence_horiz5plus=sd(incidence_horiz5plus),
#                                                                        median_incidence_horiz5plus=median(incidence_horiz5plus),
#                                                                        ll_incidence_horiz5plus=quantile(incidence_horiz5plus,this.quantile),
#                                                                        ul_incidence_horiz5plus=quantile(incidence_horiz5plus,(1-this.quantile)),
#                                                                        mean_incidence_total=mean(incidence_total),
#                                                                        sd_incidence_total=sd(incidence_total),
#                                                                        median_incidence_total=median(incidence_total),
#                                                                        ll_incidence_total=quantile(incidence_total,this.quantile),
#                                                                        ul_incidence_total=quantile(incidence_total,(1-this.quantile)),
#                                                                        mean_deaths=mean(newdeath),
#                                                                        sd_deaths=sd(newdeath),
#                                                                        median_deaths=median(newdeath),
#                                                                        ll_deaths=quantile(newdeath,this.quantile),
#                                                                        ul_deaths=quantile(newdeath,(1-this.quantile)),
#                                                                        mean_prop_treat=mean(prop_treat),
#                                                                        sd_prop_treat=sd(prop_treat),
#                                                                        median_prop_treat=median(prop_treat),
#                                                                        ll_prop_treat=quantile(prop_treat,this.quantile),
#                                                                        ul_prop_treat=quantile(prop_treat,(1-this.quantile))
#   )
#   
#   return(df_stats_by_run)
# }


#############
## Quick check - plot treatemnt for a specific run:



##df <- read_country_summary_df(ISO,years.to.get)

##df <- data.frame()
if(GETDATA==1){
  for(i_run in 1:n.runs){
    for (scenario_num in 1:n.scenarios){
      print(i_run)
      if(i_run==1 & scenario_num==1)
      {
        df <- read_country_df(ISO,i_run,scenario_num,years.to.get,header)    
      }else{
        df <- rbind(df,read_country_df(ISO,i_run,scenario_num,years.to.get,header))
      }
      
    }
  }
  saveRDS(df, paste0(graph.dir,"data",ISO,".rds"))
}else{
  df <- readRDS(paste0(graph.dir,"data",ISO,".rds"))
}
df$i_run <- factor(df$i_run, levels=seq(1,n.runs))
scenarios.to.keep <- c(2,3,4,5,6,7,8)
df <- df[df$scenario %in% scenarios.to.keep,]
# iscenario_BASE2020 = 1; %% The 'default' scenario - WUENIC 2019 BD+HepB3, 2016 treatment, no new interventions.
# iscenario_BASE2025 = 2;     %% WUENIC 2025 BD+HepB3, 2016 treatment, no new interventions. Addresses - how have changes in BD+Hep B3 coverage impacted result?
# iscenario_INFACILITYBD = 3;  %% WUENIC 2025 HepB3, 2016 treatment, BD introduced in countries where it is not already present - coverage capped at in-facility birth coverage.
# iscenario_BDWHOtarget = 4;   %% BD reaches 90% coverage by T_INTERVENTION_END
# iscenario_HepB3WHOtarget = 5; %% HepB3 reach 90% coverage by T_INTERVENTION_END
# iscenario_Treat40percent = 6; %% Treatment reach 40% coverage. Treatment rate scales up over period T_INTERVENTION_START to T_INTERVENTION_END
# iscenario_TreatWHOtarget = 7; %% Treatment reach 80% coverage. Treatment rate scales up over period T_INTERVENTION_START to T_INTERVENTION_END
# iscenario_WHOtarget = 8;     %% BD+HepB3 reach 90% coverage by T_INTERVENTION_END, treatment reaches 80% by T_INTERVENTION_END
# iscenario_MAP = 9; %% WUENIC 2025 BD+HepB3, 2016 treatment, Microarray patch introduced in T_INTERVENTION_START (increase BD coverage, but lower efficacy).
# iscenario_CPAD = 10; %% WUENIC 2025 BD+HepB3, 2016 treatment, CPAD patch introduced (increase BD but lower eff and different cost to MAP).

df$scenario <- factor(df$scenario, levels=scenarios.to.keep)
df$scenario_name <- "X"
df$scenario_name[df$scenario %in% 2] <- "Status quo"
df$scenario_name[df$scenario %in% 3] <- "BD to in-facility"
df$scenario_name[df$scenario %in% 4] <- "BD 90%"
df$scenario_name[df$scenario %in% 5] <- "HepB3 90%"
df$scenario_name[df$scenario %in% 6] <- "Treat 40%"
df$scenario_name[df$scenario %in% 7] <- "Treat 80%"
df$scenario_name[df$scenario %in% 8] <- "All WHO targets"

df$scenario_name <- factor(df$scenario_name, levels=c("Status quo","BD to in-facility","BD 90%","HepB3 90%","Treat 40%","Treat 80%","All WHO targets"))

#######################################################################################
## Make summary stats for key quantities:
#######################################################################################
ll <- 0.025
ul <- 0.975
stats_df <- df %>% group_by(Year,scenario_name) %>% 
        summarise(median_total_chronic=median(n_Chronic_sAgPos),
            ll_total_chronic=quantile(n_Chronic_sAgPos,ll),
            ul_total_chronic=quantile(n_Chronic_sAgPos,ul),
            median_prop_treat=median(prop_treat),
            ll_prop_treat=quantile(prop_treat,ll),
            ul_prop_treat=quantile(prop_treat,ul),
            ##median_Npop=median(Npop),
            ##median_cost=median(cost),
            median_incidence_neonatal=median(incidence_neonatal),
            ll_incidence_neonatal=quantile(incidence_neonatal,ll),
            ul_incidence_neonatal=quantile(incidence_neonatal,ul),
            median_incidence_horiz0to4=median(incidence_horiz0to4),
            ll_incidence_horiz0to4=quantile(incidence_horiz0to4,ll),
            ul_incidence_horiz0to4=quantile(incidence_horiz0to4,ul),
            median_incidence_horiz5plus=median(incidence_horiz5plus),
            ll_incidence_horiz5plus=quantile(incidence_horiz5plus,ll),
            ul_incidence_horiz5plus=quantile(incidence_horiz5plus,ul),
            median_incidence_total=median(incidence_neonatal+incidence_horiz0to4+incidence_horiz5plus),
            ll_incidence_total=quantile(incidence_neonatal+incidence_horiz0to4+incidence_horiz5plus,ll),
            ul_incidence_total=quantile(incidence_neonatal+incidence_horiz0to4+incidence_horiz5plus,ul),
            median_deaths=median(newdeath),
            ll_deaths=quantile(newdeath,ll),
            ul_deaths=quantile(newdeath,ul),
            median_prev_u5=median(n_U5_chronic/n_U5_all),
            ll_prev_u5=quantile(n_U5_chronic/n_U5_all,ll),
            ul_prev_u5=quantile(n_U5_chronic/n_U5_all,ul)
            )


#######################################################################################
## Plot stuff:
#######################################################################################


#######################################################################################
## Plot stuff:
#######################################################################################

pdf(paste0(graph.dir,"TreatmentCoverage",ISO,".pdf"))
ggplot(stats_df, aes(x=Year, y=100*median_prop_treat, color=scenario_name)) + 
  geom_line() +
  theme_bw() +
  ylab("Median % treated")+
  ggtitle(this.country)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20, face = "bold")) 
dev.off()

pdf(paste0(graph.dir,"Incidence",ISO,".pdf"))
ggplot(stats_df, aes(x=Year, y=median_incidence_total, color=scenario_name)) + 
  geom_line() +
  theme_bw() +
  ylab("Median incident cases")+
  ggtitle(this.country)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20, face = "bold")) 
dev.off()

pdf(paste0(graph.dir,"PrevalenceU5",ISO,".pdf"))
ggplot(stats_df, aes(x=Year, y=100*median_prev_u5, color=scenario_name)) + 
  geom_line() +
  theme_bw() +
  ylab("Median prevalence in U5s (%)")+
  ggtitle(this.country)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20, face = "bold")) 
dev.off()

pdf(paste0(graph.dir,"Deaths",ISO,".pdf"))
ggplot(stats_df, aes(x=Year, y=median_deaths, color=scenario_name)) + 
  geom_line() +
  theme_bw() +
  ylab("Median annual deaths")+
  ggtitle(this.country)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20, face = "bold")) 
dev.off()

pdf(paste0(graph.dir,"NChronic",ISO,".pdf"))
ggplot(stats_df, aes(x=Year, y=median_total_chronic, color=scenario_name)) + 
  geom_line() +
  theme_bw() +
  ylab("Median number of chronic sAg+")+
  ggtitle(this.country)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20, face = "bold")) 
dev.off()



plot_all_runs <- function(df, i_scenario, x.lab, y.lab, main.lab, xmin=NA,xmax=NA, ymin=NA, ymax=NA,hline_intercept=NA)
{
  p <- ggplot(df[df$scenario %in% i_scenario,], aes(x=Year, y=yvar, group=i_run)) + 
    geom_line(color="gray") +
    theme_bw() +
    xlab(x.lab) +
    ylab(y.lab) +
    ggtitle(main.lab)
  if(!is.na(hline_intercept)){
    p <- p + geom_hline(yintercept = hline_intercept, color="black", size=1, linetype="dotted") 
  }
  if(!is.na(xmin)){
    p <- p + coord_cartesian(xlim=c(xmin,xmax))
  }
  if(!is.na(ymin)){
    p <- p + coord_cartesian(ylim=c(ymin,ymax))
  }
  return(p)
}

this.scenario <- 2
df$yvar <- 1-df$decrease_incidentcases
pdf(paste0(graph.dir,"Reduction_incidentcases",ISO,"_scenario",as.character(this.scenario),".pdf"))
plot_all_runs(df,this.scenario,x.lab="Year", y.lab="Decrease in incident cases", main.lab=this.country, xmin=2015,xmax=2100,ymin=0,ymax=1.1,hline_intercept=0.95)
dev.off()


this.scenario <- 8
df$yvar <- df$overall_incidence
pdf(paste0(graph.dir,"Incidentcases",ISO,"_scenario",as.character(this.scenario),".pdf"))
plot_all_runs(df,this.scenario,x.lab="Year", y.lab="Number of incident cases", main.lab=this.country, xmin=2015,xmax=2100,ymin=0,ymax=1.5e5,hline_intercept=NA)
dev.off()


this.scenario <- 2
df$yvar <- df$n_TDFtreat
pdf(paste0(graph.dir,"N_treat",ISO,"_scenario",as.character(this.scenario),".pdf"))
plot_all_runs(df,this.scenario,x.lab="Year", y.lab="Number on treatment", main.lab=this.country, xmin=2015,xmax=2100,ymin=NA,ymax=NA,hline_intercept=NA)
dev.off()

this.scenario <- 2
df$yvar <- df$n_ImmReact + df$n_ChronCompartment + df$n_CompCirr + df$n_DecompCirr
pdf(paste0(graph.dir,"N_elig_butnottreat",ISO,"_scenario",as.character(this.scenario),".pdf"))
plot_all_runs(df,this.scenario,x.lab="Year", y.lab="Number eligible, not on treatment", main.lab=this.country, xmin=2015,xmax=2100,ymin=NA,ymax=NA,hline_intercept=NA)
dev.off()



this.scenario <- 2
df$yvar <- 100*df$n_U5_chronic/df$n_U5_all
pdf(paste0(graph.dir,"PrevalenceU5",ISO,"_scenario",as.character(this.scenario),".pdf"))
plot_all_runs(df,this.scenario,x.lab="Year", y.lab="Prevalence in U5s", main.lab=this.country, xmin=2015,xmax=2100,ymin=0,ymax=2,hline_intercept=NA)
dev.off()


###

this.scenario <- 2
#for(t in seq(2030,2100,10))
#{
#  print(t)
t <- 2030


## Note - have to put a "print()" wrapper around ggplot otherwise it won't write.
make_scatter_prevU5_incdec <- function(df,t,this.scenario)
{
  pdf(paste0(graph.dir,"Scatter_PrevU5_IncDecrease",ISO,"_scenario",as.character(this.scenario),"scatter",as.character(t),".pdf"))
  print(ggplot(df[df$Year %in% t & df$scenario %in% this.scenario,], aes(x = 100*n_U5_chronic/n_U5_all, y = 100*(1-decrease_incidentcases))) + 
    geom_point(size=3) +
    xlab("Prevalence of chronic infection in U5s (%)") +
    ylab("Reduction in incident cases compared to 2015 (%)") +
    theme_bw() +
    xlim(0,0.5) +
    ylim(40,100) +
    ggtitle(paste0(this.country," t=",as.character(t))) +
    theme(legend.position = "bottom",axis.text=element_text(size=16),
          axis.title=element_text(size=20,face="bold"), legend.text = element_text(size=18), 
          legend.title=element_text(size=16, face="bold")) +
    geom_hline(yintercept = 95, color="blue", size=1, linetype="dashed") +
    geom_vline(xintercept = 0.1, color="blue", size=1, linetype="dashed") )
  dev.off()
}
#}


make_scatter_prevU5_incdec(df,2030,this.scenario=2)
make_scatter_prevU5_incdec(df,2040,this.scenario=2)
make_scatter_prevU5_incdec(df,2050,this.scenario=2)
make_scatter_prevU5_incdec(df,2060,this.scenario=2)
make_scatter_prevU5_incdec(df,2070,this.scenario=2)
make_scatter_prevU5_incdec(df,2080,this.scenario=2)
make_scatter_prevU5_incdec(df,2090,this.scenario=2)
make_scatter_prevU5_incdec(df,2100,this.scenario=2)

##
## Note - have to put a "print()" wrapper around ggplot otherwise it won't write.
make_scatter_prevU5_mortality <- function(df,t,this.scenario)
{
  pdf(paste0(graph.dir,"Scatter_PrevU5_Mortality",ISO,"_scenario",as.character(this.scenario),"scatter",as.character(t),".pdf"))
  print(ggplot(df[df$Year %in% t & df$scenario %in% this.scenario,], aes(x = 100*n_U5_chronic/n_U5_all, y = 100*decrease_mortality)) + 
          geom_point(size=3) +
          xlab("Prevalence of chronic infection in U5s (%)") +
          ylab("Change in deaths compared to 2015 (%)") +
          theme_bw() +
          xlim(0,0.5) +
          ylim(0,160) +
          ggtitle(paste0(this.country," t=",as.character(t))) +
          theme(legend.position = "bottom",axis.text=element_text(size=16),
                axis.title=element_text(size=20,face="bold"), legend.text = element_text(size=18), 
                legend.title=element_text(size=16, face="bold")) +
          geom_hline(yintercept = 35, color="blue", size=1, linetype="dashed") +
          geom_vline(xintercept = 0.1, color="blue", size=1, linetype="dashed") )
  dev.off()
}
#}


make_scatter_prevU5_mortality(df,2030,this.scenario=2)
make_scatter_prevU5_mortality(df,2040,this.scenario=2)
make_scatter_prevU5_mortality(df,2050,this.scenario=2)
make_scatter_prevU5_mortality(df,2060,this.scenario=2)
make_scatter_prevU5_mortality(df,2070,this.scenario=2)
make_scatter_prevU5_mortality(df,2080,this.scenario=2)
make_scatter_prevU5_mortality(df,2090,this.scenario=2)
make_scatter_prevU5_mortality(df,2100,this.scenario=2)





## Note - have to put a "print()" wrapper around ggplot otherwise it won't write.
make_scatter_prevU5_incdec_scenarios <- function(df,t,scenario.list)
{
  plot.df <- df[df$Year %in% t & df$scenario %in% scenario.list,]
  plot.df$scenario_name <- droplevels(plot.df$scenario_name)
  pdf(paste0(graph.dir,"Scenario_Scatter_PrevU5_IncDecrease",ISO,"_",as.character(t),".pdf"))
  print(ggplot(plot.df, aes(x = 100*n_U5_chronic/n_U5_all, y = 100*(1-decrease_incidentcases), color=scenario_name)) + 
          geom_point(size=3) +
          xlab("Prevalence of chronic infection in U5s (%)") +
          ylab("Reduction in incident cases compared to 2015 (%)") +
          theme_bw() +
          xlim(0,0.5) +
          ylim(40,100) +
          ggtitle(paste0(this.country," t=",as.character(t))) +
          theme(legend.position = "bottom",axis.text=element_text(size=16),
                axis.title=element_text(size=20,face="bold"), legend.text = element_text(size=18), 
                legend.title=element_text(size=16, face="bold")) +
          geom_hline(yintercept = 95, color="blue", size=1, linetype="dashed") +
          geom_vline(xintercept = 0.1, color="blue", size=1, linetype="dashed") )
  dev.off()
}
#}
this.scenario.list <- c(2,4,8)
make_scatter_prevU5_incdec_scenarios(df,2030,this.scenario.list)
make_scatter_prevU5_incdec_scenarios(df,2040,this.scenario.list)
make_scatter_prevU5_incdec_scenarios(df,2050,this.scenario.list)
make_scatter_prevU5_incdec_scenarios(df,2060,this.scenario.list)






#######################################################################
## Attempting to animate the above:
#######################################################################


p <- ggplot(df[df$Year>2025 & df$scenario %in% this.scenario,], aes(x = 100*n_U5_chronic/n_U5_all, y = 100*decrease_mortality)) + 
       geom_point(size=3) +
       xlab("Prevalence of chronic infection in U5s (%)") +
       ylab("Change in deaths compared to 2015 (%)") +
       theme_bw() +
       xlim(0,0.5) +
       ylim(0,160) +
      labs(title = paste0(this.country,"Year: {closest_state}")) +
       theme(legend.position = "bottom",axis.text=element_text(size=16),
             axis.title=element_text(size=20,face="bold"), legend.text = element_text(size=18), 
             legend.title=element_text(size=16, face="bold")) +
       geom_hline(yintercept = 35, color="blue", size=1, linetype="dashed") +
       geom_vline(xintercept = 0.1, color="blue", size=1, linetype="dashed") 


p_animated <- p + transition_states(Year,
                                    transition_length = 2,
                                    state_length = 2) +
  enter_fade() +
  exit_fade() +
  ease_aes("linear")

anim_save(paste0(graph.dir,"scatterplot_animation_GMB_U5prev_deaths.gif"),p_animated)








#######################################################################
#######################################################################
## Unused code:
#######################################################################
#######################################################################


plot.df <- df[df$run_number %in% 1,]
p <- ggplot(plot.df, aes(x=Year, y=prop_treat, color=scenario_name)) + 
  geom_line() +
  theme_bw() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20, face = "bold")) 



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

  
  

  
    
######################################################















setwd("C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/scripts/")
current.dir <- "../HPC/"

this.country <- "GMB"


get_header <- function(filedir)
{
    header <- t(read.csv(paste0(filedir,"header.txt"),header=F,stringsAsFactors=FALSE))
    return(header)
}
    

## loadfilesv1 <- function(n.files,filedir, header)
## {
##     filelist <- list()    
##     for(i in 1:n.files)
##     {
##         infile <- read.csv(paste0(filedir,"results_THA_scenario1_default_run_",as.character(i),".csv"),header=FALSE)
##         infile <- infile[1:(nrow(infile)-1),]
##         colnames(infile) <- c("Year",header)
##         filelist[[i]] <- infile
##     }
##     return(filelist)
## }


## Read in the csv, and format as needed:
loadfiles <- function(i,country_ISO, filedir, header, scenario)
{
    output <- read.csv(paste0(filedir,"results_",country_ISO,"_scenario",as.character(scenario),"_default_run_",as.character(i),".csv"),header=FALSE)
    ## Last column is -1
    output <- output[1:(nrow(output)-1),]
    ## Add column names:
    colnames(output) <- c("Year",header)
    output$scenario <- scenario
    output$runnumber <- i
    return(output)
}



## Function to get mean of n files:
make_average <- function(country_ISO, n.files,filedir, header, scenario)
{

    output <- loadfiles(1,country_ISO, filedir, header, scenario)
    
    if(n.files>1)
    {
        for(i in 2:n.files)
        {
            infile <- loadfiles(i,country_ISO, filedir, header, scenario)
            output <- output + infile
        }
        output <- output/n.files
    }
    
    return(output)
}


load_files_as_list <- function(country_ISO, n.files,filedir, header, scenario)
{
  
  list_of_files <- vector("list", n.files)
  
  

  for(i in 1:n.files)
    {
      infile <- loadfiles(i,country_ISO, filedir, header, scenario)
      list_of_files[[i]] <- infile
    }
  
  
  return(list_of_files)
}


load_files_as_df <- function(country_ISO, n.files,filedir, header, scenario)
{
  
  df <- data.frame()
  
  
  
  for(i in 1:n.files)
  {
    infile <- loadfiles(i,country_ISO, filedir, header, scenario)
    df <- rbind(df,infile)
  }
  
  
  return(df)
}



header <- get_header(current.dir)
if(is.na(header[length(header)]))
{
    header <- header[1:(length(header)-1)]
}


##currentfiles <- loadfiles(n.runs,country_ISO, current.dir, header)

##country_files <- load_files_as_list("THA",3,current.dir, header)


country_files <- data.frame()
for(scenario in 1:n.scenarios)
{
  country_files <- rbind(country_files,load_files_as_df(this.country,n.runs,current.dir, header, scenario))
}


##t <- currentfiles[[1]]$Year


## For a specific age and sex group (summed over Treat/No treat) get the indices for a certain natural history class (e.g. alive)
get_i_nathist_by_age_sex <- function(age_gp, sex, nat_hist_state){

    ## Overall population size:
    if(nat_hist_state=="Alive")
    {
        i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D[0-9]Treat"), header), grep(paste0("Age",age_gp,sex,"_D1[02345]Treat"), header), grep(paste0("Age",age_gp,sex,"_D[0-9]NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D1[02345]NoTreat"), header)))
    ## Cumulative deaths
    }else if(nat_hist_state=="CumulativeHBVDeaths")
    {
        i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D11Treat"), header), grep(paste0("Age",age_gp,sex,"_D11NoTreat"), header)))
    }else if(nat_hist_state=="Acute")
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D1[4-5]NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D1[4-5]NoTreat"), header)))
    }else if(nat_hist_state=="HBsAgpos_eAgneg_chronic") ## sAg+, eAg-, and not in acute phase.
    {
        i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D[2-8]Treat"), header), grep(paste0("Age",age_gp,sex,"_D10Treat"), header), grep(paste0("Age",age_gp,sex,"_D1[2-3]Treat"), header), grep(paste0("Age",age_gp,sex,"_D[2-8]NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D10NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D1[2-3]NoTreat"), header)))
    }else if(nat_hist_state=="Susc") ## Compartment 1 - 
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D1Treat"), header), grep(paste0("Age",age_gp,sex,"_D1NoTreat"), header)))
    }else if(nat_hist_state=="ImmuneTolerant") ## Compartment 2 - 
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D2Treat"), header), grep(paste0("Age",age_gp,sex,"_D2NoTreat"), header)))
    }else if(nat_hist_state=="ImmuneReactive") ## Compartment 3 - 
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D3Treat"), header), grep(paste0("Age",age_gp,sex,"_D3NoTreat"), header)))
    }else if(nat_hist_state=="AsymptCarrier") ## Compartment 4 - 
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D4Treat"), header), grep(paste0("Age",age_gp,sex,"_D4NoTreat"), header)))
    }else if(nat_hist_state=="ChronicHepB") ## Compartment 5 - 
    {
        i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D5Treat"), header), grep(paste0("Age",age_gp,sex,"_D5NoTreat"), header)))
    }else if(nat_hist_state=="CompCirr") ## Compartment 6
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D6Treat"), header), grep(paste0("Age",age_gp,sex,"_D6NoTreat"), header)))
    }else if(nat_hist_state=="DecompCirr") ## Compartment 7
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D7Treat"), header), grep(paste0("Age",age_gp,sex,"_D7NoTreat"), header)))
    }else if(nat_hist_state=="HCC") ## Compartment 8
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D8Treat"), header), grep(paste0("Age",age_gp,sex,"_D8NoTreat"), header)))
    }else if(nat_hist_state=="Immune") ## Compartment 9
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D9Treat"), header), grep(paste0("Age",age_gp,sex,"_D9NoTreat"), header)))
    }else if(nat_hist_state=="OnTreatmentTDF") ## Currently getting TDF
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D10Treat"), header), grep(paste0("Age",age_gp,sex,"_D10NoTreat"), header)))
    }else if(nat_hist_state=="OnTreatment3TC") ## Currently getting 3TC
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D12Treat"), header), grep(paste0("Age",age_gp,sex,"_D12NoTreat"), header)))
    }else if(nat_hist_state=="FailedTreatment3TC") ## Failed 3TC
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D13Treat"), header), grep(paste0("Age",age_gp,sex,"_D13NoTreat"), header)))
    }else if(nat_hist_state=="NonSevereAcute") 
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D14Treat"), header), grep(paste0("Age",age_gp,sex,"_D14NoTreat"), header)))
    }else if(nat_hist_state=="SevereAcute") 
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D15Treat"), header), grep(paste0("Age",age_gp,sex,"_D15NoTreat"), header)))
      
      ####################################################################################################
      ## Now other strata of interest:
      ####################################################################################################
    }else if(nat_hist_state=="TreatStratum") ## this person can be reached by treatment progs
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D[0-9]Treat"), header), grep(paste0("Age",age_gp,sex,"_D1[02345]Treat"), header)))
    }else if(nat_hist_state=="NoTreatStratum") ## this person cannot be reached by treatment progs
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D[0-9]NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D1[02345]NoTreat"), header)))
    }else if(nat_hist_state=="OnTreat_NoTreatStratum") ## this person cannot be reached by treatment progs but is on treatment - should be zero.
    {
      i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D10Treat"), header), grep(paste0("Age",age_gp,sex,"_D10NoTreat"), header)))
    }

    ##[2:8 10 12:15]
    return(i)
}


## testing:
# states <- c("Alive", "CumulativeHBVDeaths", "Acute","HBsAgpos_eAgneg_chronic","Susc","ImmuneTolerant","ImmuneReactive","AsymptCarrier","ChronicHepB","CompCirr","DecompCirr","HCC","Immune","OnTreatmentTDF","OnTreatment3TC","FailedTreatment3TC","NonSevereAcute","SevereAcute","TreatStratum","NoTreatStratum","OnTreat_NoTreatStratum")
# i=5
# print(states[i])
# x <- get_i_nathist_by_age_sex("25_29", "M", states[i])
# colnames(country_files)[x]

## Make age group labs:
min.age <- 0
max.age <- 99
age.gp.width <- 5
n.age.groups <- (max.age+1-min.age)/age.gp.width
age_group_labs <- rep("",n.age.groups)
for (i in 1:n.age.groups)
{
    i.age <- (i-1)*age.gp.width + min.age
    age_group_labs[i] <- paste(as.character(i.age),as.character(i.age+age.gp.width-1),sep="_")
}

i_all_alive_F <- c()
i_all_alive_M <- c()
i_HBVdeath_cumulative <- c()
i_chronic_sAgPos <- c()
i_Susc <- c()
i_ImmuneTolerant <- c()
i_ImmuneReactive <- c()
i_AsymptCarrier <- c()
i_ChronicHepB <- c()
i_CompCirr <- c()
i_DecompCirr <- c()
i_HCC <- c()
i_Immune <- c()
i_OnTreatmentTDF <- c()
i_OnTreatment3TC <- c()
i_FailedTreatment3TC <- c()
i_NonSevereAcute <- c()
i_SevereAcute <- c()

## for(i in length(age_group_labs))
## {
##     i_chronic_sAgPos_byage <- append(i_chronic_sAgPos_byage,c())
##     i_HCC_byage <- append(i_HCC_byage,c())
## }

##states <- c("Alive", "CumulativeHBVDeaths", "Acute","HBsAgpos_eAgneg_chronic","Susc","ImmuneTolerant","ImmuneReactive","AsymptCarrier","ChronicHepB","CompCirr","DecompCirr","HCC","Immune","OnTreatmentTDF","OnTreatment3TC","FailedTreatment3TC","NonSevereAcute","SevereAcute","TreatStratum","NoTreatStratum","OnTreat_NoTreatStratum")
for(i in 1:length(age_group_labs))
{
    i_all_alive_F <- c(i_all_alive_F,get_i_nathist_by_age_sex(age_group_labs[i], "F", "Alive"))
    i_all_alive_M <- c(i_all_alive_M,get_i_nathist_by_age_sex(age_group_labs[i], "M", "Alive"))
    i_HBVdeath_cumulative <- c(i_HBVdeath_cumulative, get_i_nathist_by_age_sex(age_group_labs[i], "M", "CumulativeHBVDeaths"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "CumulativeHBVDeaths"))
    i_Susc <- c(i_Susc, get_i_nathist_by_age_sex(age_group_labs[i], "M", "Susc"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "Susc"))
    i_ImmuneTolerant <- c(i_ImmuneTolerant, get_i_nathist_by_age_sex(age_group_labs[i], "M", "ImmuneTolerant"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "ImmuneTolerant"))
    i_ImmuneReactive <- c(i_ImmuneReactive, get_i_nathist_by_age_sex(age_group_labs[i], "M", "ImmuneReactive"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "ImmuneReactive"))
    i_AsymptCarrier <- c(i_AsymptCarrier, get_i_nathist_by_age_sex(age_group_labs[i], "M", "AsymptCarrier"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "AsymptCarrier"))
    i_ChronicHepB <- c(i_ChronicHepB, get_i_nathist_by_age_sex(age_group_labs[i], "M", "ChronicHepB"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "ChronicHepB"))
    i_CompCirr <- c(i_CompCirr, get_i_nathist_by_age_sex(age_group_labs[i], "M", "CompCirr"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "CompCirr"))
    i_DecompCirr <- c(i_DecompCirr, get_i_nathist_by_age_sex(age_group_labs[i], "M", "DecompCirr"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "DecompCirr"))
    i_HCC <- c(i_HCC, get_i_nathist_by_age_sex(age_group_labs[i], "F", "HCC"),get_i_nathist_by_age_sex(age_group_labs[i], "M", "HCC"))
    i_Immune <- c(i_Immune, get_i_nathist_by_age_sex(age_group_labs[i], "M", "Immune"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "Immune"))
    i_OnTreatmentTDF <- c(i_OnTreatmentTDF, get_i_nathist_by_age_sex(age_group_labs[i], "M", "OnTreatmentTDF"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "OnTreatmentTDF"))
    i_OnTreatment3TC <- c(i_OnTreatment3TC, get_i_nathist_by_age_sex(age_group_labs[i], "M", "OnTreatment3TC"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "OnTreatment3TC"))
    i_FailedTreatment3TC <- c(i_FailedTreatment3TC, get_i_nathist_by_age_sex(age_group_labs[i], "M", "FailedTreatment3TC"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "FailedTreatment3TC"))
    i_NonSevereAcute <- c(i_NonSevereAcute, get_i_nathist_by_age_sex(age_group_labs[i], "M", "NonSevereAcute"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "NonSevereAcute"))
    i_SevereAcute <- c(i_SevereAcute, get_i_nathist_by_age_sex(age_group_labs[i], "M", "SevereAcute"), get_i_nathist_by_age_sex(age_group_labs[i], "F", "SevereAcute"))
    
    i_chronic_sAgPos <- c(i_chronic_sAgPos, get_i_nathist_by_age_sex(age_group_labs[i], "F", "HBsAg_chronic"),get_i_nathist_by_age_sex(age_group_labs[i], "M", "HBsAg_chronic"))
    
    if(i==1){
        i_chronic_sAgPos_byage <- data.frame(t(c(get_i_nathist_by_age_sex(age_group_labs[i], "F", "HBsAg_chronic"),get_i_nathist_by_age_sex(age_group_labs[i], "M", "HBsAg_chronic"))))
        i_HCC_byage <- data.frame(t(c(get_i_nathist_by_age_sex(age_group_labs[i], "F", "HCC"),get_i_nathist_by_age_sex(age_group_labs[i], "M", "HCC"))))
    }else{
        i_chronic_sAgPos_byage <- rbind(i_chronic_sAgPos_byage, data.frame(t(c(get_i_nathist_by_age_sex(age_group_labs[i], "F", "HBsAg_chronic"),get_i_nathist_by_age_sex(age_group_labs[i], "M", "HBsAg_chronic")))))
        i_HCC_byage <- rbind(i_HCC_byage,data.frame(t(c(get_i_nathist_by_age_sex(age_group_labs[i], "F", "HCC"),get_i_nathist_by_age_sex(age_group_labs[i], "M", "HCC")))))
    }
}
i_all_alive <- c(i_all_alive_F, i_all_alive_M)



####################################################################################
## Now add the columns to the data frame for each:
country_files$n_alive <- rowSums(country_files[,i_all_alive])
country_files$n_aliveF <- rowSums(country_files[,i_all_alive_F])
country_files$n_aliveM <- rowSums(country_files[,i_all_alive_M])

country_files$n_HBVdeath_cumulative <- rowSums(country_files[,i_HBVdeath_cumulative])
country_files$n_Susc <- rowSums(country_files[,i_Susc])
country_files$n_ImmuneTolerant <- rowSums(country_files[,i_ImmuneTolerant])
country_files$n_ImmuneReactive <- rowSums(country_files[,i_ImmuneReactive])
country_files$i_AsymptCarrier <- rowSums(country_files[,i_AsymptCarrier])
country_files$i_ChronicHepB <- rowSums(country_files[,i_ChronicHepB])
country_files$i_CompCirr <- rowSums(country_files[,i_CompCirr])
country_files$i_DecompCirr <- rowSums(country_files[,i_DecompCirr])
country_files$i_HCC <- rowSums(country_files[,i_HCC])
country_files$i_Immune <- rowSums(country_files[,i_Immune])
country_files$i_OnTreatmentTDF <- rowSums(country_files[,i_OnTreatmentTDF])
country_files$i_OnTreatment3TC <- rowSums(country_files[,i_OnTreatment3TC])
country_files$i_FailedTreatment3TC <- rowSums(country_files[,i_FailedTreatment3TC])
country_files$i_NonSevereAcute <- rowSums(country_files[,i_NonSevereAcute])
country_files$i_SevereAcute <- rowSums(country_files[,i_SevereAcute])

country_files <- country_files[,c(1,1202:ncol(country_files))]

this.scenario <- 1
this.runnumber <- 3
plot.df.wide <- country_files[country_files$scenario %in% this.scenario & country_files$runnumber %in% 1, !(colnames(country_files) %in% c("scenario","runnumber","n_aliveF","n_aliveM"))]
plot.df.long <- plot.df.wide %>% 
  pivot_longer(cols = `n_HBVdeath_cumulative`:`i_SevereAcute`,names_to = "Stage", values_to = "value")

plot.df.long %>% 
  ggplot(aes(Year, value, fill = Stage, label = Stage, color = Stage)) +
  geom_area() 

##temp <- country_files[,colnames(country_files) %in% c("runnumber","scenario","Year","n_HBVdeath_cumulative")]

##i_chronic_sAgPos
##i_chronic_sAgPos_byage
##i_HCC_byage
## Do I need to make the number of deaths annual?
# a <- reshape(data=temp, idvar = "Year", timevar = c("runnumber","scenario"), direction = "wide")
# wide.var.names <- paste0("n_HBVdeath_cumulative.",as.character(1:3))
# b <- reshape(a, 
#         direction = "long",
#         varying = list(wide.var.names),
#         v.names = "n_HBVdeath_cumulative",
#         idvar = "Year",
#         timevar = "runnumber")
# rownames(b) <- NULL
##############
##################################################################################################
## Plots:
##################################################################################################









## Make animation plot - just for run number 1:
temp.df <- country_files[country_files$runnumber %in% 1,]
plot.df <- data.frame(matrix(nrow = 0, ncol = 3))
for(y in temp.df$Year)
{
    for(i in 1:length(age_group_labs))
    {
        temp <- data.frame(y,sum(temp.df[temp.df$Year %in% y, as.numeric(i_HCC_byage[i,])]), 
                           sum(temp.df[temp.df$Year %in% y, as.numeric(i_chronic_sAgPos_byage[i,])]), age_group_labs[i], i)
        plot.df <- rbind(plot.df,temp)
    }
    
}

colnames(plot.df) <- c("Year","N_HCC","N_chronic_sAgPos","Age_group","AgeGpNum")
plot.df$Age_group <- factor(plot.df$Age_group, levels=age_group_labs)
plot.df$GotBD <- ifelse((plot.df$Year - (plot.df$AgeGpNum*5))<=1980,"No BD", "Got BD")


plot_frame <- function(frame_number) {
    i <- (frame_number-10)/5
    ggplot(plot.df[plot.df$Year %in% (1979+frame_number),], aes(Age_group, N_HCC, fill=GotBD)) + 
        geom_bar(stat = "identity") + 
        theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size=14),
          axis.text.y = element_text(size=14),
          axis.title = element_text(size=16)) +
        scale_fill_brewer(palette="Dark2")+
    ylab("Number of HCC cases") +
    xlab("Age group") +
        ##geom_vline(xintercept=frame_number, linetype="dashed", color = "red")+
     labs(title = paste0("Year:", (1979+frame_number))

          ## plot(rnorm(100), rnorm(100),
          ##      xlim = c(-3, 3), ylim = c(-3, 3),
          ##      pch = 20, col = "blue",
          ##      main = paste0("Frame ", frame_number))
     )
}

pdf("BD_2026v2.pdf")
plot_frame(47)
dev.off()

## saveGIF({
##   for (i in 1:10) {
##     plot_frame(i)
##     ani.pause(0.5)
##   }
## }, movie.name = "random_points.gif",
##         ani.width = 480,
##         ani.height = 320)



p <- ggplot(plot.df[plot.df$Year %in% seq(1980,2100),], aes(Age_group, N_HCC, fill=GotBD)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size=14),
          axis.text.y = element_text(size=14),
          axis.title = element_text(size=16)) +
        scale_fill_brewer(palette="Dark2")+
    ylab("Number of HCC cases") +
    xlab("Age group") +
    ##geom_vline(xintercept=((plot.df$Year-10)/5), linetype="dashed", color = "red")+
    labs(title = "Year: {closest_state}")

p_animated <- p + transition_states(Year,
                                    transition_length = 2,
                                    state_length = 2) +
              enter_fade() +
              exit_fade() +
              ease_aes("linear")

anim_save("p4.gif",p_animated)



## plot.it <- function(df, i_to_plot, y.name)
## {
##     df.to.plot <- data.frame(Year=df$Year, Pop=rowSums(df[,i_to_plot]))
##     ggplot(df.to.plot, aes(x=Year, y=Pop)) +
##         geom_line() +
##         theme_bw() +
##         theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) + 
##         ylab(y.name) +
##         xlab("Year") +
##         scale_y_continuous(label=comma)
##         ##geom_bar(stat="identity")
    
## df.plot <- c




## i_Susc = 1;         ## 'Susceptible', 
## i_ImmTol = 2;       ## 'HBV: Immune Tolerant',
## i_ImmReact = 3;     ## 'HBV: Immune Reactive',
## i_AsymptCarr = 4;   ## 'HBV: Asymptomatic Carrier',
## i_Chronic = 5;      ## 'HBV: Chronic Hep B',
## i_CompCirr = 6;     ## 'HBV: Comp Cirrhosis',
## i_DecompCirr = 7;   ## 'HBV: Decomp Cirrhosis',
## i_HCC = 8;          ## 'HBV: Liver Cancer',
## i_Immume = 9;       ## 'HBV: Immune (Rec. or vacc.)',
## i_TDFtreat = 10;    ## 'HBV: TDF-Treatment',
## i_HBVdeath = 11;    ## 'Prematurely dead due to HBV',
## i_3TCtreat = 12;    ## '3TC-Treatment',
## i_3TCfailed = 13;   ## 'Failed 3TC-Treatment',
## i_NonSevAcute = 14; ## 'Non-severe acute',
## i_SevereAcute = 15; ## 'Severe acute'




#######################################################################
## Unused code:
## Number in each 5 year age group:
## Add "1" as first column is year:
i_0_4F <-  1 + unique(c(grep("Age0_4F_D[0-9]Treat", header), grep("Age0_4F_D1[02345]Treat", header), grep("Age0_4F_D[0-9]NoTreat", header), grep("Age0_4F_D1[02345]NoTreat", header)))
i_0_4M <-  1 + unique(c(grep("Age0_4M_D[0-9]Treat", header), grep("Age0_4M_D1[02345]Treat", header), grep("Age0_4M_D[0-9]NoTreat", header), grep("Age0_4M_D1[02345]NoTreat", header)))
i_5_9F <-  1 + unique(c(grep("Age5_9F_D[0-9]Treat", header), grep("Age5_9F_D1[02345]Treat", header), grep("Age5_9F_D[0-9]NoTreat", header), grep("Age5_9F_D1[02345]NoTreat", header)))
i_5_9M <-  1 + unique(c(grep("Age5_9M_D[0-9]Treat", header), grep("Age5_9M_D1[02345]Treat", header), grep("Age5_9M_D[0-9]NoTreat", header), grep("Age5_9M_D1[02345]NoTreat", header)))
i_10_14F <-  1 + unique(c(grep("Age10_14F_D[0-9]Treat", header), grep("Age10_14F_D1[02345]Treat", header), grep("Age10_14F_D[0-9]NoTreat", header), grep("Age10_14F_D1[02345]NoTreat", header)))
i_10_14M <-  1 + unique(c(grep("Age10_14M_D[0-9]Treat", header), grep("Age10_14M_D1[02345]Treat", header), grep("Age10_14M_D[0-9]NoTreat", header), grep("Age10_14M_D1[02345]NoTreat", header)))
i_15_19F <-  1 + unique(c(grep("Age15_19F_D[0-9]Treat", header), grep("Age15_19F_D1[02345]Treat", header), grep("Age15_19F_D[0-9]NoTreat", header), grep("Age15_19F_D1[02345]NoTreat", header)))
i_15_19M <-  1 + unique(c(grep("Age15_19M_D[0-9]Treat", header), grep("Age15_19M_D1[02345]Treat", header), grep("Age15_19M_D[0-9]NoTreat", header), grep("Age15_19M_D1[02345]NoTreat", header)))
i_20_24F <-  1 + unique(c(grep("Age20_24F_D[0-9]Treat", header), grep("Age20_24F_D1[02345]Treat", header), grep("Age20_24F_D[0-9]NoTreat", header), grep("Age20_24F_D1[02345]NoTreat", header)))
i_20_24M <-  1 + unique(c(grep("Age20_24M_D[0-9]Treat", header), grep("Age20_24M_D1[02345]Treat", header), grep("Age20_24M_D[0-9]NoTreat", header), grep("Age20_24M_D1[02345]NoTreat", header)))
i_25_29F <-  1 + unique(c(grep("Age25_29F_D[0-9]Treat", header), grep("Age25_29F_D1[02345]Treat", header), grep("Age25_29F_D[0-9]NoTreat", header), grep("Age25_29F_D1[02345]NoTreat", header)))
i_25_29M <-  1 + unique(c(grep("Age25_29M_D[0-9]Treat", header), grep("Age25_29M_D1[02345]Treat", header), grep("Age25_29M_D[0-9]NoTreat", header), grep("Age25_29M_D1[02345]NoTreat", header)))
i_30_34F <-  1 + unique(c(grep("Age30_34F_D[0-9]Treat", header), grep("Age30_34F_D1[02345]Treat", header), grep("Age30_34F_D[0-9]NoTreat", header), grep("Age30_34F_D1[02345]NoTreat", header)))
i_30_34M <-  1 + unique(c(grep("Age30_34M_D[0-9]Treat", header), grep("Age30_34M_D1[02345]Treat", header), grep("Age30_34M_D[0-9]NoTreat", header), grep("Age30_34M_D1[02345]NoTreat", header)))
i_35_39F <-  1 + unique(c(grep("Age35_39F_D[0-9]Treat", header), grep("Age35_39F_D1[02345]Treat", header), grep("Age35_39F_D[0-9]NoTreat", header), grep("Age35_39F_D1[02345]NoTreat", header)))
i_35_39M <-  1 + unique(c(grep("Age35_39M_D[0-9]Treat", header), grep("Age35_39M_D1[02345]Treat", header), grep("Age35_39M_D[0-9]NoTreat", header), grep("Age35_39M_D1[02345]NoTreat", header)))
i_40_44F <-  1 + unique(c(grep("Age40_44F_D[0-9]Treat", header), grep("Age40_44F_D1[02345]Treat", header), grep("Age40_44F_D[0-9]NoTreat", header), grep("Age40_44F_D1[02345]NoTreat", header)))
i_40_44M <-  1 + unique(c(grep("Age40_44M_D[0-9]Treat", header), grep("Age40_44M_D1[02345]Treat", header), grep("Age40_44M_D[0-9]NoTreat", header), grep("Age40_44M_D1[02345]NoTreat", header)))
i_45_49F <-  1 + unique(c(grep("Age45_49F_D[0-9]Treat", header), grep("Age45_49F_D1[02345]Treat", header), grep("Age45_49F_D[0-9]NoTreat", header), grep("Age45_49F_D1[02345]NoTreat", header)))
i_45_49M <-  1 + unique(c(grep("Age45_49M_D[0-9]Treat", header), grep("Age45_49M_D1[02345]Treat", header), grep("Age45_49M_D[0-9]NoTreat", header), grep("Age45_49M_D1[02345]NoTreat", header)))
i_50_54F <-  1 + unique(c(grep("Age50_54F_D[0-9]Treat", header), grep("Age50_54F_D1[02345]Treat", header), grep("Age50_54F_D[0-9]NoTreat", header), grep("Age50_54F_D1[02345]NoTreat", header)))
i_50_54M <-  1 + unique(c(grep("Age50_54M_D[0-9]Treat", header), grep("Age50_54M_D1[02345]Treat", header), grep("Age50_54M_D[0-9]NoTreat", header), grep("Age50_54M_D1[02345]NoTreat", header)))
i_55_59F <-  1 + unique(c(grep("Age55_59F_D[0-9]Treat", header), grep("Age55_59F_D1[02345]Treat", header), grep("Age55_59F_D[0-9]NoTreat", header), grep("Age55_59F_D1[02345]NoTreat", header)))
i_55_59M <-  1 + unique(c(grep("Age55_59M_D[0-9]Treat", header), grep("Age55_59M_D1[02345]Treat", header), grep("Age55_59M_D[0-9]NoTreat", header), grep("Age55_59M_D1[02345]NoTreat", header)))
i_60_64F <-  1 + unique(c(grep("Age60_64F_D[0-9]Treat", header), grep("Age60_64F_D1[02345]Treat", header), grep("Age60_64F_D[0-9]NoTreat", header), grep("Age60_64F_D1[02345]NoTreat", header)))
i_60_64M <-  1 + unique(c(grep("Age60_64M_D[0-9]Treat", header), grep("Age60_64M_D1[02345]Treat", header), grep("Age60_64M_D[0-9]NoTreat", header), grep("Age60_64M_D1[02345]NoTreat", header)))
i_65_69F <-  1 + unique(c(grep("Age65_69F_D[0-9]Treat", header), grep("Age65_69F_D1[02345]Treat", header), grep("Age65_69F_D[0-9]NoTreat", header), grep("Age65_69F_D1[02345]NoTreat", header)))
i_65_69M <-  1 + unique(c(grep("Age65_69M_D[0-9]Treat", header), grep("Age65_69M_D1[02345]Treat", header), grep("Age65_69M_D[0-9]NoTreat", header), grep("Age65_69M_D1[02345]NoTreat", header)))
i_70_74F <-  1 + unique(c(grep("Age70_74F_D[0-9]Treat", header), grep("Age70_74F_D1[02345]Treat", header), grep("Age70_74F_D[0-9]NoTreat", header), grep("Age70_74F_D1[02345]NoTreat", header)))
i_70_74M <-  1 + unique(c(grep("Age70_74M_D[0-9]Treat", header), grep("Age70_74M_D1[02345]Treat", header), grep("Age70_74M_D[0-9]NoTreat", header), grep("Age70_74M_D1[02345]NoTreat", header)))
i_75_79F <-  1 + unique(c(grep("Age75_79F_D[0-9]Treat", header), grep("Age75_79F_D1[02345]Treat", header), grep("Age75_79F_D[0-9]NoTreat", header), grep("Age75_79F_D1[02345]NoTreat", header)))
i_75_79M <-  1 + unique(c(grep("Age75_79M_D[0-9]Treat", header), grep("Age75_79M_D1[02345]Treat", header), grep("Age75_79M_D[0-9]NoTreat", header), grep("Age75_79M_D1[02345]NoTreat", header)))
i_80_84F <-  1 + unique(c(grep("Age80_84F_D[0-9]Treat", header), grep("Age80_84F_D1[02345]Treat", header), grep("Age80_84F_D[0-9]NoTreat", header), grep("Age80_84F_D1[02345]NoTreat", header)))
i_80_84M <-  1 + unique(c(grep("Age80_84M_D[0-9]Treat", header), grep("Age80_84M_D1[02345]Treat", header), grep("Age80_84M_D[0-9]NoTreat", header), grep("Age80_84M_D1[02345]NoTreat", header)))
i_85_89F <-  1 + unique(c(grep("Age85_89F_D[0-9]Treat", header), grep("Age85_89F_D1[02345]Treat", header), grep("Age85_89F_D[0-9]NoTreat", header), grep("Age85_89F_D1[02345]NoTreat", header)))
i_85_89M <-  1 + unique(c(grep("Age85_89M_D[0-9]Treat", header), grep("Age85_89M_D1[02345]Treat", header), grep("Age85_89M_D[0-9]NoTreat", header), grep("Age85_89M_D1[02345]NoTreat", header)))
i_90_94F <-  1 + unique(c(grep("Age90_94F_D[0-9]Treat", header), grep("Age90_94F_D1[02345]Treat", header), grep("Age90_94F_D[0-9]NoTreat", header), grep("Age90_94F_D1[02345]NoTreat", header)))
i_90_94M <-  1 + unique(c(grep("Age90_94M_D[0-9]Treat", header), grep("Age90_94M_D1[02345]Treat", header), grep("Age90_94M_D[0-9]NoTreat", header), grep("Age90_94M_D1[02345]NoTreat", header)))
i_95_99F <-  1 + unique(c(grep("Age95_99F_D[0-9]Treat", header), grep("Age95_99F_D1[02345]Treat", header), grep("Age95_99F_D[0-9]NoTreat", header), grep("Age95_99F_D1[02345]NoTreat", header)))
i_95_99M <-  1 + unique(c(grep("Age95_99M_D[0-9]Treat", header), grep("Age95_99M_D1[02345]Treat", header), grep("Age95_99M_D[0-9]NoTreat", header), grep("Age95_99M_D1[02345]NoTreat", header)))


###########################################
  
