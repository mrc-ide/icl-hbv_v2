rm(list=ls())

library(ggplot2)
##library(gifski)
##library(gganimate)
library(scales)
library(animation)
library(gganimate)
library(tidyr)



setwd("C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/scripts/")
current.dir <- "../outputs/"

this.country <- "GMB"
n.scenarios <- 8
n.runs <- 20

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
