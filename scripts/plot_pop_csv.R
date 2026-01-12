rm(list=ls())



current.dir <- "../outputs/"
n.runs <- 1


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
loadfiles <- function(i,filedir, header)
{
    output <- read.csv(paste0(filedir,"results_THA_scenario1_default_run_",as.character(i),".csv"),header=FALSE)
    ## Last column is -1
    output <- output[1:(nrow(output)-1),]
    ## Add column names:
    colnames(output) <- c("Year",header)
    
    return(output)
}


##
make_average <- function(n.files,filedir, header)
{

    output <- loadfiles(1,filedir, header)
    
    if(n.files>1)
    {
        for(i in 2:n.files)
        {
            infile <- loadfiles(i,filedir, header)
            output <- output + infile
        }
        output <- output/n.files
    }
    
    return(output)
}


header <- get_header(current.dir)
if(is.na(header[length(header)]))
{
    header <- header[1:(length(header)-1)]
}
currentfiles <- loadfiles(n.runs,current.dir, header)


##t <- currentfiles[[1]]$Year

## Add "1" as first column is year:
i_0_4F <-  1 + unique(c(grep("Age0_4F_D[0-9]Treat", header), grep("Age0_4F_D1[02345]Treat", header), grep("Age0_4F_D[0-9]NoTreat", header), grep("Age0_4F_D1[02345]NoTreat", header)))
i_0_4M <-  1 + unique(c(grep("Age0_4M_D[0-9]Treat", header), grep("Age0_4M_D1[02345]Treat", header), grep("Age0_4M_D[0-9]NoTreat", header), grep("Age0_4M_D1[02345]NoTreat", header)))

## For a specific age and sex group (summed over Treat/No treat) get the number in a certain natural history class (e.g. alive)
get_nathist_by_age_sex <- function(age_gp, sex, nat_hist_state){

    if(nat_hist_state=="Alive")
    {
        i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D[0-9]Treat"), header), grep(paste0("Age",age_gp,sex,"_D1[02345]Treat"), header), grep(paste0("Age",age_gp,sex,"_D[0-9]NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D1[02345]NoTreat"), header)))
    }else if(nat_hist_state=="Dead")
    {
        i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D[0-9]Treat"), header), grep(paste0("Age",age_gp,sex,"_D1[02345]Treat"), header), grep(paste0("Age",age_gp,sex,"_D[0-9]NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D1[02345]NoTreat"), header)))
    }else if(nat_hist_state=="Chronic")
    {
        i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D[0-9]Treat"), header), grep(paste0("Age",age_gp,sex,"_D1[02345]Treat"), header), grep(paste0("Age",age_gp,sex,"_D[0-9]NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D1[02345]NoTreat"), header)))
    }else if(nat_hist_state=="HBsAg_incacute")
    {
        i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D[2-8]Treat"), header), grep(paste0("Age",age_gp,sex,"_D10Treat"), header), grep(paste0("Age",age_gp,sex,"_D1[2-5]Treat"), header), grep(paste0("Age",age_gp,sex,"_D[2-8]NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D10NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D1[2-5]NoTreat"), header)))
    }else if(nat_hist_state=="HBsAg_chronic")
    {
        i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D[2-8]Treat"), header), grep(paste0("Age",age_gp,sex,"_D10Treat"), header), grep(paste0("Age",age_gp,sex,"_D1[2-3]Treat"), header), grep(paste0("Age",age_gp,sex,"_D[2-8]NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D10NoTreat"), header), grep(paste0("Age",age_gp,sex,"_D1[2-3]NoTreat"), header)))
    }else if(nat_hist_state=="HCC")
    {
        i <- 1 + unique(c(grep(paste0("Age",age_gp,sex,"_D8Treat"), header), grep(paste0("Age",age_gp,sex,"_D8NoTreat"), header)))
    }



    ##[2:8 10 12:15]
    return(i)
}


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
i_chronic_sAgPos <- c()
i_HCC <- c()


## for(i in length(age_group_labs))
## {
##     i_chronic_sAgPos_byage <- append(i_chronic_sAgPos_byage,c())
##     i_HCC_byage <- append(i_HCC_byage,c())
## }
for(i in 1:length(age_group_labs))
{
    i_all_alive_F <- c(i_all_alive_F,get_nathist_by_age_sex(age_group_labs[i], "F", "Alive"))
    i_all_alive_M <- c(i_all_alive_M,get_nathist_by_age_sex(age_group_labs[i], "M", "Alive"))
    i_chronic_sAgPos <- c(i_chronic_sAgPos, get_nathist_by_age_sex(age_group_labs[i], "F", "HBsAg_chronic"),get_nathist_by_age_sex(age_group_labs[i], "M", "HBsAg_chronic"))
    i_HCC <- c(i_HCC, get_nathist_by_age_sex(age_group_labs[i], "F", "HCC"),get_nathist_by_age_sex(age_group_labs[i], "M", "HCC"))
    if(i==1){
        i_chronic_sAgPos_byage <- data.frame(t(c(get_nathist_by_age_sex(age_group_labs[i], "F", "HBsAg_chronic"),get_nathist_by_age_sex(age_group_labs[i], "M", "HBsAg_chronic"))))
        i_HCC_byage <- data.frame(t(c(get_nathist_by_age_sex(age_group_labs[i], "F", "HCC"),get_nathist_by_age_sex(age_group_labs[i], "M", "HCC"))))
    }else{
        i_chronic_sAgPos_byage <- rbind(i_chronic_sAgPos_byage, data.frame(t(c(get_nathist_by_age_sex(age_group_labs[i], "F", "HBsAg_chronic"),get_nathist_by_age_sex(age_group_labs[i], "M", "HBsAg_chronic")))))
        i_HCC_byage <- rbind(i_HCC_byage,data.frame(t(c(get_nathist_by_age_sex(age_group_labs[i], "F", "HCC"),get_nathist_by_age_sex(age_group_labs[i], "M", "HCC")))))
    }
}

df <- currentfiles

plot.df <- data.frame(matrix(nrow = 0, ncol = 2))
for(y in currentfiles$Year)
{
    temp <- c(y)
    for(i in 1:length(age_group_labs))
    {
        temp <- c(temp,sum(df[df$Year %in% y, as.numeric(i_HCC_byage[i,])]))
    }
    plot.df <- rbind(plot.df,temp)
}
colnames(plot.df) <- c("Year",age_group_labs)

i_all_alive <- c(i_all_alive_F, i_all_alive_M)

library(ggplot2)
library(scales)




p <- ggplot(df, aes(category, value, fill = category)) + 
     geom_bar(stat = "identity") + 
     scale_fill_manual(values = c("#FF9999",
                                  "#66CCCC",
                                  "#FFCC66")) + 
     theme_minimal() +
     labs(title = "Year: {closest_state}")

p_animated <- p + transition_states(year,
                                    transition_length = 1,
                                    state_length = 1) +
              enter_fade() +
              exit_fade() +
              ease_aes("linear")

p_animated

plot.it <- function(df, i_to_plot, y.name)
{
    df.to.plot <- data.frame(Year=df$Year, Pop=rowSums(df[,i_to_plot]))
    ggplot(df.to.plot, aes(x=Year, y=Pop)) +
        geom_line() +
        theme_bw() +
        theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) + 
        ylab(y.name) +
        xlab("Year") +
        scale_y_continuous(label=comma)
        ##geom_bar(stat="identity")
    
df.plot <- c




i_Susc = 1;         ## 'Susceptible', 
i_ImmTol = 2;       ## 'HBV: Immune Tolerant',
i_ImmReact = 3;     ## 'HBV: Immune Reactive',
i_AsymptCarr = 4;   ## 'HBV: Asymptomatic Carrier',
i_Chronic = 5;      ## 'HBV: Chronic Hep B',
i_CompCirr = 6;     ## 'HBV: Comp Cirrhosis',
i_DecompCirr = 7;   ## 'HBV: Decomp Cirrhosis',
i_HCC = 8;          ## 'HBV: Liver Cancer',
i_Immume = 9;       ## 'HBV: Immune (Rec. or vacc.)',
i_TDFtreat = 10;    ## 'HBV: TDF-Treatment',
i_HBVdeath = 11;    ## 'Prematurely dead due to HBV',
i_3TCtreat = 12;    ## '3TC-Treatment',
i_3TCfailed = 13;   ## 'Failed 3TC-Treatment',
i_NonSevAcute = 14; ## 'Non-severe acute',
i_SevereAcute = 15; ## 'Severe acute'

