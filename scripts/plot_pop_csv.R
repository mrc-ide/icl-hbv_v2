rm(list=ls())

library(ggplot2)
##library(gifski)
##library(gganimate)
library(scales)
library(animation)
library(gganimate)


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
i_all_alive <- c(i_all_alive_F, i_all_alive_M)



df <- currentfiles

plot.df <- data.frame(matrix(nrow = 0, ncol = 3))
for(y in currentfiles$Year)
{
    for(i in 1:length(age_group_labs))
    {
        temp <- data.frame(y,sum(df[df$Year %in% y, as.numeric(i_HCC_byage[i,])]), sum(df[df$Year %in% y, as.numeric(i_chronic_sAgPos_byage[i,])]), age_group_labs[i], i)
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

