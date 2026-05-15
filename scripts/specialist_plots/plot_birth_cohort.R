library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

##graph.dir <- "C:/Users/mpickles/OneDrive - Imperial College London/Dropbox_copy/Hepatits B/Presentations/15May2026/"

graph.dir <- "/mnt/c/Users/mpickles/OneDrive - Imperial College London/Dropbox_copy/Hepatits B/Presentations/15May2026/"

ISO <- "THA"
if(ISO=="GMB"){
    i.ISO <- "36"
}else if(ISO=="THA"){
    i.ISO <- "92"
}
                                                                                   
s <- 11 ## Birth cohort
i_run <- 1
df <- read.csv(paste0("/mnt/c/Users/mpickles/Documents/Hepatitis_B/icl-hbv/HPC/12may2026/",i.ISO,"/results_",ISO,"_scenario",as.character(s),"_default_run_",as.character(i_run),".csv"),header=FALSE)

header <- read.csv("/mnt/c/Users/mpickles/Documents/Hepatitis_B/icl-hbv/HPC/12may2026/36/header.txt",header=T,sep=",")
header <- colnames(header)

colnames(df) <- header

df <- df[df$Year %in% seq(2020,2100),]


##i_treateligible = [3 5 6 7]; %% Immune Reactive, Chronic, Comp+Decomp Cirr
##i_TDFtreat = 10;    % 'HBV: TDF-Treatment',


df$Age0_4F_treatelig <- df$Age0_4F_D3NoTreat + df$Age0_4F_D5NoTreat + df$Age0_4F_D6NoTreat + df$Age0_4F_D7NoTreat
df$Age5_9F_treatelig <- df$Age5_9F_D3NoTreat + df$Age5_9F_D5NoTreat + df$Age5_9F_D6NoTreat + df$Age5_9F_D7NoTreat
df$Age10_14F_treatelig <- df$Age10_14F_D3NoTreat + df$Age10_14F_D5NoTreat + df$Age10_14F_D6NoTreat + df$Age10_14F_D7NoTreat
df$Age15_19F_treatelig <- df$Age15_19F_D3NoTreat + df$Age15_19F_D5NoTreat + df$Age15_19F_D6NoTreat + df$Age15_19F_D7NoTreat
df$Age20_24F_treatelig <- df$Age20_24F_D3NoTreat + df$Age20_24F_D5NoTreat + df$Age20_24F_D6NoTreat + df$Age20_24F_D7NoTreat
df$Age25_29F_treatelig <- df$Age25_29F_D3NoTreat + df$Age25_29F_D5NoTreat + df$Age25_29F_D6NoTreat + df$Age25_29F_D7NoTreat
df$Age30_34F_treatelig <- df$Age30_34F_D3NoTreat + df$Age30_34F_D5NoTreat + df$Age30_34F_D6NoTreat + df$Age30_34F_D7NoTreat
df$Age35_39F_treatelig <- df$Age35_39F_D3NoTreat + df$Age35_39F_D5NoTreat + df$Age35_39F_D6NoTreat + df$Age35_39F_D7NoTreat
df$Age40_44F_treatelig <- df$Age40_44F_D3NoTreat + df$Age40_44F_D5NoTreat + df$Age40_44F_D6NoTreat + df$Age40_44F_D7NoTreat
df$Age45_49F_treatelig <- df$Age45_49F_D3NoTreat + df$Age45_49F_D5NoTreat + df$Age45_49F_D6NoTreat + df$Age45_49F_D7NoTreat
df$Age50_54F_treatelig <- df$Age50_54F_D3NoTreat + df$Age50_54F_D5NoTreat + df$Age50_54F_D6NoTreat + df$Age50_54F_D7NoTreat
df$Age55_59F_treatelig <- df$Age55_59F_D3NoTreat + df$Age55_59F_D5NoTreat + df$Age55_59F_D6NoTreat + df$Age55_59F_D7NoTreat
df$Age60_64F_treatelig <- df$Age60_64F_D3NoTreat + df$Age60_64F_D5NoTreat + df$Age60_64F_D6NoTreat + df$Age60_64F_D7NoTreat
df$Age65_69F_treatelig <- df$Age65_69F_D3NoTreat + df$Age65_69F_D5NoTreat + df$Age65_69F_D6NoTreat + df$Age65_69F_D7NoTreat
df$Age70_74F_treatelig <- df$Age70_74F_D3NoTreat + df$Age70_74F_D5NoTreat + df$Age70_74F_D6NoTreat + df$Age70_74F_D7NoTreat
df$Age75_79F_treatelig <- df$Age75_79F_D3NoTreat + df$Age75_79F_D5NoTreat + df$Age75_79F_D6NoTreat + df$Age75_79F_D7NoTreat
df$Age80_84F_treatelig <- df$Age80_84F_D3NoTreat + df$Age80_84F_D5NoTreat + df$Age80_84F_D6NoTreat + df$Age80_84F_D7NoTreat
df$Age85_89F_treatelig <- df$Age85_89F_D3NoTreat + df$Age85_89F_D5NoTreat + df$Age85_89F_D6NoTreat + df$Age85_89F_D7NoTreat
df$Age90_94F_treatelig <- df$Age90_94F_D3NoTreat + df$Age90_94F_D5NoTreat + df$Age90_94F_D6NoTreat + df$Age90_94F_D7NoTreat
df$Age95_99F_treatelig <- df$Age95_99F_D3NoTreat + df$Age95_99F_D5NoTreat + df$Age95_99F_D6NoTreat + df$Age95_99F_D7NoTreat




df$Age0_4M_treatelig <- df$Age0_4M_D3NoTreat + df$Age0_4M_D5NoTreat + df$Age0_4M_D6NoTreat + df$Age0_4M_D7NoTreat
df$Age5_9M_treatelig <- df$Age5_9M_D3NoTreat + df$Age5_9M_D5NoTreat + df$Age5_9M_D6NoTreat + df$Age5_9M_D7NoTreat
df$Age10_14M_treatelig <- df$Age10_14M_D3NoTreat + df$Age10_14M_D5NoTreat + df$Age10_14M_D6NoTreat + df$Age10_14M_D7NoTreat
df$Age15_19M_treatelig <- df$Age15_19M_D3NoTreat + df$Age15_19M_D5NoTreat + df$Age15_19M_D6NoTreat + df$Age15_19M_D7NoTreat
df$Age20_24M_treatelig <- df$Age20_24M_D3NoTreat + df$Age20_24M_D5NoTreat + df$Age20_24M_D6NoTreat + df$Age20_24M_D7NoTreat
df$Age25_29M_treatelig <- df$Age25_29M_D3NoTreat + df$Age25_29M_D5NoTreat + df$Age25_29M_D6NoTreat + df$Age25_29M_D7NoTreat
df$Age30_34M_treatelig <- df$Age30_34M_D3NoTreat + df$Age30_34M_D5NoTreat + df$Age30_34M_D6NoTreat + df$Age30_34M_D7NoTreat
df$Age35_39M_treatelig <- df$Age35_39M_D3NoTreat + df$Age35_39M_D5NoTreat + df$Age35_39M_D6NoTreat + df$Age35_39M_D7NoTreat
df$Age40_44M_treatelig <- df$Age40_44M_D3NoTreat + df$Age40_44M_D5NoTreat + df$Age40_44M_D6NoTreat + df$Age40_44M_D7NoTreat
df$Age45_49M_treatelig <- df$Age45_49M_D3NoTreat + df$Age45_49M_D5NoTreat + df$Age45_49M_D6NoTreat + df$Age45_49M_D7NoTreat
df$Age50_54M_treatelig <- df$Age50_54M_D3NoTreat + df$Age50_54M_D5NoTreat + df$Age50_54M_D6NoTreat + df$Age50_54M_D7NoTreat
df$Age55_59M_treatelig <- df$Age55_59M_D3NoTreat + df$Age55_59M_D5NoTreat + df$Age55_59M_D6NoTreat + df$Age55_59M_D7NoTreat
df$Age60_64M_treatelig <- df$Age60_64M_D3NoTreat + df$Age60_64M_D5NoTreat + df$Age60_64M_D6NoTreat + df$Age60_64M_D7NoTreat
df$Age65_69M_treatelig <- df$Age65_69M_D3NoTreat + df$Age65_69M_D5NoTreat + df$Age65_69M_D6NoTreat + df$Age65_69M_D7NoTreat
df$Age70_74M_treatelig <- df$Age70_74M_D3NoTreat + df$Age70_74M_D5NoTreat + df$Age70_74M_D6NoTreat + df$Age70_74M_D7NoTreat
df$Age75_79M_treatelig <- df$Age75_79M_D3NoTreat + df$Age75_79M_D5NoTreat + df$Age75_79M_D6NoTreat + df$Age75_79M_D7NoTreat
df$Age80_84M_treatelig <- df$Age80_84M_D3NoTreat + df$Age80_84M_D5NoTreat + df$Age80_84M_D6NoTreat + df$Age80_84M_D7NoTreat
df$Age85_89M_treatelig <- df$Age85_89M_D3NoTreat + df$Age85_89M_D5NoTreat + df$Age85_89M_D6NoTreat + df$Age85_89M_D7NoTreat
df$Age90_94M_treatelig <- df$Age90_94M_D3NoTreat + df$Age90_94M_D5NoTreat + df$Age90_94M_D6NoTreat + df$Age90_94M_D7NoTreat
df$Age95_99M_treatelig <- df$Age95_99M_D3NoTreat + df$Age95_99M_D5NoTreat + df$Age95_99M_D6NoTreat + df$Age95_99M_D7NoTreat

dim(df)
# 81 1293

df <- df %>%
  select(-contains(c("D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D11", "D12", "D13", "D14", "D15","D10Treat","Death","N_","NBD","DALY","NBirth","D1Treat","D1NoTreat","Incid")))
dim(df)
## 81 142




long_df <- pivot_longer(df, 
                                        #cols = ends_with("_treatelig"),
                        cols = starts_with("Age"), 
                        names_to = "VarName", 
                        values_to = "Value")


dim(long_df)
## 6480 3


long_df$AgeGp <- gsub("D10NoTreat","",long_df$VarName)
long_df$AgeGp <- gsub("F_","",long_df$AgeGp)
long_df$AgeGp <- gsub("M_","",long_df$AgeGp)
long_df$AgeGp <- gsub("treatelig","",long_df$AgeGp)
long_df$AgeGp <- gsub("Age","",long_df$AgeGp)
long_df$AgeGp <- gsub("_","-",long_df$AgeGp)

##age.gps <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95-99")
age.gps <- c("10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95-99")
long_df <- long_df[long_df$AgeGp %in% age.gps,]
long_df$AgeGp <- factor(long_df$AgeGp, levels = age.gps)

long_df$sex <- ifelse(grepl("F", long_df$VarName), "F", "M")
long_df$Treat <- ifelse(grepl("treatelig", long_df$VarName), "Eligible, not on treatment", "On treatment")
long_df$Treat <- factor(long_df$Treat, levels = c("On treatment","Eligible, not on treatment"))


## Get default colours and reverse them
library(scales)
hex <- rev(hue_pal()(2))

for(i in 2025:2040)
{
    pdf(paste0(graph.dir,"birthcohort_treatcoverage",ISO,"_",as.character(i),".pdf"))
    ##print(
    p <- ggplot(long_df[long_df$Year %in% i,], aes(fill=Treat, y=Value, x=AgeGp)) + 
        geom_bar(position="stack", stat="identity") +
        theme_bw() +
        xlab("") +
        ylab("Number of people") +
        ggtitle(paste0("T=",as.character(i))) +
        theme(legend.position = "bottom",axis.text=element_text(size=14),
              axis.title=element_text(size=16,face="bold"), legend.text = element_text(size=14),
              axis.text.x = element_text(angle = 90), 
              legend.title=element_text(size=14, face="bold"))+
        scale_fill_manual(values = hex)
    
    if(ISO=="GMB")
    {
        p <- p + ylim(0,1600)
    }else if(ISO=="THA")
    {
        p <- p + ylim(0,8500)
    }

    print(p)
    dev.off()
}
