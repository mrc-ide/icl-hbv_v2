## Takes WUENIC data (vaccination) and GHO treatment (supplemented by Polaris for a few countries*)
## Generates two csv files BDdatatemp.csv (!!!109 countries!!!) and treatment_cascade.csv (110 countries) 
## which are (1 row=1 country) datasets giving vaccination/% rural/% in-facility births, and % diagnosed/on treatment if diagnosed.
## These have been used in generating slides for EASL 2026, and treatment_cascade.csv will be used in modelling.

## *Neither GHO nor Polaris cover all 110 modelled countries.
## WUENIC missing Nicaragua. 
rm(list=ls())
library(stringi)
library(readxl)
library(tibble)
library(ggpubr)
library("gridExtra")
library(dplyr)

basefolder <- "C:/Users/mpickles/"
graph.dir <- "C:/Users/mpickles/OneDrive - Imperial College London/Dropbox_copy/Hepatitis B/Data/"
##basefolder <- "/mnt/c/Users/mpickles/Documents/"

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

################################################
## Load
wuenic.data.folder <- paste0(basefolder,'Documents/Hepatitis_B/icl-hbv/resources/wuenic2024data')
all_files <- list.files(path = wuenic.data.folder, full.names = TRUE)
all_filenames <- gsub(paste0(wuenic.data.folder,'/'),'',all_files)

hepb3_files <- all_filenames[stri_detect_fixed(all_filenames, "HepB3")]
hepbb_files <- all_filenames[stri_detect_fixed(all_filenames, "HepBB")]

CountryClassification <- read.csv("C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/resources/CountryClassifier.csv",header=TRUE)
CountryClassification <- CountryClassification[,colnames(CountryClassification) %in% c("ISO","WHO_region")]


## Check files all match:
##a <- gsub('HepB3','',hepb3_files)
##b <- gsub('HepBB','',hepbb_files)
##setdiff(a,b)


list_of_wuenic_countries <- gsub(' - .csv','',gsub('wuenic2024revision_coverage-data_','',gsub('HepB3','',hepb3_files)))

country_isos <- read.csv(paste0(basefolder,'Documents/Hepatitis_B/icl-hbv/resources/wuenic2024data/ListOfCountryNamesAndISOS.csv'),header=FALSE)
colnames(country_isos) <- c("ISO","Country_name")
## Remove leading whitespace:
country_isos$Country_name <- trimws(country_isos$Country_name)

country_isos$Country_name[country_isos$Country_name %in% "The former Yugoslav Republic of Macedonia"] <- "North Macedonia"
country_isos$Country_name[country_isos$Country_name %in% "Swaziland"] <- "Eswatini"

## Check I have all the countries:
found <- c()
not_found <- c()
for(i in 1:nrow(country_isos))
{
    c <- country_isos$Country_name[i]
    ##if(any(stri_detect_fixed(list_of_wuenic_countries,c))){
    if(c %in% list_of_wuenic_countries){
        found <- c(found,c)
    }else{
        not_found <- c(not_found,c)
    }
}

##not_found (I've edited North Macedonia and Eswatini above by hand so they are no longer missing).
##[1] "The former Yugoslav Republic of Macedonia" - now "North Macedonia"
##[2] "Nicaragua"
##[3] "Swaziland" - Eswatini


found <- c()
not_found <- c()
for(i in 1:length(list_of_wuenic_countries))
{
    c <- list_of_wuenic_countries[i]
    if(c %in% country_isos$Country_name){
        found <- c(found,c)
    }else{
        not_found <- c(not_found,c)
    }
}


## We use the 109 countries modelled excluding Nicaragua (as can't find WUENIC data for Nicaragua).
list_of_countries_to_use <- found


########################################################################

## Read in a transposed dataset:
read.tcsv = function(file, header=TRUE, sep=",", ...) {

  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)

  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }

  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep)
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
}


BD_model_params <- read.tcsv(paste0(basefolder,"Documents/Hepatitis_B/icl-hbv/src/raw_params/BD_table.csv"))
colnames(BD_model_params)[which(names(BD_model_params) == "ISO_code")] <- "Year"
BD_model_params$Year <- as.integer(gsub('x','',BD_model_params$Year))

Hep3_model_params <- read.tcsv(paste0(basefolder,"Documents/Hepatitis_B/icl-hbv/src/raw_params/HepB3_table.csv"))
colnames(Hep3_model_params)[which(names(Hep3_model_params) == "ISO_code")] <- "Year"
Hep3_model_params$Year <- as.integer(gsub('x','',Hep3_model_params$Year))

BD_model_params <- BD_model_params[BD_model_params$Year %in% seq(2000,2019),]
Hep3_model_params <- Hep3_model_params[Hep3_model_params$Year %in% seq(2000,2019),]


bd.data <- data.frame()
hepB3.data <- data.frame()
for (country in list_of_countries_to_use)
{
    country.ISO <- country_isos$ISO[country_isos$Country_name %in% country]
    ## Read each country's data, but just keep Country, Antigen, Year, WUENIC [estimate of coverage in %]
    bd.country.data <- read.csv(paste0(wuenic.data.folder,'/wuenic2024revision_coverage-data_',country,' - HepBB.csv'))[,seq(1,4)]
    if(nrow(bd.country.data)==0){
        bd.country.data[1,] <- c(country,"HepBB",2024,0)
    }
    bd.country.data$ISO <- country.ISO
    
    hepB3.country.data <- read.csv(paste0(wuenic.data.folder,'/wuenic2024revision_coverage-data_',country,' - HepB3.csv'))[,seq(1,4)]
    if(nrow(hepB3.country.data)==0){
        hepB3.country.data[1,] <- c(country,"HepBB",2024,0)
    }
    hepB3.country.data$ISO <- country.ISO
    
    bd.data <- rbind(bd.data, bd.country.data)
    hepB3.data <- rbind(hepB3.data, hepB3.country.data)
}

hepB3.data <- merge(hepB3.data,CountryClassification,by="ISO")
bd.data <- merge(bd.data,CountryClassification,by="ISO")
bd.data$Year <- as.numeric(bd.data$Year)
bd.data$WUENIC <- as.numeric(bd.data$WUENIC)

list.of.regions <- c("Africa","Americas","Eastern Mediterranean","Europe","South-East Asia","Western Pacific")
#r <- "Europe"
#r <- "Eastern Mediterranean"

for(r in 1:length(list.of.regions)){
  region.nospaces <- gsub(" ","_",list.of.regions[r])
  pdf(paste0(graph.dir,"WUENIC/HepB3_timetrends_panel_",region.nospaces,".pdf"),width=20,height=12)
  p <- ggplot(hepB3.data[hepB3.data$WHO_region %in% list.of.regions[r] ,], aes(x=Year, y=WUENIC)) +
    geom_line(aes(color=ISO))+
    geom_point(aes(color=ISO)) +
    ##facet_wrap(~Country_name) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.title = element_text(size=20),
          axis.text = element_text(size=18),
          plot.title = element_text(size=32, face="bold"),
          strip.text = element_text(size=12),
          panel.spacing = unit(2, "lines"),
          legend.title = element_blank(),
          legend.text = element_text(size=12)) +
    xlab("") +
    ylab("HepB3 coverage (%)") +
    ##coord_cartesian(ylim=c(0,50)) +
    ggtitle(list.of.regions[r])
  print(p)
  dev.off()
}

## Finding some outliers:
hepB3.data[hepB3.data$Year %in% seq(2010,2024) & hepB3.data$WUENIC<50,]
hepB3.data[hepB3.data$Year %in% seq(2010,2024) & hepB3.data$WUENIC<50 & hepB3.data$WHO_region %in% "Europe",]

hepB3.data[hepB3.data$Year==2024 & hepB3.data$WUENIC<60 & hepB3.data$WHO_region %in% "Eastern Mediterranean",]
hepB3.data[hepB3.data$Year==2021 & hepB3.data$WUENIC<50 ,]
hepB3.data[hepB3.data$Year %in% seq(2010,2024) & hepB3.data$ISO %in% "BIH",]
hepB3.data[hepB3.data$Year %in% seq(2010,2024) & hepB3.data$ISO %in% "SDN",]


## BD:
for(r in 1:length(list.of.regions)){
  region.nospaces <- gsub(" ","_",list.of.regions[r])
  pdf(paste0(graph.dir,"WUENIC/HepBD_timetrends_panel_",region.nospaces,".pdf"),width=20,height=12)
  p <- ggplot(bd.data[bd.data$WHO_region %in% list.of.regions[r] ,], aes(x=Year, y=WUENIC)) +
    geom_line(aes(color=ISO))+
    geom_point(aes(color=ISO)) +
    ##facet_wrap(~Country_name) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.title = element_text(size=20),
          axis.text = element_text(size=18),
          plot.title = element_text(size=32, face="bold"),
          strip.text = element_text(size=12),
          panel.spacing = unit(2, "lines"),
          legend.title = element_blank(),
          legend.text = element_text(size=12)) +
    xlab("") +
    ylab("HepB BD coverage (%)") +
    ##coord_cartesian(ylim=c(0,50)) +
    ggtitle(list.of.regions[r])
  print(p)
  dev.off()
}

bd.data[bd.data$Year==2024 & bd.data$WUENIC<60 & bd.data$WHO_region %in% "Eastern Mediterranean",]
a <- bd.data[bd.data$Year==2024 & bd.data$WUENIC<75 & bd.data$WUENIC>0,]
a[order(a$WHO_region),]
## 34 Gavi-eligible countries:
bd.data[bd.data$Year %in% 2024 & bd.data$ISO %in% 
          c("BGD","BDI","CMR","CAF","TCD","COM","COG","COD","ERI","ETH","GHA","GIN",
            "GNB","HTI","KEN","LSO","LBR","MDG","MLI","MWI","MOZ","NPL","NER",
            "RWA","SLE","SOM","SDN","SSD","TZA","TGO","UGA","YEM","ZMB","ZWE"),]

bd.data[bd.data$Year %in% 2024 & bd.data$WUENIC %in% 0 & !(bd.data$ISO %in% 
                                                             c("BGD","BDI","CMR","CAF","TCD","COM","COG","COD","ERI","ETH","GHA","GIN",
                                                               "GNB","HTI","KEN","LSO","LBR","MDG","MLI","MWI","MOZ","NPL","NER",
                                                               "RWA","SLE","SOM","SDN","SSD","TZA","TGO","UGA","YEM","ZMB","ZWE")),]
##9 countries not currently GAVI-eligible:
# 12   AGO                           Angola   HepBB 2024      0                Africa
# 131  BOL Bolivia (Plurinational State of)   HepBB 2024      0              Americas
# 482  JAM                          Jamaica   HepBB 2024      0              Americas
# 483  JOR                           Jordan   HepBB 2024      0 Eastern Mediterranean
# 562  LKA                        Sri Lanka   HepBB 2024      0       South-East Asia
# 684  MRT                       Mauritania   HepBB 2024      0                Africa
# 801  PRY                         Paraguay   HepBB 2024      0              Americas
# 877  SWZ                         Eswatini   HepBB 2024      0                Africa
# 1139 ZAF                     South Africa   HepBB 2024      0                Africa

bd.data$ISO[bd.data$Year %in% 2024 & bd.data$WUENIC %in% 0 & !(bd.data$ISO %in% 
                                                             c("BGD","BDI","CMR","CAF","TCD","COM","COG","COD","ERI","ETH","GHA","GIN",
                                                               "GNB","HTI","KEN","LSO","LBR","MDG","MLI","MWI","MOZ","NPL","NER",
                                                               "RWA","SLE","SOM","SDN","SSD","TZA","TGO","UGA","YEM","ZMB","ZWE"))]
## List of ISOs: "AGO" "BOL" "JAM" "JOR" "LKA" "MRT" "PRY" "SWZ" "ZAF"
################################################################
## Get the UN WUP data on % rural pop by country each year:
UN_rural_percent_rawdata <- read_excel(paste0(basefolder,"Documents/Hepatitis_B/icl-hbv/resources/wuenic2024data/UN_WUP2025_ruralpercentage_bycountry.xlsx"),sheet="Data")
## We want three things from this:
## 1) The current (2024 - to correspond to wuenic) % rural
## 2) The continent/region/country (location)/ISO3_Code so I can merge this with others.
## 3) Time trend % rural (in case I want to plot it) - as transpose.

## For 3), we ditch the first columns (we want to transpose, using country ISO as the column name, so each column is a country's rural %)
temp <- UN_rural_percent_rawdata[,c(4,6:ncol(UN_rural_percent_rawdata))]

UN_rural_percent_data <- UN_rural_percent_rawdata[,colnames(UN_rural_percent_rawdata) %in% c("Continent","Region","Location","ISO3_Code","2024")]
colnames(UN_rural_percent_data) <- c("Continent","Region","Country_UN","ISO","Rural_percent")

timetrend_rural_percent <- as_tibble(cbind(Year = names(temp)[2:ncol(temp)], t(temp[,2:ncol(temp)])))
colnames(timetrend_rural_percent) <- c("Year",temp$ISO3_Code)





GHOpercentdeliveredinhealthfacility_alldata <- read_excel(paste0(basefolder,"Documents/Hepatitis_B/icl-hbv/resources/wuenic2024data/GHO_prop_births_deliveredinhealthfaility_download5Feb2026.xlsx"),sheet="Data")
GHOpercentdeliveredinhealthfacility_alldata <- GHOpercentdeliveredinhealthfacility_alldata[, colnames(GHOpercentdeliveredinhealthfacility_alldata) %in% c("SpatialDimValueCode","Period","Value","IsLatestYear","ParentLocation")]
## Only take most recent year:
GHOpercentdeliveredinhealthfacility_rawdata <- GHOpercentdeliveredinhealthfacility_alldata[GHOpercentdeliveredinhealthfacility_alldata$IsLatestYear %in% TRUE, colnames(GHOpercentdeliveredinhealthfacility_alldata) %in% c("SpatialDimValueCode","Period","Value")]
colnames(GHOpercentdeliveredinhealthfacility_rawdata) <- c("ISO","Period","Percent_deliveredInHealthFacility")
## Pakistan turns up twice still in the GHO data, remove the older (2016-2020) estimate (they are the same %).
GHOpercentdeliveredinhealthfacility_rawdata <- GHOpercentdeliveredinhealthfacility_rawdata[!(GHOpercentdeliveredinhealthfacility_rawdata$ISO %in% "PAK" & GHOpercentdeliveredinhealthfacility_rawdata$Period %in% "2016-2020"),]

## From "Period" which can be a range of years, convert to a year (if a range, take the mean and round to the nearest year):
GHOpercentdeliveredinhealthfacility_rawdata$Year <- as.integer(GHOpercentdeliveredinhealthfacility_rawdata$Period)
i.to.fix <- which(is.na(GHOpercentdeliveredinhealthfacility_rawdata$Year))
for (i in i.to.fix){
    ## round() in r rounds to the closest *even* number! Which is stupid from a maths point of view.
    ## So floor(x+0.5) will round correctly
    GHOpercentdeliveredinhealthfacility_rawdata$Year[i] <- floor(mean(lapply(strsplit(GHOpercentdeliveredinhealthfacility_rawdata$Period[i],"-"), as.integer)[[1]])+0.5)
}


#######################################################################

## Pull out the most recent value:
wuenic.bd.coverage.2024 <- bd.data[bd.data$Year %in% 2024,]
colnames(wuenic.bd.coverage.2024)[which(colnames(wuenic.bd.coverage.2024) == "WUENIC")] <- "WUENIC2024"
##wuenic.bd.coverage.2024 <- wuenic.bd.coverage.2024[,!(colnames(wuenic.bd.coverage.2024) %in% "Year")]
## Now merge with the rural pervcentage:
temp <- merge(x=wuenic.bd.coverage.2024[,!(colnames(wuenic.bd.coverage.2024) %in% "Year")], y=UN_rural_percent_data, by="ISO", all.x=TRUE)

bd.coverage.2024 <- merge(x=temp, y=GHOpercentdeliveredinhealthfacility_rawdata, by="ISO", all.x=TRUE)


list.of.isos <- unique(bd.data$ISO)
df.date.first.bd <- data.frame(ISO=character(), year_first_bd=integer(), stringsAsFactors = FALSE)
for (this.iso in list.of.isos)
{
  tempdf <- as.numeric(bd.data$Year[bd.data$ISO %in% this.iso & bd.data$WUENIC>0])
  if(length(tempdf)==0)
    df.date.first.bd[nrow(df.date.first.bd) + 1,] = c(this.iso,NA)
  else
    df.date.first.bd[nrow(df.date.first.bd) + 1,] = c(this.iso,min(tempdf))
}

bd.coverage.2024 <- merge(bd.coverage.2024, df.date.first.bd, by="ISO")




##write.csv(bd.coverage.2024,paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatitis B/Presentations/15May2026/BDdatatemp.csv"), row.names = FALSE)
write.csv(bd.coverage.2024,paste0(basefolder,"Documents/Hepatitis_B/icl-hbv/resources/wuenic2024data/BDdata_infacility_rural.csv"), row.names = FALSE)



#######################################################################################
## Now extract GHO treatment data:


Prop_ontreat_of_diagnosed <- read_excel(paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatitis B/Data/GHO_treatment/Prop_ontreat_of_diagnosed.xlsx"),sheet="Prop_ontreat_of_diagnosed")

Prop_diagnosed <- read_excel(paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatitis B/Data/GHO_treatment/Prop_diagnosed.xlsx"),sheet="Prop_diagnosed")

Prop_diagnosed <- Prop_diagnosed[Prop_diagnosed$Indicator %in% "Chronic hepatitis B (HBV) diagnosis coverage as a percentage of infected (%)" & Prop_diagnosed$IsLatestYear,]
Prop_ontreat_of_diagnosed <- Prop_ontreat_of_diagnosed[Prop_ontreat_of_diagnosed$Indicator %in% "Chronic hepatitis B (HBV) treatment rate as percentage of diagnosed (%)" & Prop_ontreat_of_diagnosed$IsLatestYear,]

dim(Prop_diagnosed[Prop_diagnosed$SpatialDimValueCode %in% list.of.isos,])  ## 97
dim(Prop_ontreat_of_diagnosed[Prop_ontreat_of_diagnosed$SpatialDimValueCode %in% list.of.isos,]) ## 103


Prop_diagnosed <- Prop_diagnosed[Prop_diagnosed$SpatialDimValueCode %in% list.of.isos, colnames(Prop_diagnosed) %in% c("SpatialDimValueCode","Location","Period","Value")]
Prop_ontreat_of_diagnosed <- Prop_ontreat_of_diagnosed[Prop_ontreat_of_diagnosed$SpatialDimValueCode %in% list.of.isos, colnames(Prop_ontreat_of_diagnosed) %in% c("SpatialDimValueCode","Location","Period","Value")]

names(Prop_diagnosed)[names(Prop_diagnosed) == 'SpatialDimValueCode'] <- 'ISO'
names(Prop_diagnosed)[names(Prop_diagnosed) == 'Location'] <- 'Country_name'
names(Prop_diagnosed)[names(Prop_diagnosed) == 'Period'] <- 'Year'
names(Prop_diagnosed)[names(Prop_diagnosed) == 'Value'] <- 'Diagnosed'
## Convert from % to proportion (to be consistent with Polaris)
Prop_diagnosed$Diagnosed <- Prop_diagnosed$Diagnosed / 100

names(Prop_ontreat_of_diagnosed)[names(Prop_ontreat_of_diagnosed) == 'SpatialDimValueCode'] <- 'ISO'
names(Prop_ontreat_of_diagnosed)[names(Prop_ontreat_of_diagnosed) == 'Location'] <- 'Country_name'
names(Prop_ontreat_of_diagnosed)[names(Prop_ontreat_of_diagnosed) == 'Period'] <- 'Year'
names(Prop_ontreat_of_diagnosed)[names(Prop_ontreat_of_diagnosed) == 'Value'] <- 'OnTreatOfDiagnosed'

## Convert from % to proportion (to be consistent with Polaris)
Prop_ontreat_of_diagnosed$OnTreatOfDiagnosed <- Prop_ontreat_of_diagnosed$OnTreatOfDiagnosed / 100

## Check for missing diagnosed:
missing.isos <- setdiff(list.of.isos,Prop_diagnosed$ISO)
## Add Nicaragua if needed:
if(!("NIC" %in% list.of.isos))
  missing.isos <- c(missing.isos,"NIC")
## [1] "AZE" "BGD" "CHN" "PRK" "IND" "KIR" "MHL" "FSM" "MMR" "WSM" "TON" "TUV"
## Taken from Polaris 2025 estimates:
Polaris2025 <- read_excel(paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatitis B/Data/Polaris/Polaris Database Query – CDA Foundation.xlsx"),sheet="Cleaned")
Polaris_propdiag <- Polaris2025[Polaris2025$ISO %in% missing.isos, colnames(Polaris2025) %in% c("Country_name","ISO","Diagnosed")]
Polaris_propdiag$Year <- 2025 ## All 2025 estimates
Polaris_propdiag <- Polaris_propdiag[,c(2,1,4,3)]
## Fix by hand to make countries consistent:
Polaris_propdiag$Country_name[Polaris_propdiag$ISO %in% "FSM"] <- "Micronesia (Federated States of)"

## Merge with previous:
Prop_diagnosed <- rbind(Prop_diagnosed,Polaris_propdiag)


##
## Check for missing on treat if diagnosed:
missing.isos <- setdiff(list.of.isos,Prop_ontreat_of_diagnosed$ISO)
## Add Nicaragua if needed:
if(!("NIC" %in% list.of.isos))
  missing.isos <- c(missing.isos,"NIC")
## [1] "AZE" "BGD" "CHN" "PRK" "IND" "MMR" "NIC"
## Taken from Polaris 2025 estimates:
Polaris_proptreatdiag <- Polaris2025[Polaris2025$ISO %in% missing.isos, colnames(Polaris2025) %in% c("Country_name","ISO","Diagnosed","Annual Treated")]

names(Polaris_proptreatdiag)[names(Polaris_proptreatdiag) == 'Annual Treated'] <- 'Polaris_treatment_coverage'
Polaris_proptreatdiag$Year <- 2025 ## All 2025 estimates
Polaris_proptreatdiag$OnTreatOfDiagnosed <- Polaris_proptreatdiag$Polaris_treatment_coverage/Polaris_proptreatdiag$Diagnosed
Polaris_proptreatdiag$OnTreatOfDiagnosed[is.na(Polaris_proptreatdiag$OnTreatOfDiagnosed)] <- 0
Polaris_proptreatdiag <- Polaris_proptreatdiag[,c(2,1,5,6)]
Prop_ontreat_of_diagnosed <- rbind(Prop_ontreat_of_diagnosed,Polaris_proptreatdiag)
missing.isos <- setdiff(list.of.isos,Prop_ontreat_of_diagnosed$ISO)
dim(Prop_ontreat_of_diagnosed)

treatment_cascade <- merge(Prop_diagnosed, Prop_ontreat_of_diagnosed, by=c("ISO","Country_name"))
dim(treatment_cascade)
head(treatment_cascade)
#setdiff(list.of.isos,treatment_cascade$ISO)
#table(treatment_cascade$Year.x,treatment_cascade$Year.y)
names(treatment_cascade)[names(treatment_cascade) == 'Year.x'] <- 'Year_diag'
names(treatment_cascade)[names(treatment_cascade) == 'Year.y'] <- 'Year_treatifdiag'

treatment_cascade$treatment_coverage <- treatment_cascade$Diagnosed * treatment_cascade$OnTreatOfDiagnosed 


## Current data outliers:
## There are a number of countries with diagnosis but no treatment - keep those as-is in status quo.
## Turkmenistan (44%), Zimbabwe (29%), Kenya (30%) and Madagascar (100%) have high proportions of diagnosed on treatment.
## Turkmenistan has 6192 on treatment out of 13962 diagnosed with chronic HepB in 2022. 6192/13962 = 44%
## Zimbabwe has 2164 on treatment out of 7343 diagnosed with chronic HepB in 2022. 2164/7343 = 29.4%
## Kenya 1714 on treatment out of 5804 diagnosed with chronic HepB in 2022. 1714/5804 = 29.5%
## Madagascar 5989 on treatment out of 3191 diagnosed in 2022 >> 100% so Madagascar maybe a data glitch?
## Polaris puts annual diagnosed <1% and annual treated <1% so not helpful.
## Madagascar is low-income SSA/ East SSA. Low-income have 7% diagnosis, 1% on treatment (so 1/7=0.14 on treatment). Africa has 2% treatment,7% diagnosed. AFRO has 2% treatment, 6% diagnosed. East Africa has <1% treated, 8% diagnosed. East Africa is probably most geographically relevant and similar to low-income, so take 1/8=12.5%
##write.csv(treatment_cascade,paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatitis B/Presentations/15May2026/treatment_cascade.csv"),row.names = FALSE)
treatment_cascade$OnTreatOfDiagnosed[treatment_cascade$ISO %in% "MDG"] <- 0.125
write.csv(treatment_cascade,paste0(basefolder,"Documents/Hepatitis_B/icl-hbv/resources/treatment_cascade_GHO.csv"),row.names = FALSE)


Polaris_raw_data <- read_excel(paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatitis B/Data/Polaris/Polaris_data_from2023paperLGastHep.xlsx"),sheet="Data")
Polaris_data <- Polaris_raw_data[,colnames(Polaris_raw_data) %in% c("ISO","Diagnosed","treatment_coverage","Diagnosed_Polaris2025","AnnualTreated_Polaris2025")]

## Compare Polaris 2023 paper and GHO:
temp <- merge(treatment_cascade,Polaris_data,by="ISO")

## 
# ISO Diagnosed.x Diagnosed.y
# 12  BLZ        0.44        0.09
# 16  CHN        0.68        0.24 *** This one ***
# 32  FSM        0.17        0.06
# 33  GEO        0.26        0.23
# 51  KIR        0.07        0.01
# 59  MHL        0.70        0.57
# 63  MNG        0.46        0.10
# 78  RWA        0.67        0.41
# 84  SOM        0.00        0.07
# 96  TON        0.34        0.28
# 98  TUV        0.23        0.05
temp[which(abs(temp$Diagnosed.x-temp$Diagnosed.y)>0.02),colnames(temp) %in% c("ISO","Diagnosed.x","Diagnosed.y")]

## treatment_coverage.x is treatment of sAg+; treatment_coverage.y is treatemnt of eligible; so expect treatment_coverage.y to be about 3x treatment_coverage.x. Roughly speaking Rwanda (22% of eligibles vs 4% of sAg+) is slightly off.
## The only major difference is china (30% of sAg+ is >>>> 15% of eligibles)
temp[which(abs(temp$treatment_coverage.x-temp$treatment_coverage.y)>0.05),colnames(temp) %in% c("ISO","treatment_coverage.x","treatment_coverage.y")]

## Given the large difference for China (and given Polaris's 2025 dataset - the excel above - gives China at 68% diagnosis and 30% treatemnt I will use that.)
## Comparing the 2023 and 2025 Polaris data, the differences between 2023 Polaris and GHO seem to be smaller/gone by 2025 Polaris.
## So use 2025 Polaris.

plot(Polaris_data$treatment_coverage,Polaris_data$AnnualTreated_Polaris2025)
setdiff(list.of.isos,Prop_ontreat_of_diagnosed$ISO)

## Polaris 2025 diagnosed data seems more up-to-date than Polaris 2023 paper.
## However, annual_treated_Polaris2025 seems to be smaller than 2023 value - it looks to me like the denominators are different (the 2025 probably uses all chronic sAg+).
## So use 2023 treatment_coverage. *EXCEPT FOR CHINA* where the treatment/diagnosis is 15% in 2023 Polaris but treatment/sAg+ is 30% in GHO and 2025 Polaris.

treatment_data_final <- Polaris_data[,colnames(Polaris_data) %in% c("ISO","Diagnosed_Polaris2025","treatment_coverage")]
treatment_data_final$treatment_coverage[treatment_data_final$ISO %in% "CHN"] <- Polaris_data$AnnualTreated_Polaris2025[Polaris_data$ISO %in% "CHN"]

colnames(treatment_data_final)[which(names(treatment_data_final) == "Diagnosed_Polaris2025")] <- "prop_diagnosed"
colnames(treatment_data_final)[which(names(treatment_data_final) == "treatment_coverage")] <- "treatment_coverage_of_eligible"
## Make order: ISO, prop_diagnosed, treatment_coverage_of_eligible:
treatment_data_final <- treatment_data_final[,c(1,3,2)]
head(treatment_data_final)

## Note that some final fixes are needed:
treatment_data_final[treatment_data_final$treatment_coverage_of_eligible > treatment_data_final$prop_diagnosed,]
# ISO   prop_diagnosed treatment_coverage_of_eligible
# <chr>          <dbl>                          <dbl>
#   1 AFG             0.05                          0.09 
# 2 EGY             0.06                          0.1  
# 3 KEN             0                             0.004
# 4 MWI             0.01                          0.03 
# 5 ZWE             0                             0.006
## AFG - 2025 dataset gives diag 5% treat <1%, so take treatment_coverage_of_eligible = 0.5%.
## EGY - 2025 dataset gives diag 5% treat 3%, so take treatment_coverage_of_eligible = 3%.
## KEN - 2025 dataset gives diag and treat both <1%, so take prop_diagnosed = 1% (so % of diagnosed on treat is not 100%).
## MWI - 2025 dataset gives diag 1% treat <1%, so take treatment_coverage_of_eligible = 0.5%.
## ZWE - 2025 dataset gives diag and treat both <1%, so take prop_diagnosed = 1% (so % of diagnosed on treat is not 100%).

treatment_data_final$treatment_coverage_of_eligible[treatment_data_final$ISO %in% "AFG"] <- 0.005
treatment_data_final$treatment_coverage_of_eligible[treatment_data_final$ISO %in% "EGY"] <- 0.03
treatment_data_final$prop_diagnosed[treatment_data_final$ISO %in% "KEN"] <- 0.01
treatment_data_final$treatment_coverage_of_eligible[treatment_data_final$ISO %in% "MWI"] <- 0.005
treatment_data_final$prop_diagnosed[treatment_data_final$ISO %in% "ZWE"] <- 0.01

## Double-check:
treatment_data_final[treatment_data_final$treatment_coverage_of_eligible > treatment_data_final$prop_diagnosed,]

write.csv(treatment_data_final,paste0(basefolder,"Documents/Hepatitis_B/icl-hbv/resources/treatment_cascade.csv"),row.names = FALSE)


# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("AZE","Azerbaijan",2025,5)
# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("BGD","Bangladesh",2025,1)
# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("CHN","China",2025,68)
# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("PRK","Korea, Dem. People's Rep.",2025,11)
# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("IND","India",2025,3)
# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("KIR","Kiribati",2025,7)
# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("MHL","Marshall Islands",2025,70)
# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("FSM","Micronesia, Fed. Sts.",2025,17)
# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("MMR","Myanmar",2025,0)
# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("WSM","Samoa",2025,3)
# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("TON","Tonga",2025,34)
# Prop_diagnosed[nrow(Prop_diagnosed) + 1,] = c("TUV","Tuvalu",2025,23)







ANC_coverage_raw <- read_excel(paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatitis B/Data/ANC_coverage/WB_ANC_coverage_to2025.xlsx"),sheet="Data")
ANC.colnames <- colnames(ANC_coverage_raw)
i.ISO <- which(ANC.colnames=="Country Code")
ANC_coverage_processed <- data.frame(ISO=character(), year_of_last_ANC_coverage_data=integer(),ANC_coverage=numeric())
for(i in 1:nrow(ANC_coverage_raw))
{
  i.last <- last(which(ANC_coverage_raw[i,]!=""))
  print(ANC.colnames[i.last])
  ANC_coverage_processed[nrow(ANC_coverage_processed) + 1,] = c(ANC_coverage_raw[i,i.ISO], as.integer(ANC.colnames[i.last]),as.numeric(ANC_coverage_raw[i,i.last]))
}


if("NIC" %in% list.of.isos){
  ANC_coverage_processed <- ANC_coverage_processed[ANC_coverage_processed$ISO %in% list.of.isos,]
}else{
  ANC_coverage_processed <- ANC_coverage_processed[ANC_coverage_processed$ISO %in% c(list.of.isos,"NIC"),]
}

dim(ANC_coverage_processed)
write.csv(ANC_coverage_processed,paste0(basefolder,"Documents/Hepatitis_B/icl-hbv/resources/ANC_coverage_data.csv"), row.names = FALSE)



## Compare HepB3 and ANC most recent

ANC.HepB3 <- merge(ANC_coverage_processed,hepB3.data[hepB3.data$Year %in% 2024,], by="ISO")

pdf(paste0(graph.dir,"WUENIC/HepB3_ANC_coverage.pdf"),width=6,height=6)
ggplot(ANC.HepB3, aes(x=ANC_coverage, y=WUENIC)) +
  geom_point(aes(color=WHO_region)) +
  geom_abline(slope=1,intercept=0) +
  geom_text(aes(label=ifelse(WUENIC<70,as.character(ISO),'')),hjust=0,vjust=0) + 
  xlab("WB ANC coverage, most recent year (%)") +
  ylab("WUENIC Hep B3 coverage, 2024 (%)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank())
dev.off()
