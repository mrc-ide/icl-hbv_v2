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


basefolder <- "C:/Users/mpickles/"
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




##write.csv(bd.coverage.2024,paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatits B/Presentations/15May2026/BDdatatemp.csv"), row.names = FALSE)
write.csv(bd.coverage.2024,paste0(basefolder,"Documents/Hepatitis_B/icl-hbv/resources/wuenic2024data/BDdata_infacility_rural.csv"), row.names = FALSE)



#######################################################################################
## Now extract GHO treatment data:


Prop_ontreat_of_diagnosed <- read_excel(paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatits B/Data/GHO_treatment/Prop_ontreat_of_diagnosed.xlsx"),sheet="Prop_ontreat_of_diagnosed")

Prop_diagnosed <- read_excel(paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatits B/Data/GHO_treatment/Prop_diagnosed.xlsx"),sheet="Prop_diagnosed")

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
Polaris2025 <- read_excel(paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatits B/Data/Polaris/Polaris Database Query – CDA Foundation.xlsx"),sheet="Cleaned")
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


##write.csv(treatment_cascade,paste0(basefolder,"OneDrive - Imperial College London/Dropbox_copy/Hepatits B/Presentations/15May2026/treatment_cascade.csv"),row.names = FALSE)
write.csv(treatment_cascade,paste0(basefolder,"Documents/Hepatitis_B/icl-hbv/resources/treatment_cascade.csv"),row.names = FALSE)


setdiff(list.of.isos,Prop_ontreat_of_diagnosed$SpatialDimValueCode)


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










