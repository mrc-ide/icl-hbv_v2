rm(list=ls())
library(readxl)

a <- read.csv("C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/src/raw_params/ListOfISOs.txt",sep=" ",header=FALSE)
ListOfISOsModelled <- a$V2
rm(a)


setwd("/Users/mpickles/OneDrive - Imperial College London/Dropbox_copy/Hepatits B/")

gbd <- read_excel("ListOfCountryRegions.xlsx", sheet = "GBD_processed")
who.df <- read_excel("ListOfCountryRegions.xlsx", sheet = "WHO_regions")
wb <- read_excel("ListOfCountryRegions.xlsx", sheet = "WorldBankIncomeGroups")

## Check for missing regions:
gbd[is.na(gbd$GBD_region),]
gbd[is.na(gbd$GBD_region_aggregated),]
who.df[is.na(who.df$ParentLocation),]
wb[is.na(wb$WB_Region),]


## Check for duplicates:
who.df$Location[duplicated(who.df$Location)]
gbd$Country_WHOformat[duplicated(gbd$Country_WHOformat)]
wb$Code[duplicated(wb$Code)] ## WB ISO
who.df$SpatialDimValueCode[duplicated(who.df$SpatialDimValueCode)] ## WHO ISO

setdiff(who.df$SpatialDimValueCode,wb$Code)
## character(0)
setdiff(wb$Code,who.df$SpatialDimValueCode)
##[1] "ASM" "ABW" "BMU" "VGB" "CYM" "CHI" "CUW" "FRO" "PYF" "GIB" "GRL" "GUM" "HKG" "IMN" "XKX" "LIE" "MAC" "MCO" "NCL"
##[20] "MNP" "SMR" "SXM" "MAF" "TWN" "TCA" "VIR"
# wb$Country[wb$Code %in% setdiff(wb$Code,who.df$SpatialDimValueCode)]
# [1] "American Samoa"            "Aruba"                     "Bermuda"                   "British Virgin Islands"   
# [5] "Cayman Islands"            "Channel Islands"           "Curaçao"                   "Faroe Islands"            
# [9] "French Polynesia"          "Gibraltar"                 "Greenland"                 "Guam"                     
# [13] "Hong Kong SAR, China"      "Isle of Man"               "Kosovo"                    "Liechtenstein"            
# [17] "Macao SAR, China"          "Monaco"                    "New Caledonia"             "Northern Mariana Islands" 
# [21] "San Marino"                "Sint Maarten (Dutch part)" "St. Martin (French part)"  "Taiwan, China"            
# [25] "Turks and Caicos Islands"  "Virgin Islands (U.S.)"   


setdiff(gbd$Country_WHOformat,who.df$Location)
## "Greenland" "Bermuda"

setdiff(who.df$Location,gbd$Country_WHOformat)
## character(0)


################################################################################################################


## Tidy up data frames and merge:
gbd$Country <- gbd$Country_WHOformat
gbd <- gbd[,colnames(gbd) %in% c("Country","GBD_region","GBD_region_aggregated")]

who.df$Country <- who.df$Location
who.df$ISO <- who.df$SpatialDimValueCode
who.df$WHO_region <- who.df$ParentLocation
who.df <- who.df[,colnames(who.df) %in% c("Country","ISO","WHO_region")]


wb$ISO <- wb$Code
wb$WB_income_group <- wb$`Income group`
wb <- wb[,colnames(wb) %in% c("ISO","WB_Region","WB_income_group ")]

CountryClassifierDF <- merge(who.df, wb, by=c("ISO"))

CountryClassifierDF <- merge(CountryClassifierDF, gbd, by=c("Country"))

setdiff(CountryClassifierDF$ISO,ListOfISOsModelled)
## 82 countries

setdiff(ListOfISOsModelled, CountryClassifierDF$ISO)
## character(0)


CountryClassifierDF <- CountryClassifierDF[CountryClassifierDF$ISO %in% ListOfISOsModelled,]
write.csv(CountryClassifierDF, file="C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/resources/CountryClassifier.csv", row.names=FALSE)

