rm(list=ls())
library(stringi)
library(readxl)
library(tibble)
##library(ggpubr)
##library("gridExtra")

## Script takes GHO data on in-facility births from
## https://www.who.int/data/gho/data/indicators/indicator-details/GHO/institutional-births-(-)
## Last updated 2024-08-06
## Downloaded 6 Feb 2026
## Births delivered in a health facility (facility/institutional births), proportion (%)

GHOpercentdeliveredinhealthfacility_alldata <- read_excel("C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/resources/wuenic2024data/GHO_prop_births_deliveredinhealthfaility_download5Feb2026.xlsx",sheet="Data")
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

GHOpercentdeliveredinhealthfacility_rawdata$prop_delivered_in_facility <- GHOpercentdeliveredinhealthfacility_rawdata$Percent_deliveredInHealthFacility/100

GHOpercentdeliveredinhealthfacility_rawdata <- GHOpercentdeliveredinhealthfacility_rawdata[,colnames(GHOpercentdeliveredinhealthfacility_rawdata) %in% c("ISO","prop_delivered_in_facility")]
write.csv(x=GHOpercentdeliveredinhealthfacility_rawdata,file="C:/Users/mpickles/Documents/Hepatitis_B/icl-hbv/resources/GHO_in_facility_births_summary.csv",row.names=FALSE, quote = FALSE)
