#################
# EMR catchment and selection bias
# Citation: Goldstein ND, Kahal D, Testa K, Burstyn I. Inverse probability weighting for selection bias in a Delaware community health center electronic medical record study of community deprivation and hepatitis C prevalence. Manuscript in preparation.
# 4/3/19 -- Neal Goldstein
#################

### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable
library(tidycensus) #retrieve ACS data, note if error installing on MacOS see: https://github.com/r-quantities/units/issues/1
library(maptools) #maptools provides the spatial recognition
library(lme4) #mixed effects modeling, dichotomous
library(blme) #bayesian nonlinear mixed effects (bglmer)
library(RColorBrewer)
library(boot) #bootstrapping for IPW CIs


### READ DATA ###

#WFH patient data: exports from FQHC
load("Retrospective_redcap_NE.Rdata")
retro_NE = retro
load("Retrospective_redcap_4th.Rdata")
retro_4th = retro
rm(retro)

#zip codes in DE: Delaware from U.S. Census Tiger website, New Castle County from DE open data portal https://data-nccde.opendata.arcgis.com/datasets/download-zipcodes-shapefile
de_zcta = readShapePoly("tl_2010_10_zcta510/tl_2010_10_zcta510")
ncc_zcta = readShapePoly("Zipcodes_SHP/Zipcodes")

#zctas for crosswalk: from UDS mapper https://www.udsmapper.org/zcta-crosswalk.cfm
zcta = read.csv("zip_to_zcta_2019.csv", as.is=T, stringsAsFactors=F)

#retrieve census variables of interest: using tidycensus but could manually obtain from FactFinder
census_api_key("PASTE API KEY HERE")
population = get_acs(geography="zcta", table="B01003", year=2017, output="wide")
income = get_acs(geography="zcta", table="B19013", year=2017, output="wide")
fertility = get_acs(geography="zcta", variables="B13016_001", year=2017, output="wide")
race = get_acs(geography="zcta", table="B02001", year=2017, output="wide")
insurance_employment = get_acs(geography="zcta", table="B27011", year=2017, output="wide")

#census-based deprivation index: details provided here https://towardsdatascience.com/a-census-based-deprivation-index-using-r-7aa738da697c, https://www.ncbi.nlm.nih.gov/pubmed/17031568
deprivation = get_acs(geography="zcta", variables=c("B17001_002", "B17001_001", "B06009_002" , "B06009_001","B09008_011","B09008_001","B08124_002", "B08124_001", "B25014_005","B25014_006",  "B25014_007","B25014_011", "B25014_012", "B25014_013","B25014_001", "B19058_002", "B19058_001","C23002C_021", "C23002D_008","C23002C_017", "C23002D_003","B19001_002", "B19001_003", "B19001_004","B19001_005", "B19001_006", "B19001_001"), output="wide", year=2017)

#bus stops in NCC Delaware: from DE open data portal https://firstmap-delaware.opendata.arcgis.com/datasets/delaware-bus-stops/data
busstops = read.csv("Delaware_Bus_Stops.csv", as.is=T, stringsAsFactors=F)
busstops = busstops[busstops$COUNTY=="New Castle County", ]


### CODE and LINK DATA ###

#based on end point of 2018-10-31 for NE and 2019-06-30 for 4th
retro_NE$days_since_last_visit = as.numeric(as.Date("2018-10-31") - retro_NE$last_visit_date)
retro_4th$days_since_last_visit = as.numeric(as.Date("2019-06-30") - retro_4th$last_visit_date)

#white vs non-white
retro_NE$race_2cat = ifelse(retro_NE$race==0, 0, 1)
retro_4th$race_2cat = ifelse(retro_4th$race==0, 0, 1)

#public and no insurance vs. insured
retro_NE$insurance_2cat = ifelse(retro_NE$insurance==1, 1, 0)
retro_4th$insurance_2cat = ifelse(retro_4th$insurance==1, 1, 0)

#recode census variables
population$census_population = population$B01003_001E
income$census_income = income$B19013_001E
fertility$census_fertility = fertility$B13016_001E
race$census_white_percent = race$B02001_002E/race$B02001_001E*100
race$census_black_percent = race$B02001_003E/race$B02001_001E*100
race$census_other_percent = 100 - (race$census_white_percent + race$census_black_percent)
insurance_employment$census_employed_percent = insurance_employment$B27011_003E/insurance_employment$B27011_002E*100
insurance_employment$census_insured_percent = (insurance_employment$B27011_004E+insurance_employment$B27011_009E+insurance_employment$B27011_014E)/insurance_employment$B27011_001E*100
insurance_employment$census_insured_private_percent = (insurance_employment$B27011_005E+insurance_employment$B27011_010E+insurance_employment$B27011_015E)/insurance_employment$B27011_001E*100
insurance_employment$census_insured_govt_percent = (insurance_employment$B27011_006E+insurance_employment$B27011_011E+insurance_employment$B27011_016E)/insurance_employment$B27011_001E*100

#create deprivation index: https://towardsdatascience.com/a-census-based-deprivation-index-using-r-7aa738da697c, https://www.ncbi.nlm.nih.gov/pubmed/17031568
deprivation$pct_poverty = deprivation$B17001_002E / deprivation$B17001_001E
deprivation$pct_noHS = deprivation$B06009_002E / deprivation$B06009_001E
deprivation$pct_FHH = deprivation$B09008_011E / deprivation$B09008_001E
deprivation$pct_mgmt = deprivation$B08124_002E / deprivation$B08124_001E 
deprivation$pct_crowd = (deprivation$B25014_005E + deprivation$B25014_006E + deprivation$B25014_007E + deprivation$B25014_011E + deprivation$B25014_012E + deprivation$B25014_013E) / deprivation$B25014_001E
deprivation$pct_pubassist = deprivation$B19058_002E / deprivation$B19058_001E
deprivation$pct_unempl = (deprivation$C23002C_021E + deprivation$C23002D_008E) / (deprivation$C23002C_017E + deprivation$C23002D_003E)
deprivation$pct_under30K = ((deprivation$B19001_002E + deprivation$B19001_003E + deprivation$B19001_004E + deprivation$B19001_005E + deprivation$B19001_006E) / deprivation$B19001_001E)
deprivation_matrix = as.matrix(deprivation[, c("pct_poverty","pct_noHS","pct_FHH","pct_mgmt","pct_crowd","pct_pubassist", "pct_unempl","pct_under30K")])
deprivation_matrix[is.nan(deprivation_matrix)] = 0
deprivation$census_NDI = principal(deprivation_matrix,nfactors = 1)$scores 

#join individual data from both sites
retro_NE$clinic = "NE"
retro_4th$clinic = "4th"
retro_4th$redcap_id = retro_4th$redcap_id + 10000
retro = rbind(retro_NE, retro_4th)

#clean zip codes
retro$zip = as.character(retro$zip)
retro$zip[nchar(retro$zip)<5] = paste("0",retro$zip[nchar(retro$zip)<5],sep="")
retro$zip[nchar(retro$zip)>5] = substr(retro$zip[nchar(retro$zip)>5],1,5)

#crosswalk to zctas
retro$zcta = NA
for (i in 1:nrow(retro)) {
  z = zcta$ZCTA[which(zcta$ZIP_CODE==as.numeric(retro$zip[i]))]
  retro$zcta[i] = ifelse(length(z)>0, z, NA)
}
rm(i,z)

#merge census data to individual data
retro = merge(retro, population[,c("GEOID","census_population")], by.x="zcta", by.y="GEOID", all.x=T, all.y=F)
retro = merge(retro, income[,c("GEOID","census_income")], by.x="zcta", by.y="GEOID", all.x=T, all.y=F)
retro = merge(retro, fertility[,c("GEOID","census_fertility")], by.x="zcta", by.y="GEOID", all.x=T, all.y=F)
retro = merge(retro, race[,c("GEOID","census_white_percent","census_black_percent","census_other_percent")], by.x="zcta", by.y="GEOID", all.x=T, all.y=F)
retro = merge(retro, insurance_employment[,c("GEOID","census_employed_percent","census_insured_percent","census_insured_private_percent","census_insured_govt_percent")], by.x="zcta", by.y="GEOID", all.x=T, all.y=F)
retro = merge(retro, deprivation[,c("GEOID","census_NDI")], by.x="zcta", by.y="GEOID", all.x=T, all.y=F)

#map bus stop zip codes
busstops$zcta = NA
for (i in 1:nrow(busstops)) {
  
  #create a spatial point object from current coordinate
  spatial_pt = SpatialPoints(data.frame(x = as.numeric(busstops$X[i]), y = as.numeric(busstops$Y[i])))
  
  #map the point to the index
  map_index = over(spatial_pt, de_zcta)
  
  busstops$zcta[i] = as.character(map_index$ZCTA5CE10)
}
rm(i,spatial_pt,map_index)

#join to retro count of bus stops per zcta
retro$bus_stops = NA
unique_zcta = unique(retro$zcta)
for (i in 1:length(unique_zcta)) {
  n_stops = sum(busstops$zcta==unique_zcta[i])
  retro$bus_stops[retro$zcta==unique_zcta[i]] = ifelse(length(n_stops)>0, n_stops, NA)
}
rm(i,unique_zcta,n_stops)

#create a Wilmington DE dataset, zip codes from: https://tools.usps.com/zip-code-lookup.htm?bycitystate
de_wilmington = data.frame("zip"=c("19801","19802","19803","19804","19805","19806","19807","19808","19809","19810"),"zcta"=NA,"n_patients_4th"=NA,"n_patients_NE"=NA,"n_patients_total"=NA,"population"=NA, stringsAsFactors=F)
for (i in 1:nrow(de_wilmington)) {
  
  #zcta crosswalk
  z = zcta$ZCTA[which(zcta$ZIP_CODE==as.numeric(de_wilmington$zip[i]))]
  de_wilmington$zcta[i] = ifelse(length(z)>0, z, NA)
  
  #number of patients
  de_wilmington$n_patients_4th[i] = sum(retro$zcta==de_wilmington$zcta[i] & retro$clinic=="4th", na.rm=T)
  de_wilmington$n_patients_NE[i] = sum(retro$zcta==de_wilmington$zcta[i] & retro$clinic=="NE", na.rm=T)
  de_wilmington$n_patients_total[i] = sum(retro$zcta==de_wilmington$zcta[i], na.rm=T)
  
  #census population
  de_wilmington$population[i] = population$census_population[population$GEOID==de_wilmington$zcta[i]]
}
rm(i,z)

#create an NCC DE dataset, zip codes from New Castle County shapefile
de_ncc = data.frame("zip"=unique(ncc_zcta$ZIPCODE[order(ncc_zcta$ZIPCODE)]),"zcta"=NA,"n_patients_4th"=NA,"n_patients_NE"=NA,"n_patients_total"=NA,"population"=NA,"n_stops"=NA, stringsAsFactors=F)
for (i in 1:nrow(de_ncc)) {
  
  #zcta crosswalk
  z = zcta$ZCTA[which(zcta$ZIP_CODE==as.numeric(de_ncc$zip[i]))]
  de_ncc$zcta[i] = ifelse(length(z)>0, z, NA)
  
  #number of patients
  de_ncc$n_patients_4th[i] = sum(retro$zcta==de_ncc$zcta[i] & retro$clinic=="4th", na.rm=T)
  de_ncc$n_patients_NE[i] = sum(retro$zcta==de_ncc$zcta[i] & retro$clinic=="NE", na.rm=T)
  de_ncc$n_patients_total[i] = sum(retro$zcta==de_ncc$zcta[i], na.rm=T)
  
  #number of bus stops
  n_stops = sum(busstops$zcta==de_ncc$zcta[i])
  de_ncc$n_stops[i] = ifelse(length(n_stops)>0, n_stops, NA)
  
  #census population
  de_ncc$population[i] = population$census_population[population$GEOID==de_ncc$zcta[i]]
}
rm(i,z,n_stops)


### SAVE ANALYSIS DATASET ###

#rm(deprivation_matrix,deprivation,fertility,income,insurance_employment,population,race,retro_NE,retro_4th,zcta,busstops)
save.image("Retrospective_analysis.Rdata")
write.csv(de_ncc, "NCC_data.csv",na="",row.names=F)
write.csv(retro, "Patient_data.csv",na="",row.names=F)


### LOAD DATA ###

load("/Retrospective_analysis.Rdata")

#subset to NCC patients for analyses
retro_ncc = retro[retro$zcta %in% de_ncc$zcta, ]


### ADD QGIS CALCULATIONS ###

#distance:
#since this is based on zcta, first random points corresponding to each patient are generated (which acts as their address)
#then this is randomly joined to the individual patient, therefore each patient is randomly located within a given zcta

distances_4th = read.csv("Patients_4th_distance.csv", as.is=T, stringsAsFactors=F)
distances_ne = read.csv("Patients_ne_distance.csv", as.is=T, stringsAsFactors=F)

retro_ncc$distance_to_clinic = NA
for (i in 1:nrow(retro_ncc)) {
  if (retro_ncc$clinic[i]=="4th") {
    rand_distances = distances_4th$HubDist[distances_4th$zcta==as.numeric(retro_ncc$zcta[i])]
    retro_ncc$distance_to_clinic[i] = ifelse(length(rand_distances)>0, sample(rand_distances,1), NA)
  } else if (retro_ncc$clinic[i]=="NE") {
    rand_distances = distances_ne$HubDist[distances_ne$zcta==as.numeric(retro_ncc$zcta[i])]
    retro_ncc$distance_to_clinic[i] = ifelse(length(rand_distances)>0, sample(rand_distances,1), NA)
  }
}
rm(i,rand_distances)

#busstop density per sq mi
busstops_density = read.csv("busstop_density.csv", as.is=T, stringsAsFactors=F)

#omit the zcta w/ no patients
busstops_density = busstops_density[complete.cases(busstops_density), ]

#note: one zcta appears twice due to geography (19732), just take first value
retro_ncc$bus_stops_density = NA
for (i in 1:nrow(retro_ncc)) {
  retro_ncc$bus_stops_density[i] = round(busstops_density$n_stops_sqmi[busstops_density$zcta==retro_ncc$zcta[i]],1)[1]
}


### CATCHMENT AREA (CA) ANALYSIS ###

#distance-based CA (QGIS analysis)
retro_ncc$CA_distance_10mi = ifelse(retro_ncc$distance_to_clinic<=10, 1, 0)

#cumulative distance proportion CA (QGIS analysis)
retro_ncc$CA_distance_75per[retro_ncc$clinic=="4th"] = ifelse(retro_ncc$distance_to_clinic[retro_ncc$clinic=="4th"]<=quantile(retro_ncc$distance_to_clinic[retro_ncc$clinic=="4th"], probs=0.75, na.rm=T), 1, 0)
retro_ncc$CA_distance_80per[retro_ncc$clinic=="4th"] = ifelse(retro_ncc$distance_to_clinic[retro_ncc$clinic=="4th"]<=quantile(retro_ncc$distance_to_clinic[retro_ncc$clinic=="4th"], probs=0.80, na.rm=T), 1, 0)
retro_ncc$CA_distance_90per[retro_ncc$clinic=="4th"] = ifelse(retro_ncc$distance_to_clinic[retro_ncc$clinic=="4th"]<=quantile(retro_ncc$distance_to_clinic[retro_ncc$clinic=="4th"], probs=0.90, na.rm=T), 1, 0)
retro_ncc$CA_distance_75per[retro_ncc$clinic=="NE"] = ifelse(retro_ncc$distance_to_clinic[retro_ncc$clinic=="NE"]<=quantile(retro_ncc$distance_to_clinic[retro_ncc$clinic=="NE"], probs=0.75, na.rm=T), 1, 0)
retro_ncc$CA_distance_80per[retro_ncc$clinic=="NE"] = ifelse(retro_ncc$distance_to_clinic[retro_ncc$clinic=="NE"]<=quantile(retro_ncc$distance_to_clinic[retro_ncc$clinic=="NE"], probs=0.80, na.rm=T), 1, 0)
retro_ncc$CA_distance_90per[retro_ncc$clinic=="NE"] = ifelse(retro_ncc$distance_to_clinic[retro_ncc$clinic=="NE"]<=quantile(retro_ncc$distance_to_clinic[retro_ncc$clinic=="NE"], probs=0.90, na.rm=T), 1, 0)

max(retro_ncc$distance_to_clinic[retro_ncc$CA_distance_75per==1])
max(retro_ncc$distance_to_clinic[retro_ncc$CA_distance_80per==1])
max(retro_ncc$distance_to_clinic[retro_ncc$CA_distance_90per==1])

#Poisson model of spatial scan statistic (SaTScan analysis)
cluster_zip = c("19801","19802","19806","19805","19809","19804","19710","19706","19731","19733","19720")
retro_ncc$CA_poisson = ifelse(retro_ncc$zcta %in% cluster_zip, 1, 0)

#check if CA is related to exposure and outcome; evidence of possible selection bias
summary(lm(census_NDI ~ CA_distance_90per, data=retro_ncc))
summary(glm(hepc ~ CA_distance_90per, data=retro_ncc, family=binomial()))

#predictors of CA
#summary(glm(CA_distance_10mi ~ last_visit_age, family=binomial(), data=retro_ncc))
summary(glm(CA_distance_10mi ~ as.factor(sex), family=binomial(), data=retro_ncc))
#summary(glm(CA_distance_10mi ~ as.factor(race_2cat), family=binomial(), data=retro_ncc))
#summary(glm(CA_distance_10mi ~ as.factor(ethnicity), family=binomial(), data=retro_ncc))
summary(glm(CA_distance_10mi ~ as.factor(insurance_2cat), family=binomial(), data=retro_ncc))
summary(glm(CA_distance_10mi ~ n_visits, family=binomial(), data=retro_ncc))
#summary(glm(CA_distance_10mi ~ bus_stops, family=binomial(), data=retro_ncc))

summary(glm(CA_distance_80per ~ last_visit_age, family=binomial(), data=retro_ncc))
#summary(glm(CA_distance_80per ~ as.factor(sex), family=binomial(), data=retro_ncc))
summary(glm(CA_distance_80per ~ as.factor(race_2cat), family=binomial(), data=retro_ncc))
summary(glm(CA_distance_80per ~ as.factor(ethnicity), family=binomial(), data=retro_ncc))
summary(glm(CA_distance_80per ~ as.factor(insurance_2cat), family=binomial(), data=retro_ncc))
summary(glm(CA_distance_80per ~ n_visits, family=binomial(), data=retro_ncc))
summary(glm(CA_distance_80per ~ bus_stops_density, family=binomial(), data=retro_ncc))

#tertiles for modeling
#retro_ncc$n_visits_quantile = cut(retro_ncc$n_visits, quantile(retro_ncc$n_visits, probs=c(0,0.33,0.67,1)), include.lowest=T)
#retro_ncc$bus_stops_quantile = cut(retro_ncc$bus_stops, quantile(retro_ncc$bus_stops, probs=c(0,0.33,0.67,1), na.rm=T), include.lowest=T)

#CA model using cumulative distance proportion
#model = glmer(CA_distance_10mi ~ (1 | zcta) + as.factor(sex) + as.factor(insurance_2cat) + as.factor(n_visits_quantile), family=binomial(), data=retro_ncc, control=glmerControl(optimizer="bobyqa"))
model = glmer(CA_distance_75per ~ (1 | zcta) + scale(last_visit_age) + as.factor(race_2cat) + as.factor(ethnicity) + as.factor(insurance_2cat) + scale(n_visits) + bus_stops_density, family=binomial(), data=retro_ncc, control=glmerControl(optimizer="bobyqa"))
model = glmer(CA_distance_80per ~ (1 | zcta) + scale(last_visit_age) + as.factor(race_2cat) + as.factor(ethnicity) + as.factor(insurance_2cat) + scale(n_visits) + bus_stops_density, family=binomial(), data=retro_ncc, control=glmerControl(optimizer="bobyqa"))
model = glmer(CA_distance_90per ~ (1 | zcta) + scale(last_visit_age) + as.factor(race_2cat) + as.factor(ethnicity) + as.factor(insurance_2cat) + scale(n_visits) + bus_stops_density, family=binomial(), data=retro_ncc, control=glmerControl(optimizer="bobyqa"))
summary(model)
exp(fixef(model))
exp(confint.merMod(model, method="Wald"))

#PS: probablity of being in catchment (ie propensity score)
#IPW: weighting for selection bias correction; see https://doi.org/10.1177%2F2167696815621645
#note: be sure to update catchment variables in these code depending on CA model
retro_ncc$PS[complete.cases(retro_ncc[,c("CA_distance_75per","zcta","race_2cat","ethnicity","insurance_2cat","n_visits","bus_stops")])] = predict(model, type="response")
retro_ncc$IPW = ifelse(retro_ncc$CA_distance_75per==1, 1/retro_ncc$PS, 1/(1-retro_ncc$PS))
retro_ncc$IPW_stabilized = ifelse(retro_ncc$CA_distance_75per==1, mean(retro_ncc$CA_distance_75per)/retro_ncc$PS, (1-mean(retro_ncc$CA_distance_75per))/(1-retro_ncc$PS))

#cap extreme weights by 2nd and 98th percentiles for sensitivity analysis
IPW_quantile = as.numeric(quantile(retro_ncc$IPW,c(.02,.98), na.rm=T))
retro_ncc$IPW_trimmed = ifelse(retro_ncc$IPW < IPW_quantile[1], IPW_quantile[1], retro_ncc$IPW)
retro_ncc$IPW_trimmed = ifelse(retro_ncc$IPW > IPW_quantile[2], IPW_quantile[2], retro_ncc$IPW)
IPW_quantile = as.numeric(quantile(retro_ncc$IPW_stabilized,c(.02,.98), na.rm=T))
retro_ncc$IPW_stabilized_trimmed = ifelse(retro_ncc$IPW_stabilized < IPW_quantile[1], IPW_quantile[1], retro_ncc$IPW_stabilized)
retro_ncc$IPW_stabilized_trimmed = ifelse(retro_ncc$IPW_stabilized > IPW_quantile[2], IPW_quantile[2], retro_ncc$IPW_stabilized)
rm(IPW_quantile)


### EXPLORATORY ANALYSES ###

#general characteristics
nrow(retro)
nrow(retro_ncc)
describe(retro_ncc$last_visit_age); IQR(retro_ncc$last_visit_age)
CrossTable(retro_ncc$sex)
CrossTable(retro_ncc$race_2cat)
CrossTable(retro_ncc$ethnicity)
CrossTable(retro_ncc$insurance_2cat)
CrossTable(retro_ncc$clinic)
describe(retro_ncc$n_visits); IQR(retro_ncc$n_visits)
describe(retro_ncc$days_since_last_visit); IQR(retro_ncc$days_since_last_visit)
CrossTable(retro_ncc$hepc)
describe(retro_ncc$distance_to_clinic); IQR(retro_ncc$distance_to_clinic)
describe(retro_ncc$bus_stops_density); IQR(retro_ncc$bus_stops_density)
describe(retro_ncc$census_NDI)

sum(retro_ncc$CA_distance_75per==1)
describe(retro_ncc$last_visit_age[retro_ncc$CA_distance_75per==1]); IQR(retro_ncc$last_visit_age[retro_ncc$CA_distance_75per==1])
CrossTable(retro_ncc$sex[retro_ncc$CA_distance_75per==1])
CrossTable(retro_ncc$race_2cat[retro_ncc$CA_distance_75per==1])
CrossTable(retro_ncc$ethnicity[retro_ncc$CA_distance_75per==1])
CrossTable(retro_ncc$insurance_2cat[retro_ncc$CA_distance_75per==1])
CrossTable(retro_ncc$clinic[retro_ncc$CA_distance_75per==1])
describe(retro_ncc$n_visits[retro_ncc$CA_distance_75per==1]); IQR(retro_ncc$n_visits[retro_ncc$CA_distance_75per==1])
describe(retro_ncc$days_since_last_visit[retro_ncc$CA_distance_75per==1]); IQR(retro_ncc$days_since_last_visit[retro_ncc$CA_distance_75per==1])
CrossTable(retro_ncc$hepc[retro_ncc$CA_distance_75per==1])
describe(retro_ncc$distance_to_clinic[retro_ncc$CA_distance_75per==1]); IQR(retro_ncc$distance_to_clinic[retro_ncc$CA_distance_75per==1])
describe(retro_ncc$bus_stops_density[retro_ncc$CA_distance_75per==1]); IQR(retro_ncc$bus_stops_density[retro_ncc$CA_distance_75per==1])
describe(retro_ncc$census_NDI[retro_ncc$CA_distance_75per==1])

sum(retro_ncc$CA_distance_80per==1)
describe(retro_ncc$last_visit_age[retro_ncc$CA_distance_80per==1]); IQR(retro_ncc$last_visit_age[retro_ncc$CA_distance_80per==1])
CrossTable(retro_ncc$sex[retro_ncc$CA_distance_80per==1])
CrossTable(retro_ncc$race_2cat[retro_ncc$CA_distance_80per==1])
CrossTable(retro_ncc$ethnicity[retro_ncc$CA_distance_80per==1])
CrossTable(retro_ncc$insurance_2cat[retro_ncc$CA_distance_80per==1])
CrossTable(retro_ncc$clinic[retro_ncc$CA_distance_80per==1])
describe(retro_ncc$n_visits[retro_ncc$CA_distance_80per==1]); IQR(retro_ncc$n_visits[retro_ncc$CA_distance_80per==1])
describe(retro_ncc$days_since_last_visit[retro_ncc$CA_distance_80per==1]); IQR(retro_ncc$days_since_last_visit[retro_ncc$CA_distance_80per==1])
CrossTable(retro_ncc$hepc[retro_ncc$CA_distance_80per==1])
describe(retro_ncc$distance_to_clinic[retro_ncc$CA_distance_80per==1]); IQR(retro_ncc$distance_to_clinic[retro_ncc$CA_distance_80per==1])
describe(retro_ncc$bus_stops_density[retro_ncc$CA_distance_80per==1]); IQR(retro_ncc$bus_stops_density[retro_ncc$CA_distance_80per==1])
describe(retro_ncc$census_NDI[retro_ncc$CA_distance_80per==1])

sum(retro_ncc$CA_distance_90per==1)
describe(retro_ncc$last_visit_age[retro_ncc$CA_distance_90per==1]); IQR(retro_ncc$last_visit_age[retro_ncc$CA_distance_90per==1])
CrossTable(retro_ncc$sex[retro_ncc$CA_distance_90per==1])
CrossTable(retro_ncc$race_2cat[retro_ncc$CA_distance_90per==1])
CrossTable(retro_ncc$ethnicity[retro_ncc$CA_distance_90per==1])
CrossTable(retro_ncc$insurance_2cat[retro_ncc$CA_distance_90per==1])
CrossTable(retro_ncc$clinic[retro_ncc$CA_distance_90per==1])
describe(retro_ncc$n_visits[retro_ncc$CA_distance_90per==1]); IQR(retro_ncc$n_visits[retro_ncc$CA_distance_90per==1])
describe(retro_ncc$days_since_last_visit[retro_ncc$CA_distance_90per==1]); IQR(retro_ncc$days_since_last_visit[retro_ncc$CA_distance_90per==1])
CrossTable(retro_ncc$hepc[retro_ncc$CA_distance_90per==1])
describe(retro_ncc$distance_to_clinic[retro_ncc$CA_distance_90per==1]); IQR(retro_ncc$distance_to_clinic[retro_ncc$CA_distance_90per==1])
describe(retro_ncc$bus_stops_density[retro_ncc$CA_distance_90per==1]); IQR(retro_ncc$bus_stops_density[retro_ncc$CA_distance_90per==1])
describe(retro_ncc$census_NDI[retro_ncc$CA_distance_90per==1])

#prevalence by distance thresholds
CrossTable(retro_ncc$hepc, retro_ncc$CA_distance_75per, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(retro_ncc$hepc, retro_ncc$CA_distance_80per, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(retro_ncc$hepc, retro_ncc$CA_distance_90per, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

#prevalence of HCV by PS and IPW
IPW_quantile = as.numeric(quantile(retro_ncc$IPW,c(.05,.95), na.rm=T))
sum(retro_ncc$hepc[retro_ncc$IPW <= IPW_quantile[1]], na.rm=T) / length(retro_ncc$hepc[retro_ncc$IPW <= IPW_quantile[1]]) * 100
sum(retro_ncc$hepc[retro_ncc$IPW >= IPW_quantile[2]], na.rm=T) / length(retro_ncc$hepc[retro_ncc$IPW >= IPW_quantile[2]]) * 100

uniqueIPW = unique(retro_ncc$IPW)
uniqueIPW = uniqueIPW[order(uniqueIPW)]
prev_below = NA #in catchment
prev_at_or_above = NA #out of catchment
for (i in 2:(length(uniqueIPW)-1)) {
  prev_below = c(prev_below, sum(retro_ncc$hepc[retro_ncc$IPW < uniqueIPW[i]], na.rm=T) / length(retro_ncc$hepc[retro_ncc$IPW < uniqueIPW[i]]) * 100)
  prev_at_or_above = c(prev_at_or_above, sum(retro_ncc$hepc[retro_ncc$IPW >= uniqueIPW[i]], na.rm=T) / length(retro_ncc$hepc[retro_ncc$IPW >= uniqueIPW[i]]) * 100)
}
rm(i)

describe(prev_below)
describe(prev_at_or_above)

plot(prev_below, type="l", ylim=c(0,2), ylab="HCV prevalence", xlab="IPW percentile (below/above)", xaxt="n", lwd=3, col=gray.colors(2)[1])
axis(side=1, at=seq(2,(length(uniqueIPW)+10), by=length(uniqueIPW)/5), labels=c("0/100","20/80","40/60","60/40","80/20","100/0"))
lines(rev(prev_at_or_above), type="l", lwd=3, col=gray.colors(2)[2])
legend("top", legend=c("Below","Above"), lwd=3, col=gray.colors(2), horiz=T)

model = glmer(hepc ~ (1 | zip) + scale(IPW), family=binomial(), data=retro_ncc)
summary(model)

#association of IPW with ADI and HCV
boxplot(retro_ncc$IPW ~ retro_ncc$hepc, xlab="HCV dx", ylab="IPW", main="Relationship between IPW and HCV")
plot(retro_ncc$IPW, retro_ncc$census_NDI, xlab="IPW", ylab="Area Deprivation Index", main="Relationship between IPW and ADI")
abline(lm(retro_ncc$census_NDI ~ retro_ncc$IPW), col="red")

#ZIP code estimates
zcta_estimates = data.frame("zcta"=unique(retro_ncc$zcta),"adi"=NA,"total_n"=NA,"total_ipw"=NA,"nohcv_n"=NA,"nohcv_ipw"=NA,"yeshcv_n"=NA,"yeshcv_ipw"=NA, "prevalence"=NA, "prevalence_weighted"=NA, "median_ipw"=NA, "median_distance"=NA, "cohort_insurance"=NA, "census_insurance"=NA, stringsAsFactors=F)
for (i in 1:nrow(zcta_estimates)) {
  zcta_estimates$adi[i] = retro_ncc$census_NDI[retro_ncc$zcta==zcta_estimates$zcta[i]][1]
  
  zcta_estimates$total_n[i] = nrow(retro_ncc[retro_ncc$zcta==zcta_estimates$zcta[i], ])
  zcta_estimates$total_ipw[i] = paste(round(median(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i]], na.rm=T), 4), " (",
                                      round(median(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i]], na.rm=T) - IQR(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i]], na.rm=T), 4), ", ",
                                      round(median(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i]], na.rm=T) + IQR(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i]], na.rm=T), 4),
                                      ")", sep="")
  
  zcta_estimates$nohcv_n[i] = nrow(retro_ncc[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==0, ])
  zcta_estimates$nohcv_ipw[i] = paste(round(median(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==0], na.rm=T), 4), " (",
                                      round(median(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==0], na.rm=T) - IQR(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==0], na.rm=T), 4), ", ",
                                      round(median(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==0], na.rm=T) + IQR(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==0], na.rm=T), 4),
                                      ")", sep="")
  
  zcta_estimates$yeshcv_n[i] = nrow(retro_ncc[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==1, ])
  zcta_estimates$yeshcv_ipw[i] = paste(round(median(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==1], na.rm=T), 4), " (",
                                       round(median(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==1], na.rm=T) - IQR(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==1], na.rm=T), 4), ", ",
                                       round(median(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==1], na.rm=T) + IQR(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i] & retro_ncc$hepc==1], na.rm=T), 4),
                                       ")", sep="")
  
  zcta_estimates$prevalence[i] = sum(retro_ncc$hepc[retro_ncc$zcta==zcta_estimates$zcta[i]], na.rm=T) / nrow(retro_ncc[retro_ncc$zcta==zcta_estimates$zcta[i], ]) * 100
  #zcta_estimates$prevalence_weighted[i] = sum(retro_ncc$hepc[retro_ncc$zcta==zcta_estimates$zcta[i]] * retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i]], na.rm=T)
  
  zcta_estimates$median_ipw[i] = median(retro_ncc$IPW[retro_ncc$zcta==zcta_estimates$zcta[i]], na.rm=T)
  
  zcta_estimates$median_distance[i] = median(retro_ncc$distance_to_clinic[retro_ncc$zcta==zcta_estimates$zcta[i]], na.rm=T)
  
  zcta_estimates$cohort_insurance[i] = sum(retro_ncc$insurance_2cat[retro_ncc$zcta==zcta_estimates$zcta[i]]==0, na.rm=T) / nrow(retro_ncc[retro_ncc$zcta==zcta_estimates$zcta[i], ]) * 100
  zcta_estimates$census_insurance[i] = retro_ncc$census_insured_govt_percent[retro_ncc$zcta==zcta_estimates$zcta[i]][1]
}
rm(i)

#prevalence of HCV by zip code
boxplot(zcta_estimates$prevalence) #one potential outlier: 19938
boxplot(zcta_estimates$prevalence[zcta_estimates$zcta!=19938])

#insurance by zip code
describe(zcta_estimates$cohort_insurance); IQR(zcta_estimates$cohort_insurance, na.rm=T)
describe(zcta_estimates$census_insurance); IQR(zcta_estimates$census_insurance, na.rm=T)
boxplot(zcta_estimates$cohort_insurance, zcta_estimates$census_insurance)
zcta_estimates$insurance_diff = zcta_estimates$cohort_insurance - zcta_estimates$census_insurance
plot(zcta_estimates$median_distance, zcta_estimates$insurance_diff)
boxplot(zcta_estimates$insurance_diff, breaks="fd")
describe(zcta_estimates$insurance_diff[zcta_estimates$zcta!=19938]); IQR(zcta_estimates$insurance_diff[zcta_estimates$zcta!=19938], na.rm=T)

#associations of ZIP code ADI and HCV prevalence

#unweighted
plot(zcta_estimates$adi, zcta_estimates$prevalence)
abline(lm(zcta_estimates$prevalence ~ zcta_estimates$adi))

#remove high prev(HCV) outlier: 19938
plot(zcta_estimates$adi[zcta_estimates$zcta!=19938], zcta_estimates$prevalence[zcta_estimates$zcta!=19938])
abline(lm(zcta_estimates$prevalence[zcta_estimates$zcta!=19938] ~ zcta_estimates$adi[zcta_estimates$zcta!=19938]))

#weighted prevalence
zcta_estimates$prevalence_weighted = sum(zcta_estimates$prevalence)*zcta_estimates$prevalence*zcta_estimates$median_ipw/sum(zcta_estimates$median_ipw,na.rm=T)

#weighted
plot(zcta_estimates$adi, zcta_estimates$prevalence_weighted)
abline(lm(zcta_estimates$prevalence_weighted ~ zcta_estimates$adi, weights=zcta_estimates$median_ipw))

#remove high prev(HCV) outlier: 19938
plot(zcta_estimates$adi[zcta_estimates$zcta!=19938], zcta_estimates$prevalence_weighted[zcta_estimates$zcta!=19938])
abline(lm(zcta_estimates$prevalence_weighted[zcta_estimates$zcta!=19938] ~ zcta_estimates$adi[zcta_estimates$zcta!=19938], weights=zcta_estimates$median_ipw[zcta_estimates$zcta!=19938]))

#explore heterogeniety by IPW
plot(zcta_estimates$median_ipw, zcta_estimates$prevalence)
abline(lm(zcta_estimates$prevalence ~ zcta_estimates$median_ipw))

#remove high prev(HCV) outlier: 19938
plot(zcta_estimates$median_ipw[zcta_estimates$zcta!=19938], zcta_estimates$prevalence[zcta_estimates$zcta!=19938])
abline(lm(zcta_estimates$prevalence[zcta_estimates$zcta!=19938] ~ zcta_estimates$median_ipw[zcta_estimates$zcta!=19938]))

#plot of ADI vs prevalence HCV with selection probabities w/outlier removed
# par(xpd=T)
# plot(zcta_estimates$adi, zcta_estimates$prevalence, cex=, xlab="ZIP code ADI", ylab="Observed prevalence of HCV (%)", xlim=c(-2,3), ylim=c(0,4))
# legend(x=-1.7,y=5.5, legend=c("  1", "  1-1.2", "1.2-2.2     ", "  >2.2"), pch=1, pt.cex=, x.intersp=2, bty="n", cex=0.9, title="IPW", horiz=T)
# legend(x=-1,y=4.7, legend=c("Unweighted","Weighted"), lty=c(1,2), bty="n", cex=0.9, horiz=T)
# par(xpd=F)
# abline(lm(zcta_estimates$prevalence[zcta_estimates$zcta!=19938] ~ zcta_estimates$adi[zcta_estimates$zcta!=19938]), lty=1)
# abline(lm(zcta_estimates$prevalence_weighted[zcta_estimates$zcta!=19938] ~ zcta_estimates$adi[zcta_estimates$zcta!=19938]), lty=2)

#plot of ADI vs prevalence HCV with distance w/outlier removed: Figure 2C
zcta_estimates$bubble_size  = as.integer(cut(zcta_estimates$median_distance, quantile(zcta_estimates$median_distance, na.rm=T), include.lowest=TRUE))
zcta_estimates$bubble_size = ifelse(zcta_estimates$bubble_size==1, 0.7, ifelse(zcta_estimates$bubble_size==2, 1.2, ifelse(zcta_estimates$bubble_size==3, 2.8, zcta_estimates$bubble_size)))
plot(zcta_estimates$adi, zcta_estimates$prevalence, cex=zcta_estimates$bubble_size, xlab="ZIP code ADI", ylab="Observed prevalence of HCV (%)", xlim=c(-2,3), ylim=c(0,4))
abline(lm(zcta_estimates$prevalence[zcta_estimates$zcta!=19938] ~ zcta_estimates$adi[zcta_estimates$zcta!=19938]), lty=2)
par(xpd=T)
legend(x=-1.2,y=5.0, legend=c("<4", "4 to 6", "6 to 12   ", " >12"), pch=1, pt.cex=c(0.7,1.2,2.8,4), cex=0.9, title="Distance to FQHC (mi)", horiz=T)
par(xpd=F)

summary(glm(hepc ~ census_NDI, data=retro_ncc, family=binomial()))
confint(glm(hepc ~ census_NDI, data=retro_ncc, family=binomial()))

summary(glm(hepc ~ census_NDI, data=subset(retro_ncc,CA_distance_75per==0), family=binomial()))
confint(glm(hepc ~ census_NDI, data=subset(retro_ncc,CA_distance_75per==0), family=binomial()))

summary(glm(hepc ~ census_NDI, data=subset(retro_ncc,CA_distance_75per==1), family=binomial()))
confint(glm(hepc ~ census_NDI, data=subset(retro_ncc,CA_distance_75per==1), family=binomial()))

#distance vs prevalence w/outlier removed: Figure 2A
plot(zcta_estimates$median_distance, zcta_estimates$prevalence, xlab="Distance to FQHC (mi)", ylab="Observed prevalence of HCV (%)", ylim=c(0,4), pch=16)
abline(lm(zcta_estimates$prevalence[zcta_estimates$zcta!=19938] ~ zcta_estimates$median_distance[zcta_estimates$zcta!=19938]), lty=2)
summary(glm(hepc ~ distance_to_clinic, data=retro_ncc, family=binomial()))
confint(glm(hepc ~ distance_to_clinic, data=retro_ncc, family=binomial()))

#distance vs ADI w/outlier removed: Figure 2B
plot(zcta_estimates$median_distance, zcta_estimates$adi, xlab="Distance to FQHC (mi)", ylab="ZIP code ADI", pch=16)
abline(lm(zcta_estimates$adi[zcta_estimates$zcta!=19938] ~ zcta_estimates$median_distance[zcta_estimates$zcta!=19938]), lty=2)
summary(lm(census_NDI ~ distance_to_clinic, data=retro_ncc))
confint(lm(census_NDI ~ distance_to_clinic, data=retro_ncc))

#note: running influence.measures() on the three ecological models in Figure 2 corroborates that 19938 has excessive influences on the data


### EXAMPLE ANALYSIS: NAIVE and IPW ###

#individual crude associations
describeBy(retro_ncc$last_visit_age, retro_ncc$hepc); wilcox.test(retro_ncc$last_visit_age ~ retro_ncc$hepc)
CrossTable(retro_ncc$sex, retro_ncc$hepc, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(retro_ncc$race, retro_ncc$hepc, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(retro_ncc$ethnicity, retro_ncc$hepc, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(retro_ncc$homeless_any_visit, retro_ncc$hepc, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(retro_ncc$language, retro_ncc$hepc, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(retro_ncc$insurance, retro_ncc$hepc, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(retro_ncc$n_visits, retro_ncc$hepc); wilcox.test(retro_ncc$n_visits ~ retro_ncc$hepc)
describeBy(retro_ncc$days_since_last_visit, retro_ncc$hepc); wilcox.test(retro_ncc$days_since_last_visit ~ retro_ncc$hepc)
CrossTable(retro_ncc$hiv, retro_ncc$hepc, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

#ZIP code crude associations
describeBy(retro_ncc$census_population, retro_ncc$hepc); wilcox.test(retro_ncc$census_population ~ retro_ncc$hepc)
describeBy(retro_ncc$census_income, retro_ncc$hepc); wilcox.test(retro_ncc$census_income ~ retro_ncc$hepc)
describeBy(retro_ncc$census_fertility, retro_ncc$hepc); wilcox.test(retro_ncc$census_fertility ~ retro_ncc$hepc)
describeBy(retro_ncc$census_white_percent, retro_ncc$hepc); wilcox.test(retro_ncc$census_white_percent ~ retro_ncc$hepc)
describeBy(retro_ncc$census_black_percent, retro_ncc$hepc); wilcox.test(retro_ncc$census_black_percent ~ retro_ncc$hepc)
describeBy(retro_ncc$census_other_percent, retro_ncc$hepc); wilcox.test(retro_ncc$census_other_percent ~ retro_ncc$hepc)
describeBy(retro_ncc$census_employed_percent, retro_ncc$hepc); wilcox.test(retro_ncc$census_employed_percent ~ retro_ncc$hepc)
describeBy(retro_ncc$census_insured_percent, retro_ncc$hepc); wilcox.test(retro_ncc$census_insured_percent ~ retro_ncc$hepc)
describeBy(retro_ncc$census_insured_private_percent, retro_ncc$hepc); wilcox.test(retro_ncc$census_insured_private_percent ~ retro_ncc$hepc)
describeBy(retro_ncc$census_insured_govt_percent, retro_ncc$hepc); wilcox.test(retro_ncc$census_insured_govt_percent ~ retro_ncc$hepc)
describeBy(retro_ncc$census_NDI, retro_ncc$hepc); wilcox.test(retro_ncc$census_NDI ~ retro_ncc$hepc)

#check for clustering
model = glmer(hepc ~ (1 | zip), family=binomial(), data=retro_ncc)
summary(model)

#estimates of error terms
print(VarCorr(model),comp=c("Variance","Std.Dev."))
getME(model,"theta") #std dev

#compute median OR, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
#can calculate based on std dev or variance
exp(0.95*getME(model,"theta"))
exp(0.95*sqrt(as.numeric(print(VarCorr(model),comp="Variance"))))

#fully adjusted naive model
model = bglmer(hepc ~ (1 | zcta) + scale(last_visit_age) + as.factor(sex) + as.factor(race_2cat) + as.factor(insurance_2cat) + census_NDI, family=binomial(), data=retro_ncc, control=glmerControl(optimizer="bobyqa"), fixef.prior=normal)

summary(model)
exp(fixef(model))
exp(confint.merMod(model, method="Wald"))

#fully adjusted IPW model
model = bglmer(hepc ~ (1 | zcta) + scale(last_visit_age) + as.factor(sex) + as.factor(race_2cat) + as.factor(insurance_2cat) + census_NDI, family=binomial(), data=retro_ncc, weights=IPW_stabilized, control=glmerControl(optimizer="bobyqa"), fixef.prior=normal)

summary(model)
exp(fixef(model))

#bootstrap for fixed effect CIs since IPW; adapted from https://cran.r-project.org/web/packages/ipw/ipw.pdf (p9)
boot.fun = function(dat, index){
  exp(fixef(bglmer(hepc ~ (1 | zcta) + scale(last_visit_age) + as.factor(sex) + as.factor(race_2cat) + as.factor(insurance_2cat) + census_NDI, family=binomial(), data=dat[index, ], weights=IPW_stabilized, control=glmerControl(optimizer="bobyqa"), fixef.prior=normal)))
}

bootres = boot(retro_ncc, boot.fun, 1000, parallel="multicore", ncpus=6)
boot.ci(bootres, type = "basic", t0=bootres$t0[1], t=bootres$t[,1]) #intercept
boot.ci(bootres, type = "basic", t0=bootres$t0[2], t=bootres$t[,2]) #age
boot.ci(bootres, type = "basic", t0=bootres$t0[3], t=bootres$t[,3]) #sex
boot.ci(bootres, type = "basic", t0=bootres$t0[4], t=bootres$t[,4]) #race
boot.ci(bootres, type = "basic", t0=bootres$t0[5], t=bootres$t[,5]) #insurance
boot.ci(bootres, type = "basic", t0=bootres$t0[6], t=bootres$t[,6]) #adi

# 
# #plot of ADI vs HCV with selection probabities from naive model
# retro_ncc$hepc_pred[complete.cases(retro_ncc[,c("zcta","last_visit_age","sex","race_2cat","insurance_2cat","census_NDI")])] = predict(model, type="response")
# plot_data = data.frame("zcta"=unique(retro_ncc$zcta), "census_NDI"=NA, "hepc_pred"=NA, "IPW"=NA)
# for (i in 1:nrow(plot_data)) {
#   plot_data$census_NDI[i] = mean(retro_ncc$census_NDI[retro_ncc$zcta==plot_data$zcta[i]], na.rm=T)
#   plot_data$hepc_pred[i] = mean(retro_ncc$hepc_pred[retro_ncc$zcta==plot_data$zcta[i]], na.rm=T)
#   plot_data$IPW[i] = mean(retro_ncc$IPW[retro_ncc$zcta==plot_data$zcta[i]], na.rm=T)
# }
# 
# plot(plot_data$census_NDI, plot_data$hepc_pred, cex=quantile(plot_data$IPW,na.rm=T), xlab="ZIP code ADI", ylab="Predicted probability of HCV", xlim=c(-2,3), ylim=c(0,0.1))
# legend_sizes = quantile(plot_data$IPW,na.rm=T)[2-4]
# legend_sizes[5] = legend_sizes[4]
# legend_sizes[4] = 0
# legend("topright", legend=c("1", "1-1.2", "1.2-2.2", "", ">2.2"), pch=1, pt.cex=legend_sizes, x.intersp=2, bty="n", cex=0.9, title="IPW")
# lines(x=c(2.1,2.1),y=c(0.045,0.11))
# lines(x=c(2.1,3.2),y=c(0.045,0.045))

