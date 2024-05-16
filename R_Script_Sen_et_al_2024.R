################################################################################
# ---------------------------------------------------------------------------- #
# ------------- R script: Benthic remineralization under future -------------- #
# --------------Arctic conditions: evaluating the potential for -------------- #
# -----------changes in carbon sequestration in warming sediments ------------ #
# ------------------------------ Sen et al. 2024 ----------------------------- #
################################################################################
# Script creator: Eric Jorda Molina - 05.02.2024 #
# ---------------------------------------------------------------------------- #
# Retrieving environmental data Nansen Legacy

# Sediment properties (TOC%, N% and C/N)
# Q1: doi: https://doi.org/10.21335/NMDC-1821375519
# Q2: doi: https://doi.org/10.21335/NMDC-350572235
# Q3: doi: https://doi.org/10.21335/NMDC-490057692
# Q4: doi: https://doi.org/10.21335/NMDC-799257283

# Download all data sets for stations P1, P4, P6 and P7 for cruises Q1 (March),
# Q2 (May), Q3 (August) and Q4 (December) to the working directory where this
# R script is save (the files and the script should be in the same folder)

# Set a working directory path relative to where the R script file is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Install and load packages

library(ncdf4) # package for netcdf manipulation
library(dplyr)

################################################################################
######-------------------- SEDIMENT PROPERTIES DATA --------------------########
################################################################################
# Download files to folder where R script is saved
files <- list.files(pattern="*.nc", full.names=TRUE,
                    recursive=FALSE)
# Select only files for sediment properties
files_sedproperties<-files[grepl("properties",files)]

# Create a dataframe with variables to be extracted
sedproperties_data<-data.frame(matrix(ncol = 7, nrow = 0))
colnames(sedproperties_data)<-c("Cruise","Station","BC","Depth","TOC","TC","TN")
# Create loop to extract variable from each file
for (f in files_sedproperties) {
  print(f)
  opened<-nc_open(f)
  depths <- opened$dim[[1]]$vals
  # We split the names of the files to extract the cruise, station, replicate information
  samp <- strsplit(f, "/")[[1]][2]
  samp1<- strsplit(samp, "_")[[1]][1]
  samp2<- strsplit(samp1,"(?<=[0-9])(?=[^0-9])",perl=T)[[1]]

  #########
  # We only extract for the 0-1 cm interval depth (which is the first row) and
  # we take 3 replicates (3 columns) when there are values for each of them and,
  # if not, we exclude them and we only take 2 replicates.
  if (ncol(ncvar_get(opened,"total_carbon"))==3) {
    TOC<-ncvar_get(opened,"total_organic_carbon")[1,1:3]
    TC<-ncvar_get(opened,"total_carbon")[1,1:3]
    TN<-ncvar_get(opened,"total_nitrogen")[1,1:3]


    new_row_1<- c(samp2[1], samp2[2], 1, depths[1],TOC[1], TC[1], TN[1])
    new_row_2<- c(samp2[1], samp2[2], 2, depths[1],TOC[2], TC[2], TN[2])
    new_row_3<- c(samp2[1], samp2[2], 3, depths[1],TOC[3], TC[3], TN[3])
    sedproperties_data <- rbind(sedproperties_data, new_row_1,new_row_2,new_row_3)

  } else {
    TOC<-ncvar_get(opened,"total_organic_carbon")[1,1:2]
    TC<-ncvar_get(opened,"total_carbon")[1,1:2]
    TN<-ncvar_get(opened,"total_nitrogen")[1,1:2]


    new_row_1<- c(samp[1], samp2[2],1,depths[1],TOC[1], TC[1], TN[1])
    new_row_2<- c(samp[1], samp2[2],2, depths[1],TOC[2], TC[2], TN[2])

    sedproperties_data <- rbind(sedproperties_data, new_row_1,new_row_2)
  }
}
head(sedproperties_data)
colnames(sedproperties_data)<-c("Cruise","Station","BC","Depth","TOC","TC","TN")

# For some "Cruise" names it did not completely work so we fix them manually
sedproperties_data$Cruise[sedproperties_data$Cruise == "Q1P7_sediment_properties.nc"] <- "Q1"
# Assign variables to numeric
sedproperties_data$TOC<-as.numeric(sedproperties_data$TOC)
sedproperties_data$TC<-as.numeric(sedproperties_data$TC)
sedproperties_data$TN<-as.numeric(sedproperties_data$TN)
str(sedproperties_data)
#Calculate the average and sd between box core replicates
sedproperties_data_avg<-sedproperties_data%>%
  group_by(Cruise,Station)%>%
  dplyr::summarise(
    across(TOC:TN, list(mean = mean, sd = sd),na.rm = TRUE))

# Calculate C/N ratio
sedproperties_data_avg$CN<-sedproperties_data_avg$TOC_mean/sedproperties_data_avg$TN_mean

# Final dataframe for sediment properties
sedproperties_data_avg



################################################################################
######--------------------- SEDIMENT PIGMENTS DATA ---------------------########
################################################################################
# Sediment pigments (Chl_a, Phaeopigments)
# Q1: doi: https://doi.org/10.11582/2024.00047
# Q2: doi: https://doi.org/10.11582/2024.00046
# Q3: doi: https://doi.org/10.11582/2024.00045
# Q4: doi: https://doi.org/10.11582/2024.00044

# Download files to folder where R script is saved
files <- list.files()
files_pigments<-files[grepl("pigment",files)]

# Create a dataframe with variables to be extracted
pigments_data<-data.frame(matrix(ncol = 6, nrow = 0))
colnames(pigments_data)<-c("Cruise","Station","BC","Depth","Chlorophyll_a",
                           "Phaeopigments")
for (f in files_pigments) {
  opened<-nc_open(f)
  depths <- opened$dim[[1]]$vals
  # We split the names of the files to extract the cruise, station, replicate information
  samp <- strsplit(f, "_")[[1]][7]
  samp2<- strsplit(samp, "\\.")[[1]][1]
  samp3<- strsplit(samp2,"(?<=[0-9])(?=[^0-9])",perl=T)[[1]]
  # We want the 0-1 (0.5) and 1-2 (1.5) cm depth intervals for chla and phaeopigments
  chlorophyll_a_d1<-ncvar_get(opened,"chlorophyll_a_concentration")[1]
  phaeopigments_d1<-ncvar_get(opened,"phaeopigments_concentration")[1]
  chlorophyll_a_d2<-ncvar_get(opened,"chlorophyll_a_concentration")[2]
  phaeopigments_d2<-ncvar_get(opened,"phaeopigments_concentration")[2]
  new_row_d1 <- c(samp3[1], samp3[2], samp3[3], depths[1] ,chlorophyll_a_d1,phaeopigments_d1)
  new_row_d2 <- c(samp3[1], samp3[2], samp3[3], depths[2] ,chlorophyll_a_d2,phaeopigments_d2)
  pigments_data <- rbind(pigments_data, new_row_d1, new_row_d2)
}
head(pigments_data)
colnames(pigments_data) <- c("Cruise","Station","BC","Depth","Chlorophyll_a","Phaeopigments")
head(pigments_data)
pigments_data$Chlorophyll_a<-as.numeric(pigments_data$Chlorophyll_a)
pigments_data$Phaeopigments<-as.numeric(pigments_data$Phaeopigments)

# Sum of 0-1 and 1-2 cm layers for each Chlorophyll_a and Phaeopigments
pigments_data_sum_0_2<-pigments_data%>%
  group_by(Cruise,Station,BC)%>%
  dplyr::summarise(
    across(c(Chlorophyll_a,Phaeopigments),sum))
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#


################################################################################
####### - Macrofauna dataset Nansen Legacy - Jorda-Molina et al.,2024 - #######
### - Retrieved from NMDC - doi: https://doi.org/10.21335/NMDC-1152502405 - ####
################################################################################
# Script creator: Eric Jorda  Molina - 05.02.2024 #
# ---------------------------------------------------------------------------- #

# Loading of R packages:
library(tidyr)
library(dplyr)
library(ggplot2)
library(upstartr)
library(ggh4x)
library(stringr)
library(utils)

# First download the macrofauna dataset.
# The dataset is also accessible in GBIF as a Darwin Core file:
# https://gbif.imr.no/ipt/resource?r=aen_benthic_macrofauna#anchor-description

# Import the occurrence.txt file
# We have to take into account the spaces in the values of some columns by using
# comment.char=""!
macrofauna<-read.table("occurrence.txt",sep="\t",  skip = 0, header = TRUE,
                       comment.char = "", check.names = FALSE)
# We exclude formanifera from the dataset:
macrofauna<-macrofauna[!grepl("Foraminifera", macrofauna$occurrenceID),]
# The respiration incubations were only carried out in stations P1, P4, P6 and
# P7, so we filter by those stations
macrofauna_inc<-filter(macrofauna, locationID =="P1"|locationID =="P4"|
                         locationID =="P6"|locationID =="P7")

# Now we select only the cores in which we assessed the ambient respiration
# rates (T1 treatment) and both T1 and T3 for Q4_P6
macrofauna_T1<-macrofauna_inc%>%
  filter(str_detect(occurrenceID,"_T1|AeN Q4_P6"))
# We remove Q1_P7 because we do not have respiration data for that sample
macrofauna_T1<-macrofauna_T1%>%
  filter(!str_detect(occurrenceID,"AeN Q1_P7"))


# Now we split the matrialSampleID column into new columns:
# "Cruise", "Station", "BoxCore" and "Replicate"
macrofauna_T1<-macrofauna_T1%>%
  separate(materialSampleID, c("Cruise", "Station","BoxCore","Replicate"), "_")

# Now we subset for the abundance variable (individualCount in the dataset):

Abundance<-macrofauna_T1[c("Cruise","Station","Replicate","individualCount",
                           "scientificName","identificationQualifier")]
Abundance<-Abundance[!(Abundance$Cruise=="AeN Q4"&Abundance$Station=="P6"&
                         grepl("^T3", Abundance$Replicate)),]


# And we subset for the biomass variable (organismQuantity in the dataset):

Biomass<-macrofauna_T1[c("Cruise","Station","Replicate","organismQuantity",
                           "scientificName","identificationQualifier")]

Biomass<-Biomass[!(Biomass$Cruise=="AeN Q4"&Biomass$Station=="P6"&
                     grepl("^T1", Biomass$Replicate)),]

# Now we remove the fragments of organisms for the abundance subset
Abundance_clean<-Abundance[!grepl("fragments",
                                  Abundance$identificationQualifier),]



library(writexl)
write_xlsx(Abundance_clean,"./Abundance_macrofauna.xlsx")
write_xlsx(Biomass,"./Biomass_macrofauna.xlsx")

#####---------- CONVERTING DATAFRAME TO MATRIX FORM FOR ABUNDANCE ----------####
################################################################################
### Converting Abundance dataset into matrix form considering each replicate ###
# The first step is to ensure that there are no duplicates of taxa
# (there shouldn't be).
# This way, if there are, they are summed for each station/replicate
AeN_abu_matrix<-Abundance_clean%>%
  pivot_wider(id_cols=c(Cruise,Station,Replicate),names_from="scientificName",
          values_from="individualCount",values_fn=list("individualCount"=sum))
AeN_abu_matrix[is.na(AeN_abu_matrix)] = 0

AeN_abu_matrix<-as.data.frame(AeN_abu_matrix)
AeN_abu_matrix<-AeN_abu_matrix %>% unite("ID",c(Cruise,Station,Replicate),
                                         remove=T)
row.names(AeN_abu_matrix)=AeN_abu_matrix$ID
AeN_abu_matrix<-AeN_abu_matrix[,-1]


# -----------------------------------------------------------------------------#
AeN_abu_matrix$Sample=row.names(AeN_abu_matrix)
#First we put the species from the matrix form into the table form again so that
#  all species are the same in all samples before we do the average per species
table_AeN_matrix<-pivot_longer(AeN_abu_matrix,cols=-Sample,names_to="Species",
                               values_to="Abundance")
table_AeN_matrix[c("Season", "Station","Replicate")] <-
  str_split_fixed(table_AeN_matrix$Sample, "_", 3)
table_AeN_matrix$Season_Station<-paste(table_AeN_matrix$Season,
                                       table_AeN_matrix$Station,sep="_")

# Now we do the average by species for Station and Season #
AeN_avg_abu<-table_AeN_matrix%>%
  group_by(Season,Station,Species)%>%
  dplyr::summarise(
    across(Abundance,~mean(.))
  )%>%glimpse


# Select the five most abundant species per season/station
five_most_abundant<-AeN_avg_abu %>%
  group_by(Season,Station)%>%
  arrange(desc(Abundance)) %>%
  slice(1:5)


#####---------- CONVERTING DATAFRAME TO MATRIX FORM FOR BIOMASS ------------####
################################################################################
#### Converting Biomass dataset into matrix form considering each replicate ####
AeN_bio_matrix<-Biomass%>%
  pivot_wider(id_cols=c(Cruise,Station,Replicate),names_from="scientificName",
              values_from="organismQuantity",
              values_fn=list("organismQuantity"=sum))
AeN_bio_matrix[is.na(AeN_bio_matrix)] = 0

AeN_bio_matrix<-as.data.frame(AeN_bio_matrix)
AeN_bio_matrix<-AeN_bio_matrix %>% unite("ID",c(Cruise,Station,Replicate),
                                         remove=T)
row.names(AeN_bio_matrix)=AeN_bio_matrix$ID
AeN_bio_matrix<-AeN_bio_matrix[,-1]



# -----------------------------------------------------------------------------#
AeN_bio_matrix$Sample=row.names(AeN_bio_matrix)
# First we put the species from the matrix form into the table form again so
# that all speciesare the same in all samples before we do the average
table_AeN_matrix<-pivot_longer(AeN_bio_matrix,cols=-Sample,names_to="Species",
                               values_to="Weight")
table_AeN_matrix[c("Season", "Station","Replicate")] <-
  str_split_fixed(table_AeN_matrix$Sample, "_", 3)
table_AeN_matrix$Season_Station<-paste(table_AeN_matrix$Season,
                                       table_AeN_matrix$Station,sep="_")

# Now we do the average by species for Station and Season #
AeN_avg_bio<-table_AeN_matrix%>%
  group_by(Season,Station,Species)%>%
  dplyr::summarise(
    across(Weight,~mean(.))
  )%>%glimpse

# Retrieve the 5 top taxa for biomass
five_most_biomass<-AeN_avg_bio %>%
  group_by(Season,Station)%>%
  arrange(desc(Weight)) %>%
  slice(1:5)


################################################################################
###### ------- Plotting the five most abundant and biomass taxa -------- #######

plot_five_abundance<-ggplot(five_most_abundant,aes(y=Abundance,
        x=reorder_within(Species,Abundance,list(Station,Season)),fill=Season))+
  geom_bar(position=position_dodge(width=1),stat="identity")+
 facet_grid2(Season~Station, scales = "free",independent="y",drop=TRUE)+
  coord_flip()+
  scale_x_reordered()

plot_five_biomass<-ggplot(five_most_biomass,aes(y=Weight,
          x=reorder_within(Species,Weight,list(Station,Season)),fill=Season))+
  geom_bar(position=position_dodge(width=1),stat="identity")+
  facet_grid2(Season~Station, scales = "free",independent="y",drop=TRUE)+
  coord_flip()+
  scale_x_reordered()
# ---------------------------------------------------------------------------- #

### Now calculate the sum for abundance and biomass for each station/season ###
total_abundance<-AeN_avg_abu%>%
  group_by(Season,Station)%>%
  dplyr::summarise(
    across(Abundance,~sum(.))
  )%>%glimpse


total_biomass<-AeN_avg_bio%>%
  group_by(Season,Station)%>%
  dplyr::summarise(
    across(Weight,~sum(.))
  )%>%glimpse


write_xlsx(total_abundance,"./Total_abundance_macrofauna.xlsx")
write_xlsx(total_biomass,"./Total_Biomass_macrofauna.xlsx")


################################################################################
### ---------------- MULTIPLE REGRESSION WITH SOD RATES -------------------- ###
################################################################################
#use averaged abundance and biomass from macrofauna dataset above and sediment
#parameters to make the multiple regression with SOD rates

#multiple regression:
library(tidyverse)
# csv file submitted
mult_regresssion <- read.csv("./mult_regresssion.csv", sep = ";", header = TRUE)

#model coefficients:
model <- lm(rate ~ chla + phaeo + TOC + N + C.N + abundance + biomass + temp, data = mult_regresssion)
summary(model)


###############################################################################
#### ---------- PLOT FOR LOG RESPONSE RESPIRATION RATES --------------- #######
###############################################################################
# All calculations for log response ratios are done in excel
# Respiration results #

resp_treatments<-read.table("Treatments_rep.txt",header=T,sep="\t",dec=",")
resp_treatments$Treatment<-as.factor(resp_treatments$Treatment)

resp_treatments$Treatment  <- factor(resp_treatments$Treatment, levels =
                                       c("WARM","ALGAE","ALGAE_WARM"))
resp_treat_plot<-ggplot(resp_treatments, aes(x=Treatment, y=RR., color =
                                               Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(resp_treatments,mapping=aes(ymin = X95.CIlower,
                                            ymax = X95.CIupper),
                width = 0.25, size = 1) +
  scale_colour_manual(values = c("firebrick2", "forestgreen", "deepskyblue2")) +
  theme(axis.title.x = element_blank()) +
  facet_grid(Season~Station,scales="free_y")+expand_limits(y = 0)+
  geom_hline(yintercept=100)+
  scale_y_continuous(trans='log10')
