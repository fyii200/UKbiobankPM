###############################################################
#   Script for preprocessing variables analysed in the study  #
###############################################################
# Author : Fabian Yii                          
# Email  : fabian.yii@ed.ac.uk or fslyii@hotmail.com                       
# 2024

rm(list=ls())
library(stringr)



##### Set-up
# Get project directory
projectDir  <- dirname(getwd())

# Read participant, blood measure and retinal feature data
d           <- read.csv(file.path(projectDir, "participantData", "cleaned_data_long_PM_cohort.csv"))
blood       <- read.csv(file.path(projectDir, "participantData", "cleaned_data_bloodMeasuresAndBMI_all.csv"))
retina      <- read.csv(file.path(projectDir, "participantData", "fullDerivedData.csv"))

# Merge dataframes
d <- merge(d, blood, by = c("id", "age", "YOB", "sex", "townsend"))
d <- merge(d, retina, by = "fundus")



##### Preprocess demographic or socioeconomic variables 
# Binarise education level into "Post-secondary" and "≤ Secondary"
for(i in 1:nrow(d)){
  splitText        <- unlist(strsplit(d$edu[i], ""))
  if(sum(splitText=="|")>0){
    splitInd       <- which(splitText == "|")[1]
    d$edu[i]       <- str_c(splitText[1:splitInd-1], collapse = "") } }
d$eduBinary        <- factor(ifelse(d$edu == "College or University degree" | d$edu == "Other professional qualifications eg: nursing, teaching" | d$edu == "NVQ or HND or HNC or equivalent", "Post-secondary", "≤ Secondary"))

# Group "Indian" and "Pakistani" together
d$ethnic           <- ifelse(d$ethnic == "Indian" | d$ethnic == "Pakistani", "Indian/Pakistani", d$ethnic)                      

# Binarise ethnicity into "White" and "non-White"
d$ethnicBinary     <- ifelse(d$ethnic=="British" | d$ethnic=="Irish" | d$ethnic=="White" | d$ethnic=="Any other white background", "White", "non-White")



##### Preprocess fundus imaging features 
# Compute retinal arteriovenous ratio (AVR)
d$AVR              <- d$CRAE_Knudtson / d$CRVE_Knudtson

# Compute and binarise OD tilt using the widely adopted cutoff of 1.3
d$ODtilt           <- (d$major_length/d$minor_length) >= 1.3 

# Compute and binarise OD orientation
d$horODorientation <- ifelse(abs(d$orientation) < 45, TRUE, FALSE)



##### Preprocess lifestyle variables
# Binarise smoking frequency into "Never/previous" and "Current"
d$smokeBinary      <- ifelse(d$smoke == "Prefer not to answer" | d$smoke == "" | d$smoke == "Previous" | d$smoke == "Never", "Never/previous", "Current")

# Binarise alcohol consumption into "frequent" and "infrequent", using "daily or almost daily" as cut-off
d$alcoholFreqBin   <- factor(ifelse(d$alcoholFreq == "Daily or almost daily", "Frequent", "Infrequent"))



##### Preprocess health-related variables
# Cardiovascular disease defined based on linked healthcare data
d$cardiovascularDisease                                <- d$myocardialInfarction | d$cardiomyopathy | d$ischaemicHeartDisease | d$cardiacArrest | d$heartFailure | d$multipleValvularHeartDisease | d$atherosclerosis | d$otherHeartDiseases

# Hypertension defined based on linked healthcare data as well as systolic and diastolic blood pressure
d$hypertensionCombined                                 <- d$SysBP >= 140 | d$DiasBP >= 90 | d$hypertension
d[is.na(d$hypertensionCombined),]$hypertensionCombined <- d[is.na(d$hypertensionCombined),]$hypertension

# Diabetes defined based on linked healthcare data and random glucose (mmol/L) 
# Cut-off: 200 mg/dL (200*0.0555 = 11.1mmol/L)
d$diabetesCombined                                     <- d$glucose>11.1 | d$diabetes
d[is.na(d$diabetesCombined),]$diabetesCombined         <- d[is.na(d$diabetesCombined),]$diabetes



##### Remove ungradable images, as determined by the reading centre and during adjudication 
# Read PM grading data (3024 eyes of 2000 unique participants)
grading                      <- read.csv(file.path(projectDir, "gradedData", "cleanedGrading.csv"))

# Create a column indicating if a given image has been adjudicated
grading$adjudicated          <- ifelse(grading$adjudicatedMM_MJQ != "" | grading$adjudicatedMM_IJCM != "" | grading$adjudicatedPlusLesions != "", TRUE, FALSE)

# Image quality = "Include" if graded as "Good" or "Usable"
grading$imageQuality_grader1 <- ifelse(grading$imageQuality_grader1 != "Reject", "Include", "Reject")
grading$imageQuality_grader2 <- ifelse(grading$imageQuality_grader2 != "Reject", "Include", "Reject")

# Image quality: 7 images were unanimously rejected by the two graders, leaving 3017 images
rejectBool                   <- grading$imageQuality_grader1 == "Reject" & grading$imageQuality_grader2=="Reject"
cleanData                    <- grading[!rejectBool,] 
cleanData                    <- merge(cleanData, d, by = c("fundus", "eye"))

# Image quality: 11 images were further rejected during clinical adjudciation, leaving 3006 images
cleanData                    <- cleanData[cleanData$adjudicatedMM_IJCM != "Reject" & cleanData$adjudicatedMM_MJQ != "Reject",]



##### Add final grading
# Create a new column indicating final MM grade for each image
cleanData$finalMM                <- ""
for(i in 1:nrow(cleanData)){
  if(cleanData$adjudicated[i]==FALSE){
    cleanData$finalMM[i]         <- cleanData$MM_grader1[i]
  } else{
    cleanData$finalMM[i]         <- ifelse(cleanData$adjudicatedMM_IJCM[i] == "", cleanData$adjudicatedMM_MJQ[i], cleanData$adjudicatedMM_IJCM[i])
    } 
  }

# "plus" lesion 1: Myopic CNV
cleanData$myopicCNV_grader1      <- ifelse(is.na(cleanData$myopicCNV_grader1), "NA", cleanData$myopicCNV_grader1)
cleanData$myopicCNV_grader2      <- ifelse(is.na(cleanData$myopicCNV_grader2), "NA", cleanData$myopicCNV_grader2)
cleanData$finalMyopicCNV         <- ""
for(i in 1:nrow(cleanData)){
  if(cleanData$adjudicated[i] == TRUE & cleanData$adjudicatedMyopicCNV[i] == "yes"){
    cleanData$finalMyopicCNV[i]  <- "yes"
  } else if(cleanData$adjudicated[i] == TRUE & cleanData$adjudicatedMyopicCNV[i] == "no" | cleanData$adjudicatedMyopicCNV[i] == ""){
    cleanData$finalMyopicCNV[i]  <- "no"
  } else if(cleanData$adjudicated[i] == FALSE & cleanData$myopicCNV_grader1[i] == cleanData$myopicCNV_grader2[i]) { 
    cleanData$finalMyopicCNV[i]  <- cleanData$myopicCNV_grader1[i]
  } else{
    cleanData$finalMyopicCNV[i]  <- "Need adjudication"} }

# "plus" lesion 2: Fuchs spots
cleanData$FS_grader1             <- ifelse(is.na(cleanData$FS_grader1), "NA", cleanData$FS_grader1)
cleanData$FS_grader2             <- ifelse(is.na(cleanData$FS_grader2), "NA", cleanData$FS_grader2)
cleanData$finalFS                <- ""
for(i in 1:nrow(cleanData)){
  if(cleanData$adjudicated[i] == TRUE & cleanData$adjudicatedFS[i] == "yes"){
    cleanData$finalFS[i]         <- "yes"
  } else if(cleanData$adjudicated[i] == TRUE & cleanData$adjudicatedFS[i] == "no" | cleanData$adjudicatedFS[i] == ""){
    cleanData$finalFS[i]         <- "no"
  } else if(cleanData$adjudicated[i] == FALSE & cleanData$FS_grader1[i] == cleanData$FS_grader2[i]) { 
    cleanData$finalFS[i]         <- cleanData$FS_grader1[i]
  } else{
    cleanData$finalFS[i]         <- "Need adjudication" } }

# "plus" lesion 3: Lacquer cracks
cleanData$LC_grader1             <- ifelse(is.na(cleanData$LC_grader1), "NA", cleanData$LC_grader1)
cleanData$LC_grader2             <- ifelse(is.na(cleanData$LC_grader2), "NA", cleanData$LC_grader2)
cleanData$finalLC                <- ""
for(i in 1:nrow(cleanData)){
  if(cleanData$adjudicated[i] == TRUE & cleanData$adjudicatedLC[i] == "yes"){
    cleanData$finalLC[i]         <- "yes"
  } else if(cleanData$adjudicated[i] == TRUE & cleanData$adjudicatedLC[i] == "no" | cleanData$adjudicatedLC[i] == ""){
    cleanData$finalLC[i]         <- "no"
  } else if(cleanData$adjudicated[i] == FALSE & cleanData$LC_grader1[i] == cleanData$LC_grader2[i]) { 
    cleanData$finalLC[i]         <- cleanData$LC_grader1[i]
  } else{
    cleanData$finalLC[i]         <- "Need adjudication" } }

# Posterior staphyloma
cleanData$staphyloma_grader1     <- ifelse(is.na(cleanData$staphyloma_grader1), "NA", cleanData$staphyloma_grader1)
cleanData$staphyloma_grader2     <- ifelse(is.na(cleanData$staphyloma_grader2), "NA", cleanData$staphyloma_grader2)
cleanData$finalStaphyloma        <- ""
for(i in 1:nrow(cleanData)){
  if(cleanData$adjudicated[i] == TRUE & cleanData$adjudicatedStaphyloma[i] == "yes"){
    cleanData$finalStaphyloma[i] <- "yes"
  } else if(cleanData$adjudicated[i] == TRUE & cleanData$adjudicatedStaphyloma[i] == "no" | cleanData$adjudicatedStaphyloma[i] == ""){
    cleanData$finalStaphyloma[i] <- "no"
  } else if(cleanData$adjudicated[i] == FALSE & cleanData$staphyloma_grader1[i] == cleanData$staphyloma_grader2[i]) { 
    cleanData$finalStaphyloma[i] <- cleanData$staphyloma_grader1[i]
  } else{
    cleanData$finalStaphyloma[i] <- "Need adjudication" } }

# Create new labels indicating the presence of MM, any of the 3 plus lesion(s) and staphyloma, and PM
cleanData$finalMMbinary          <- ifelse(cleanData$finalMM != "No MM - 0" & cleanData$finalMM != "Fundus tessellation only - 1", TRUE, FALSE)
cleanData$finalPlusLesionBinary  <- ifelse(cleanData$finalMyopicCNV == "no" & cleanData$finalFS == "no" & cleanData$finalLC == "no" & cleanData$finalStaphyloma == "no", FALSE, TRUE)
cleanData$PMbinary_eyeLevel      <- cleanData$finalMMbinary | cleanData$finalPlusLesionBinary



##### Check inter-grader agreement for PM binary classification
MMgrader1 <- cleanData$MM_grader1 != "No MM - 0" & cleanData$MM_grader1 != "Fundus tessellation only - 1" 
PMgrader1 <- MMgrader1 | cleanData$FS_grader1 == "yes" | cleanData$LC_grader1 == "yes" | cleanData$myopicCNV_grader1 == "yes" | cleanData$staphyloma_grader1 == "yes"
MMgrader2 <- cleanData$MM_grader2 != "No MM - 0" & cleanData$MM_grader2 != "Fundus tessellation only - 1"
PMgrader2 <- MMgrader2 | cleanData$FS_grader2 == "yes" | cleanData$LC_grader2 == "yes" | cleanData$myopicCNV_grader2 == "yes" | cleanData$staphyloma_grader2 == "yes"
sum(PMgrader1 == PMgrader2, na.rm = TRUE) / nrow(cleanData)



##### Write dataframe
gradedPMfull           <- merge(d, cleanData[,c(1,3:29,217:225)], all=TRUE, by="fundus")
names(gradedPMfull)[1] <- "name"
rejectNames            <- setdiff(grading$fundus, cleanData$fundus)
for(name in rejectNames){
  gradedPMfull[which(gradedPMfull$name == name),]$PMbinary_eyeLevel <- "Reject"
}
write.csv(gradedPMfull, "gradedPMfull.csv")



























