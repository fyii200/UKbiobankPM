###############################################################
# Script for calculating the prevalence of pathologic myopia  #
# and exploring its risk factors along with fundus biomarkers #
###############################################################
# Author : Fabian Yii                          
# Email  : fabian.yii@ed.ac.uk or fslyii@hotmail.com                       
# 2024

# install.packages("epiR", type="binary")     # Package for computing 95%$ CI for prevalance      
# install.packages("vcd")                     # Package for computing kappa stats
rm(list=ls())
library(ggplot2)
library(cowplot)
library(lme4)
library(sjPlot)
library(epiR)
library(stringr)
library(car)
library(vcd)
source("utils.R")



##### Set-up
# Get project directory
projectDir      <- dirname(getwd())

# Read data
d               <- read.csv(file.path(projectDir, "gradedData", "gradedPMfull.csv"))
d$sleepDuration <- as.numeric(d$sleepDuration)
population      <- subset(d, SER <= -5)
sample          <- subset(d, !is.na(PMbinary_eyeLevel) & PMbinary_eyeLevel != "Reject")
rejected        <- subset(d, PMbinary_eyeLevel == "Reject")



##### Factorise myopic maculopathy categories 
for(i in c("finalMM", "MM_grader1", "MM_grader2")){
  sample[,i]           <- factor(sample[,i])
  if(i == "finalMM"){
    levels(sample[,i]) <- c("M1", "M4", "M2 (macular)", "M0", "M2 (non-macular)", "M3") } 
  else{
    levels(sample[,i]) <- c("NA", "M1", "M4", "M2 (macular)", "M0", "M2 (non-macular)", "M3") 
    } 
  sample[,i]           <- factor(sample[,i], c("M0", "M1", "M2 (non-macular)", "M2 (macular)", "M3", "M4")) 
}



##### Descriptive statistics
# Summary statistics of key variables
printSummary(sample)                                   # Overall (normal + PM)
printSummary(subset(sample, PMbinary_eyeLevel==FALSE)) # Normal eyes
printSummary(subset(sample, PMbinary_eyeLevel==TRUE))  # Eyes with PM
printSummary(rejected)                                 # Ungradable fundus photographs

# Plot: VA by MM category  
ggplot(sample, aes(x = finalMM, y = VA, group = finalMM)) + geom_boxplot()  + labs(x = "\nMyopic maculopathy category") + theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
for(i in levels(sample$finalMM)){
  VAquantiles <- quantile(subset(sample, finalMM==i)$VA, na.rm=TRUE)
  medianVA    <- as.numeric(VAquantiles[3])
  IQRVA       <- VAquantiles[4] - VAquantiles[2]
  print(paste0(i, " median [IQR] VA :", medianVA, " [", IQRVA, "]")) }

# Weighted kappa score (inter-rater agreement)
graderPair <- interaction(sample$initials_grader1, sample$initials_grader2)
kappas     <- c()
for(i in levels(graderPair)){
  rows <- which(graderPair==i)
  if(length(rows)>0){
    grading  <- table(sample[rows, "MM_grader1"], sample[rows, "MM_grader2"])
    kappaRes <- Kappa(grading)
    kappaRes <- as.numeric(kappaRes$Weighted[1])
    weight   <- length(rows) / length(graderPair)
    kappas   <- c(kappas, kappaRes*weight) } 
  }
paste("Weighted inter-grader kappa:", round(sum(kappas),2))

# Contingency table: grader 1 vs grader 2
sum(sample$adjudicated) / nrow(sample)      # 480 images (16%) were adjudicated
table(sample$MM_grader1, sample$MM_grader2) # 241 (50%) & 132 (27.5%) due to M1/non-macular M2 & non-macular M2/macular M2 disagreements



##### Prevalence at individual level
# POPULATION: 3821 participants had high myopia in at least one eye
npopEyes                        <- nrow(population)
npopIndividuals                 <- length(unique(population$id))
sample$PMbinary_individualLevel <- NA
sample$PMbinary_eyeLevel        <- as.logical(sample$PMbinary_eyeLevel)
for(ID in unique(sample$id)){
  sample[which(sample$id==ID),]$PMbinary_individualLevel <- ifelse(sum(subset(sample, id==ID)$PMbinary_eyeLevel)>0, TRUE, FALSE) }
N_cases                         <- sum(sample[!duplicated(sample$id),]$PMbinary_individualLevel)
N_all                           <- length(unique(sample$id))
tmp                             <- as.matrix(cbind(N_cases, N_all))

# Overall prevalence: 41.7% [95% CI: 39.5% to 43.9%]; 831 of 1994 participants affected in at least one eye 
table(sample[!duplicated(sample$id),]$PMbinary_individualLevel)
epi.conf(tmp, ctype = "prevalence", method = "exact", N = npopIndividuals, design = 1, conf.level = 0.95) * 100

# Prevalence by sex
resultPlaceholder               <- printPrevalence(sample, sample$sex, npopIndividuals, eyeLevel = FALSE)

# Prevalence by ethnicity
resultPlaceholder               <- printPrevalence(sample, sample$ethnicBinary, npopIndividuals, eyeLevel = FALSE) # White vs non-White 
resultPlaceholder               <- printPrevalence(sample, sample$ethnic, npopIndividuals, eyeLevel = FALSE)       # Individual ethnic groups



##### Prevalence at eye level
# POPULATION: 5778 eyes had high myopia
npopEyes                                                     <- nrow(population)
N_cases                                                      <- sum(sample$PMbinary_eyeLevel) 
N_all                                                        <- nrow(sample)
tmp                                                          <- as.matrix(cbind(N_cases, N_all))

# Prevalence: 37.9% [95% CI: 36.1% to 39.6%]; 1138 out of 3006 eyes have pathologic myopic
table(sample$PMbinary_eyeLevel)               
epi.conf(tmp, ctype = "prevalence", method = "exact", N = npopEyes, design = 1, conf.level = 0.95) * 100

# Breakdown by MM category, plus lesions and posterior staphyloma
table(sample$finalMM)         # 90 (3.0%), 1780 (59.2%), 1107 (36.8%), 24 (0.8%) and 5 (0.2%) eyes out of 3006 have MM 0, MM 1, MM 2, MM 3 & MM 4. Of the 1107 eyes with M2, 648 (58.5%) and 459 (41.5%) have the non-macular and macular subtypes.  
table(sample$finalMyopicCNV)  # None of the eyes have myopic CNV
table(sample$finalFS)         # 4 eyes have Fuchs spot
table(sample$finalLC)         # 7 eyes have Lacquer cracks
table(sample$finalStaphyloma) # 1 eyes have suspected staphyloma

# PM prevalance in different SER groups
sample$SERgroup                                              <- "-5 to -5.99"
sample[which(sample$SER <= -6 & sample$SER > -7),]$SERgroup  <- "-6 to -6.99"
sample[which(sample$SER <= -7 & sample$SER > -8),]$SERgroup  <- "-7 to -7.99"
sample[which(sample$SER <= -8 & sample$SER > -9),]$SERgroup  <- "-8 to -8.99"
sample[which(sample$SER <= -9 & sample$SER > -10),]$SERgroup <- "-9 to -9.99"
sample[which(sample$SER <= -10),]$SERgroup                   <- "≤ -10"
sample$SERgroup                                              <- factor(sample$SERgroup)
resultPlaceholder                                            <- printPrevalence(sample, sample$SERgroup, npopEyes, eyeLevel = TRUE)

# PM prevalence in different age groups
sample$ageGroup                                              <- "40 to 49"
sample[which(sample$age >= 50 & sample$SER < 60),]$ageGroup  <- "50 to 59"
sample[which(sample$age >= 60 & sample$SER <= 70),]$ageGroup <- "60 to 70"
sample$ageGroup                                              <- factor(sample$ageGroup)
resultPlaceholder                                            <-  printPrevalence(sample, sample$ageGroup, npopIndividuals, eyeLevel = FALSE)

# PM prevalence in different SER by age groups
fig2A                                                        <- plotInteraction(sample, sample$SERgroup, sample$ageGroup, npopEyes, "Age (y)", c(15, 80), xlab = "\n", ylab = "\n")

# PM prevalence in different SER by sex groups
fig2B                                                        <- plotInteraction(sample, sample$SERgroup, sample$sex, npopEyes, "Sex", c(25, 75), xlab = "\n", ylab = "\n")

# PM prevalence in different SER by ethnic groups
fig2C                                                        <- plotInteraction(sample, sample$SERgroup, sample$ethnicBinary, npopEyes, "Ethnicity", c(0, 70), xlab = "\n", ylab = "\n") # White vs non-White by SER

# PM prevalence in different SER by Townsend groups
sample$townsendGroup                                         <- createQuantiles(sample, sample$townsend)
sample$townsendGroup                                         <- factor(sample$townsendGroup) 
levels(sample$townsendGroup)                                 <- c("1 (least)", "2", "3", "4")
fig2D                                                        <- plotInteraction(sample, sample$SERgroup, factor(sample$townsendGroup), npopEyes, " Deprivation\n quantile", c(15, 85), xlab = "\n", ylab = "\n")

# PM prevalence in different SER by FD groups
sample$FDgroup                                               <- createQuantiles(sample, sample$FD_combined)
sample$FDgroup                                               <- factor(sample$FDgroup) 
levels(sample$FDgroup)                                       <- c("1 (least complex)", "2", "3", "4")
fig3A                                                        <- plotInteraction(sample, sample$SERgroup, factor(sample$FDgroup), npopEyes, "Fractal dimension quantile", c(10, 90), xlab = "\n", ylab = "\n")

# PM prevalence in different SER by AVR groups
sample$AVRgroup                                              <- createQuantiles(sample, sample$AVR)
sample$AVRgroup                                              <- factor(sample$AVRgroup) 
levels(sample$AVRgroup)                                      <- c("1 (lowest ratio)", "2", "3", "4")
fig3B                                                        <- plotInteraction(sample, sample$SERgroup, factor(sample$AVRgroup), npopEyes, "Arteriovenous ratio quantile", c(15, 80), xlab = "\n", ylab = "\n")

# PM prevalence in different SER by OD orientation groups
rows                                                         <- !is.na(sample$horODorientation)
group1                                                       <- sample$SERgroup[rows]
group2                                                       <- ifelse(sample$horODorientation, "Horizontal", "Vertical")[rows]
fig3C                                                        <- plotInteraction(sample[rows,], group1, group2, npopEyes, "Optic disc orientation", c(20, 80), xlab = "\n", ylab = "\n")

# Combined plot (Figure 2): stratified by demographic/socioeconomic variables
combinedPlot <- cowplot::plot_grid(fig2A, fig2B, fig2C, fig2D, ncol = 2, axis = "tblr", align = "v")
ggdraw() + 
  draw_plot(combinedPlot) + 
  draw_label("SER (D)", x = 0.5, y = 0, vjust = -0.5, angle = 0, size = 15, fontface = "bold") +
  draw_label("Pevalence (%)", x = 0, y = 0.5, vjust = 1.2, angle = 90, size = 15, fontface = "bold") 
ggsave(file.path(projectDir, "figures", "fig2.png"), width = 10.5, height = 8)

# Combined plot (Figure 3): stratified by fundus biomarkers
combinedPlot <- cowplot::plot_grid(fig3A, fig3B, fig3C, ncol = 1, axis = "tblr", align = "v")
ggdraw() + 
  draw_plot(combinedPlot) + 
  draw_label("SER (D)", x = 0.5, y = 0, vjust = -0.4, angle = 0, size = 15, fontface = "bold") +
  draw_label("Pevalence (%)", x = 0, y = 0.5, vjust = 1.2, angle = 90, size = 15, fontface = "bold")
ggsave(file.path(projectDir, "figures", "fig3.png"), width = 7, height = 10)



##### Risk factors and fundus biomarkers
# Baseline model (SER, age & sex)
sample$sex    <- relevel(factor(sample$sex), "Male") # Male as reference
baselineModel <- glmer(PMbinary_eyeLevel ~ 
                         scale(SER) + 
                         scale(age) + 
                         sex + 
                         (1|id), 
                       family = binomial, 
                       data   = sample) 
tab_model(baselineModel)

# Demographic & lifestyle factors
fitModel(sample, "alcoholFreqBin")              # Alcohol consumption
fitModel(sample, "townsend")                    # Townsend deprivation index    
fitModel(sample, "ethnicBinary")                # Ethnic group
fitModel(sample, "smokeBinary")                 # Smoking status
fitModel(sample, "eduBinary")                   # Highest education level
fitModel(sample, "sleepDuration")               # Sleep duration

# General health indicators
fitModel(sample, "hypertensionCombined")        # Hypertension
fitModel(sample, "cardiovascularDisease")       # Cardiovascular disease
fitModel(sample, "diabetesCombined")            # Diabetes
fitModel(sample, "BMI")                         # Body mass index
fitModel(sample, "totalCholesterol")            # Total serum cholesterol

# Ocular factors
fitModel(sample, "glaucoma")                    # Glaucoma
fitModel(sample, "IOP")                         # IOP

# Fundus biomarkers
fitModel(sample, "FD_combined")                 # Vessel fractal dimension
fitModel(sample, "AVR")                         # Retinal arteriovenous ratio
fitModel(sample, "conc_rp_artery")              # Temporal arterial concavity
fitModel(sample, "conc_rp_vein")                # Temporal venous concavity
fitModel(sample, "Tortuosity_density_combined") # Vessel tortuosity
fitModel(sample, "ODtilt")                      # Optic disc tilt
fitModel(sample, "horODorientation")            # Optic disc orientation

# Final multivariable model
finalModel <- glmer(PMbinary_eyeLevel ~ 
                      scale(SER)         + 
                      scale(age)         + 
                      sex                + 
                      scale(townsend)    + 
                      ethnicBinary       + 
                      scale(FD_combined) + 
                      scale(AVR)         + 
                      horODorientation   + 
                      (1|id), 
                      family  = binomial, 
                      data    = sample, 
                      control = glmerControl(optimizer = "bobyqa", 
                                             optCtrl   = list(maxfun=2e5)) )
tab_model(finalModel)     # Display model
vif(finalModel)           # Variance inflation factor
