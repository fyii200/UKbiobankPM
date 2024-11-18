##############################################################
# Script defining the various functions used in "analysis.R  #
##############################################################
# Author : Fabian Yii                          
# Email  : fabian.yii@ed.ac.uk or fslyii@hotmail.com                       
# 2024



## FUNCTION 1: print summary statistics (mean ± SD for continuous variables).
##   Eye-level variables are summarised at eye level, while person-specific
##   variables are summarised at individual level.
printSummary <- function(data){
  IndividualVariables <- c("age", "sex", "townsend", "eduBinary", "ethnicBinary", "smokeBinary", "alcoholFreqBin", "sleepDuration", "hypertension", "cardiovascularDisease", "diabetes", "BMI", "totalCholesterol")
  eyeVariables        <- c("SER", "glaucoma", "VA", "IOP", "FD_combined", "Tortuosity_density_combined", "conc_rp_artery", "conc_rp_vein", "AVR", "ODtilt", "horODorientation")
  
  # If person-level variable make sure only one row is included
  for(i in IndividualVariables){
    indData <- data[!duplicated(data$id),]
    if(class(indData[[i]]) == "character" | class(indData[[i]]) == "logical"){
      print(paste("===============================", i, "==============================="))
      print(table(indData[[i]])) }
    else{
      meanVal <- round(mean(indData[[i]], na.rm=TRUE), 5)
      sdVal   <- round(sd(indData[[i]], na.rm=TRUE), 5)
      print(paste("===============================", i, "==============================="))
      print(paste(meanVal, "±", sdVal)) } }
  
  # Eye-level variable
  for(i in eyeVariables){
    if(class(data[[i]]) == "character" | class(data[[i]]) == "logical"){
      print(paste("===============================", i, "==============================="))
      print(table(data[[i]])) }
    else{
      meanVal <- round(mean(data[[i]], na.rm=TRUE), 5)
      sdVal   <- round(sd(data[[i]], na.rm=TRUE), 5)
      print(paste("===============================", i, "==============================="))
      print(paste(meanVal, "±", sdVal)) } }
}



## FUNCTION 2: Divide continuous variable into quantiles.
##   Return dataframe with a new column indicating quantile membership.
createQuantiles <- function(data, variable){
  data$group <- 1
  quantiles <- quantile(variable, seq(0.25, 0.75, 0.25), na.rm=TRUE)
  for(i in 2:4){
    if(i != 4){
      lowerLim <- as.numeric(quantiles[i-1])
      upperLim <- as.numeric(quantiles[i])
      data[which(variable >= lowerLim & variable < upperLim),]$group <- i } 
    else{ 
      data[which(variable >= as.numeric(quantiles[i-1])),]$group <- i } 
  }
  return(data$group)
}



## FUNCTION 3: Compute and print prevalence ± 95% CI by "group".
##   Return dataframe with new columns added indicating the
##   group-specific estimated prevalence ± 95% CI.
printPrevalence <- function(data, group, Npop, eyeLevel){
  data$group       <- factor(group)
  data$PMprev      <- NA
  data$PMprevLower <- NA
  data$PMprevUpper <- NA
  for(i in levels(data$group)){
    if(eyeLevel){
      N_cases    <- sum(data[which(data$group == i),]$PMbinary_eyeLevel)
      N_all      <- sum(data$group == i)
    } else{
      uniqueData <- data[!duplicated(data$id) & data$group==i,]
      N_cases    <- sum(uniqueData$PMbinary_individualLevel)
      N_all      <- nrow(uniqueData)
    }
    tmp         <- as.matrix(cbind(N_cases, N_all))
    prev        <- epi.conf(tmp, 
                            ctype      ="prevalence", 
                            method     = "exact", 
                            N          = Npop, 
                            design     = 1, 
                            conf.level = 0.95) * 100
    print(paste0("PM prevalence is ", 
                 round(prev[[1]],1), "% [", 
                 round(prev[[2]],1), " to ", 
                 round(prev[[3]],1), "] in ", 
                 i, " group")) 
    
    data[which(data$group==i),]$PMprev      <- as.numeric(prev[[1]]) 
    data[which(data$group==i),]$PMprevLower <- as.numeric(prev[[2]])
    data[which(data$group==i),]$PMprevUpper <- as.numeric(prev[[3]]) }
  return(data)
  }



## FUNCTION 4: Compute, print and plot prevalence ± 95% CI by interaction group.
plotInteraction <- function(data, group1, group2, NpopEyes, legendTitle, ylim = c(20,80), xlab = "\nSpherical equivalent refraction", ylab = "Prevalence (%)\n"){
  
  # Compute and print prevalence by interaction group
  print("=============== By interaction group ==============")
  data$interacted   <- interaction(group1, group2)
  data              <- printPrevalence(data, data$interacted, NpopEyes, eyeLevel = TRUE)
  
  # PLOT
  bgCol     <- rgb(0.6,0.1,0, alpha=0.01)
  fgCol     <- rgb(0.8,0.4,0, alpha=0.03)
  cbPalette <- c("#E69F00", "#009E73", "#D55E00", "#CC79A7") # colour blind palette
  ggplot(data) + 
    
    theme_minimal() +
    
    geom_errorbar(aes(x     = group1, 
                      y     = PMprev, 
                      color = group2, 
                      ymin  = PMprevLower, 
                      ymax  = PMprevUpper), 
                  width     = 0.3,
                  size      = 1.02,
                  position  = position_dodge(.4)) +
    
    geom_point(aes(x     = group1, 
                   y     = PMprev, 
                   color = group2),
               size = 2.5,
               position  = position_dodge(.4)) +
    
    theme(text               = element_text(family = "Times New Roman"),
          panel.background   = element_rect(fill = bgCol, color = NA),
          panel.border       = element_blank(),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          axis.text.x        = element_text(size = 10, face = "bold"),
          axis.text.y        = element_text(size = 12),
          legend.text        = element_text(size = 12),
          legend.title       = element_text(size = 13),
          legend.margin      = margin(t = -8, b  = -15),
          legend.key.size    = unit(1, "cm"),
          legend.spacing.y   = unit(0, "pt"),
          legend.position    = "top") +
    
    scale_y_continuous(breaks = seq(ylim[1], ylim[2], 10), 
                       limits = ylim) +
    
    scale_color_manual(name   = legendTitle,
                       values = cbPalette) +
    
    labs(x = xlab, 
         y = ylab)
}



## FUNCTION 5: Fit and display a random-effects logistic regression model, 
##   with "variable" as the independent variable and pathologic myopia as
##   the dependent variable, controlling for SER, age and sex. 
fitModel <- function(data, variable, scale=TRUE){
  # If scaling is required
  if(scale){
    # Apply scaling if variable is numeric
    if(is.numeric(data[,variable])){
      formula  <- as.formula(paste("PMbinary_eyeLevel ~ scale(SER) + scale(age) + sex +", paste0("scale(", variable, ")"), "+ (1|id)")) 
    } else{
      formula  <- as.formula(paste("PMbinary_eyeLevel ~ scale(SER) + scale(age) + sex +", variable, "+ (1|id)"))
    }
  }
  # If scaling is not required
  else{
    formula <- as.formula(paste("PMbinary_eyeLevel ~ SER + age + sex +", variable, "+ (1|id)"))
  }
  model <- glmer(formula, 
                 family  = binomial, 
                 data    = sample,
                 control = glmerControl(optimizer = "bobyqa", 
                                        optCtrl   = list(maxfun=2e5)) )
  return(tab_model(model))
}

