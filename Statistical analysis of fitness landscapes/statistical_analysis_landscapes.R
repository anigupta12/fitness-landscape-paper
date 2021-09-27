## Loading the required packages
library(dplyr)
library(ggplot2)
library(writexl)
library(car)
library(usdm)
require(graphics)
library(Hmisc)
library(MASS)
library(corrgram)
library(corrplot)
library("corrgram")

## RUN 'data_format_for_R_design_matrix.m' to get the input files '606.csv' and 'malT.csv' for this

# Load design matrix from design_matrix_for_R_FINAL_v2.mat -> output of matlab file 'data_format_for_R_design_matrix.m'
# library(R.matlab) 
# data <- readMat('design_matrix_for_R_FINAL_v2.mat')

# Loading identity of genotypes
wt <- read.csv('606.csv', header = FALSE)
malT <- read.csv('malT.csv', header = FALSE)
data$design.matrix.genotype.wt <- wt
data$design.matrix.genotype.malT <- malT

# Creating design matrix in R format
fitness_input = rbind(data$design.matrix.fitness.wt,data$design.matrix.fitness.malT)
genotypeID_input = rbind(data$design.matrix.genotypeID.wt,data$design.matrix.genotypeID.malT)
host_input = rbind(data$design.matrix.host.wt,data$design.matrix.host.malT)
genotype = rbind(data$design.matrix.genotype.wt,data$design.matrix.genotype.malT)
design_matrix <- data.frame("fitness" = fitness_input, "genotypeID" = genotypeID_input, "genotype" = genotype, "host" = host_input)

# Checking design matrix
head(design_matrix)
str(design_matrix)
tail(design_matrix)

## Create all two-way interactions with 10 sites for mutation and 1 host type -> 55 terms + 11 main interaction terms
X_matrix <- lapply(design_matrix[,-c(1,2)], factor)  # convert all columns into factors
two_way_model_matrix <- model.matrix(object = ~.^2, data = X_matrix)
two_way_model_matrix <- two_way_model_matrix[,-c(1)]
str(two_way_model_matrix)
head(two_way_model_matrix)
ncol(two_way_model_matrix)

# Creating a matrix with only 10 sites of mutations
X_site_matrix <- lapply(design_matrix[,-c(1,2,13)], factor)  # removing intercept, genotypeID and fitness column
str(X_site_matrix) 

# to get all two-way interaction between 10 sites -> 45 terms, 10C2 = 45
two_way_site_model_matrix <- model.matrix(object = ~.^2, data = X_site_matrix)
two_way_site_model_matrix <- two_way_site_model_matrix[,-c(1:11)]  # removing all one-way main effect terms and intercept term
str(two_way_site_model_matrix)
head(two_way_site_model_matrix)

# to get all the two-way sites with host -> 45 interaction terms,  10C2 * 1 = 45
model_matrix_2site_1host <- model.matrix(object = ~(.):two_way_model_matrix[,"hostwt"], data = data.frame(two_way_site_model_matrix))
model_matrix_2site_1host <- model_matrix_2site_1host[,-c(1)]  # removing intercept term
head(model_matrix_2site_1host)
str(model_matrix_2site_1host)

# model matrix without fitness -> 111 terms [[ 11 one-way (10 site + 1 host), 55 two-way terms (11C2 = 55), 45 two-way site by 1 host terms (10C2 * 1 = 45) ]]
model_matrix_final <- data.frame(cbind(two_way_model_matrix,model_matrix_2site_1host))
str(model_matrix_final, list.len=ncol(model_matrix_final))
ncol(model_matrix_final)

### FINAL data frame to use linear regression on ###
df_model_matrix_final <- data.frame(cbind("fitness" = design_matrix$fitness,model_matrix_final))
str(df_model_matrix_final)
str(df_model_matrix_final, list.len=ncol(model_matrix_final))

# LM fit on full model with all the interaction terms -> 111 terms [[ 11 one-way (10 site + 1 host), 55 two-way terms (11C2 = 55), 45 two-way site by 1 host terms (10C2 * 1 = 45) ]]
fit_obj_lm <- lm(fitness ~., data = df_model_matrix_final)
summary(fit_obj_lm)
summary.aov(fit_obj_lm)
# to get the R^2
summary(fit_obj_lm)$adj.r.squared 

# Using step() function to choose the best fit model based on minimizing AIC
fit_obj_step<-step(fit_obj_lm)  
summary(fit_obj_step)
summary.aov(fit_obj_step)
summary(fit_obj_step)$adj.r.squared 
# compare this model's R^2 with the full model's R^2


# # Using stepAIC() function.. We get the same answer as above
# fit_obj_stepAIC<-stepAIC(fit_obj_lm)
# summary(fit_obj_stepAIC)
# summary.aov(fit_obj_stepAIC)
# summary(fit_obj_stepAIC)$adj.r.squared
# 
# # Using step() function with direction.. We get the same result as normal step function
fit_obj_step_direction<-step(fit_obj_lm, upper=~., scope=~., direction = c("both"))
summary(fit_obj_step_direction)
summary.aov(fit_obj_step_direction)
summary(fit_obj_step_direction)$adj.r.squared

influence.measures(fit_obj_step_direction)
residuals(fit_obj_step)


# checking that R^2, residuals and the sum of squares are consistent
SSTotal <- var( df_model_matrix_final$fitness ) * (nrow(df_model_matrix_final)-1)
SSE     <- sum( fit_obj_step$resid^2 )
anova_output_step <- anova(fit_obj_step)
SS <- anova_output_step[,2]
sum(SS)
sum_SS_3way <- sum(SS[48:76])
sum_SS_1_way <- sum(SS[1:8])

# to get all the terms present in the final model 
interaction_names_step <- row.names(anova(fit_obj_step))
interaction_names_step_matrix <- matrix(row.names(anova(fit_obj_step)))

## To extract formula of the best fit model using step() funtion
formula_step <- formula(fit_obj_step)
# To add one term in that formula 
formula_step_wSite7 <- update(formula_step, ~.+ genotype.V71)
fit_obj_step_wSite7 <-  lm(formula_step_wSite7,data = df_model_matrix_final)
summary.aov(fit_obj_step_wSite7)
summary(fit_obj_step_wSite7)$adj.r.squared 

##
#vif(fit_obj_step)

#######  Benjamini Hochberg method
## Writing a function to get final test decisions based on controlling fdr using Benjamini Hochberg method
fdr_decision_cal<-function(p_val,alpha){
  p_val_ordered<-sort(p_val)
  fdr_decision<-rep(0,length(p_val_ordered))
  flag_atleast_1_discovery<-0
  for (i in length(p_val_ordered):1){
    if(p_val_ordered[i]<=((i/length(p_val_ordered))*alpha)){
      i_max<-i
      flag_atleast_1_discovery<-1
      fdr_decision[order(p_val)[1:i_max]]<-1
      break()
    }
  }
  return(list(fdr_decision,flag_atleast_1_discovery))
}

# Extracting p-values of all the terms present in the final model
p_val_final_model <- summary.aov(fit_obj_step)[[1]][["Pr(>F)"]] 
p_val_final_model <- p_val_final_model[-c(length(p_val_final_model))] # removing the last entry which is NA for the Residuals in the output of summary.aov
str(p_val_final_model)

## Correcting at 0.1 significance level (normally we use 0.1 after fdr correction)
final_decison_vec<-fdr_decision_cal(p_val = p_val_final_model,alpha = 0.1)[[1]]

## Extracting the genotypes with significance differences in host
fdr_passed_final_terms <- interaction_names_step_matrix[which(final_decison_vec==1)]

# Recording output
# write.table(fdr_passed_final_terms, "FDR passed significant terms in AIC reduced model.csv", sep="\t")
# write.table(interaction_names_step_matrix, "significant terms in AIC reduced model.csv", sep="\t")

#######

## Checking VIF with only main effects

X_matrix_model <- model.matrix(object=~., data=X_matrix)
X_matrix_model <- X_matrix_model[,-c(1)]
X_matrix_wFitnness <- data.frame(cbind("fitness" = design_matrix$fitness,data.frame(X_matrix_model)))
str(X_matrix_wFitnness)
fit_main_lm <- lm(fitness ~., data = X_matrix_wFitnness)
summary(fit_main_lm)
summary.aov(fit_main_lm)
vif(fit_main_lm)

#
rcorr_matrix <- rcorr(as.matrix(X_matrix_model))
str(rcorr_matrix)
corrgram_matrix <- corrgram(data.frame(X_matrix_model))
str(corrgram_matrix)

#
rcorr_model_matrix_final <- rcorr(as.matrix(model_matrix_final))
str(rcorr_model_matrix_final)
corrgram_model_matrix_final <- corrgram(data.frame(model_matrix_final))
str(corrgram_model_matrix_final)
library("xlsx")
write.xlsx(corrgram_model_matrix_final)


# to check none of the two sites are perfectly coorelated in nthe design_matrix
design_site_host <- design_matrix[,-c(1,2)]
for(i in 1:(ncol(design_site_host)-1)){
  for(j in (i+1):ncol(design_site_host)){
    print(prod(design_site_host[,i],design_site_host[,j]))
  }
}

### From Type III SS in R in  http://md.psych.bio.uni-goettingen.de/mv/unit/lm_cat/lm_cat_unbal_ss_explained.html 

options(contrasts = c("contr.sum","contr.poly"))
model_after_drop1 <- drop1(fit_obj_lm, .~., test="F")
