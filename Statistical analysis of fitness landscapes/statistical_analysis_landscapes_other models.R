# LM on combined fitness landscape data but with using less complex models instead of full-factorial model with GxG, GxGxE, GxE

#write.table(mydata, "c:/mydata.txt", sep="\t")

# LM fit with only main interactions of host and mutations -> total of 11 terms
X_site_matrix_wHost <- lapply(design_matrix[,-c(1,2)], factor)  # removing intercept, genotypeID and fitness column
str(X_site_matrix_wHost) 
model_matrix_1way_final <- model.matrix(object = ~., data = X_site_matrix_wHost)
model_matrix_1way_final <- model_matrix_1way_final[,-c(1)]
head(model_matrix_1way_final)
model_matrix_1way_final <- data.frame(model_matrix_1way_final)
df_model_matrix_1way_final <- data.frame(cbind("fitness" = design_matrix$fitness,model_matrix_1way_final))
fit_obj_lm_1way <- lm(fitness ~., data = df_model_matrix_1way_final)
summary(fit_obj_lm_1way)
summary.aov(fit_obj_lm_1way)
# to get the R^2
summary(fit_obj_lm_1way)$adj.r.squared 


# LM with only main interactions, and GbyG, and GbyE
model_matrix_2way_final <- data.frame(two_way_model_matrix)
df_model_matrix_2way_final <- data.frame(cbind("fitness" = design_matrix$fitness,model_matrix_2way_final))
fit_obj_lm_2way <- lm(fitness ~., data = df_model_matrix_2way_final)
summary(fit_obj_lm_2way)
summary.aov(fit_obj_lm_2way)
# to get the R^2
summary(fit_obj_lm_2way)$adj.r.squared 

fit_obj_step_2way<-step(fit_obj_lm_2way) 
summary(fit_obj_step_2way)$adj.r.squared 

SSTotal_2way <- var( df_model_matrix_2way_final$fitness ) * (nrow(df_model_matrix_2way_final)-1)
SSE_2way    <- sum( fit_obj_lm_2way$resid^2 )
anova_output_step_2way <- anova(fit_obj_lm_2way)
SS_2way <- anova_output_step[,2]
sum(SS_2way)
#sum_SS_3way <- sum(SS[48:76])
#sum_SS_1_way <- sum(SS[1:8])



############ For analysis of individual landscapes ###############

#data$design.matrix.genotype.wt <- wt
#data$design.matrix.genotype.malT <- malT

# Creating design matrix for wt and malT from combined design_matrix

design_matrix_wt <- design_matrix[1:495,]
design_matrix_wt <- design_matrix_wt[,-c(2,13)]  # removing genotypeID and fitness columns
design_matrix_wt$genotype<-as.factor(design_matrix_wt$genotype)

design_matrix_malT <- design_matrix[496:803,]
design_matrix_malT <- design_matrix_malT[,-c(2,13)]

### for wt ###

X_matrix_wt <- lapply(design_matrix_wt[,-c(1)], factor)  # convert all columns into factors except fitness
two_way_model_matrix_wt <- model.matrix(object = ~.^2, data = X_matrix_wt)
two_way_model_matrix_wt <- two_way_model_matrix_wt[,-c(1)]
str(two_way_model_matrix_wt)
head(two_way_model_matrix_wt)
ncol(two_way_model_matrix_wt)

model_matrix_wt_final <- data.frame(two_way_model_matrix_wt)
str(model_matrix_wt_final, list.len=ncol(model_matrix_wt_final))
ncol(model_matrix_wt_final)

# FINAL data frame for wt #
df_model_matrix_wt_final <- data.frame(cbind("fitness" = design_matrix_wt$fitness,model_matrix_wt_final))
str(df_model_matrix_wt_final)
str(df_model_matrix_wt_final, list.len=ncol(model_matrix_wt_final))

fit_obj_lm_wt <- lm(fitness ~., data = df_model_matrix_wt_final)
summary(fit_obj_lm_wt)
summary.aov(fit_obj_lm_wt)
# to get the R^2
summary(fit_obj_lm_wt)$adj.r.squared 

interaction_names_wt <- row.names(anova(fit_obj_lm_wt))
interaction_names_matrix_wt <- matrix(row.names(anova(fit_obj_lm_wt)))

##### THE STEP function only removes 2 terms so not using it in final paper ####
# Using step() function to choose the best fit model based on minimizing AIC ##### THE STEP function only removes 2 terms so not using it ####
fit_obj_step_wt<-step(fit_obj_lm_wt)  
summary(fit_obj_step_wt)
summary.aov(fit_obj_step_wt)         
summary(fit_obj_step_wt)$adj.r.squared 
# compare this model's R^2 with the full model's R^2
# vif(fit_obj_step_wt)
# to get all the terms present in the final model 
interaction_names_step_wt <- row.names(anova(fit_obj_step_wt))
interaction_names_step_matrix_wt <- matrix(row.names(anova(fit_obj_step_wt)))

#######  Benjamini Hochberg method
## Writing a function to get final test decisions based on controlling fdr using Benjamini Hochberg method (I wrote this function down to use in one of my other projects and it has been verified multiple times!!)
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
p_val_final_model_wt <- summary.aov(fit_obj_lm_wt)[[1]][["Pr(>F)"]] 
p_val_final_model_wt <- p_val_final_model_wt[-c(length(p_val_final_model))] # removing the last entry which is NA for the Residuals in the output of summary.aov
str(p_val_final_model_wt)

## Correcting at 0.1 significance level (normally we use 0.1 after fdr correction)
final_decison_vec_wt<-fdr_decision_cal(p_val = p_val_final_model_wt,alpha = 0.1)[[1]]

## Extracting the genotypes with significance differences in host
fdr_passed_final_terms_wt <- interaction_names_matrix_wt[which(final_decison_vec_wt==1)]

write.table(fdr_passed_final_terms_wt, "FDR passed significant terms in WT landscape regression.csv", sep="\t")

# checking that R^2, residuals and the sum of squares are consistent
SSTotal_wt <- var( df_model_matrix_wt_final$fitness ) * (nrow(df_model_matrix_wt_final)-1)
SSE_wt     <- sum( fit_obj_lm_wt$resid^2 )
anova_output_wt <- anova(fit_obj_lm_wt)
SS_wt <- anova_output_wt[,2]
sum(SS_wt)
sum_SS_2way_wt <- sum(SS_wt[11:55])
sum_SS_1way_wt <- sum(SS_wt[1:10])


### for malT ###

X_matrix_malT <- lapply(design_matrix_malT[,-c(1)], factor)  # convert all columns into factors except fitness
two_way_model_matrix_malT <- model.matrix(object = ~.^2, data = X_matrix_malT)
two_way_model_matrix_malT <- two_way_model_matrix_malT[,-c(1)]
str(two_way_model_matrix_malT)
head(two_way_model_matrix_malT)
ncol(two_way_model_matrix_malT)

model_matrix_malT_final <- data.frame(two_way_model_matrix_malT)
str(model_matrix_malT_final, list.len=ncol(model_matrix_malT_final))
ncol(model_matrix_malT_final)

# FINAL data frame for malT #
df_model_matrix_malT_final <- data.frame(cbind("fitness" = design_matrix_malT$fitness,model_matrix_malT_final))
str(df_model_matrix_malT_final)
str(df_model_matrix_malT_final, list.len=ncol(model_matrix_malT_final))

fit_obj_lm_malT <- lm(fitness ~., data = df_model_matrix_malT_final)
summary(fit_obj_lm_malT)
summary.aov(fit_obj_lm_malT)
# to get the R^2
summary(fit_obj_lm_malT)$adj.r.squared 

interaction_names_malT <- row.names(anova(fit_obj_lm_malT))
interaction_names_matrix_malT <- matrix(row.names(anova(fit_obj_lm_malT)))

##### THE STEP function only removes 2 terms so not using it in final paper ####
# Using step() function to choose the best fit model based on minimizing AIC ##### THE STEP function only removes 2 terms so not using it ####
fit_obj_step_malT<-step(fit_obj_lm_malT)  
summary(fit_obj_step_malT)
summary.aov(fit_obj_step_malT)
summary(fit_obj_step_malT)$adj.r.squared 
# compare this model's R^2 with the full model's R^2
# vif(fit_obj_step_malT)
# to get all the terms present in the final model 
interaction_names_step_malT <- row.names(anova(fit_obj_step_malT))
interaction_names_step_matrix_malT <- matrix(row.names(anova(fit_obj_step_malT)))

#######  Benjamini Hochberg method
## Writing a function to get final test decisions based on controlling fdr using Benjamini Hochberg method (I wrote this function down to use in one of my other projects and it has been verified multiple times!!)
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
p_val_final_model_malT <- summary.aov(fit_obj_lm_malT)[[1]][["Pr(>F)"]] 
p_val_final_model_malT <- p_val_final_model_malT[-c(length(p_val_final_model))] # removing the last entry which is NA for the Residuals in the output of summary.aov
str(p_val_final_model_malT)

## Correcting at 0.1 significance level (normally we use 0.1 after fdr correction)
final_decison_vec_malT<-fdr_decision_cal(p_val = p_val_final_model_malT,alpha = 0.1)[[1]]

## Extracting the genotypes with significance differences in host
fdr_passed_final_terms_malT <- interaction_names_matrix_malT[which(final_decison_vec_malT==1)]

write.table(fdr_passed_final_terms_malT, "FDR passed significant terms in malT landscape regression.csv", sep="\t")

# checking that R^2, residuals and the sum of squares are consistent
SSTotal_malT <- var( df_model_matrix_malT_final$fitness ) * (nrow(df_model_matrix_malT_final)-1)
SSE_malT     <- sum( fit_obj_lm_malT$resid^2 )
anova_output_malT <- anova(fit_obj_lm_malT)
SS_malT <- anova_output_malT[,2]
sum(SS_malT)
sum_SS_2way_malT <- sum(SS_malT[11:55])
sum_SS_1way_malT <- sum(SS_malT[1:10])


######

AvsAB_606 <- c(2.937604532, 1.845510535,0.984737269)
t.test(AvsAB_606,var.equal = FALSE, conf.level = 0.95)
wilcox.test(AvsAB_606)

AvsAB_EcC4 <- c(-1.530188541, -1.028707979, -1.349502281) 
t.test(AvsAB_EcC4,var.equal = FALSE, conf.level = 0.95)


c126vsA_606 <- c(11.86839376, 11.46645951, 11.7377501)
t.test(c126vsA_606,var.equal = FALSE, conf.level = 0.95)
c126vsA_EcC4 <- c(-0.09680756,-0.052056362,-0.108544833)
t.test(c126vsA_EcC4,var.equal = FALSE, conf.level = 0.95)

EcC4vs606_wc126 <- c(0.040660043, 0.006570665, 0.051165343)
t.test(EcC4vs606_wc126,var.equal = FALSE, conf.level = 0.95)
t.test(EcC4vs606_wc126,var.equal = TRUE, conf.level = 0.95)
t.test(EcC4vs606_wc126,var.equal = FALSE, alternative = c("g"), conf.level = 0.95)
