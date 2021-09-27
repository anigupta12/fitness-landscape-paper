library(dplyr)
library(ggplot2)
library(writexl)
library(car)
library(usdm)
require(graphics)
library(Hmisc)
library("corrgram")
library(tidyr)
library(lme4)
library(lmerTest)

neutral1_data <- read.csv(file = "min_neutral1.csv", header = TRUE)
neutral2_data <- read.csv(file = "min_neutral2.csv", header = TRUE)

joint_dataset <- inner_join(neutral1_data, neutral2_data)

joint_long_dataset <- joint_dataset %>% 
  gather(N1R1:N2R4, key="Condition", value="Value")

joint_long_dataset$Time = substr(joint_long_dataset$Condition,1,2)
joint_long_dataset$Measure = substr(joint_long_dataset$Condition,3,4)


# test ####
m0 = lmer(Value ~ 1 + (1|Measure)+(1|ID), joint_long_dataset)
m1 = lmer(Value ~ Time + (1|Measure)+(1|ID), joint_long_dataset)# varying intercept

anova(m0,m1)
summary(m1)
m3 = lmer(Value ~ Time + (Time|Measure)+(1|ID), joint_long_dataset) # varying slope
lmer(Value ~ Time + (Time|Flask)+(Time|Panda), joint_long_dataset)