rm(list=ls())

library(Ecdat)
library(tidry)

data(RetSchool)

summary(RetSchool)

RetSchool <- drop_na(RetSchool)

View(RetSchool)

wage <- RetSchool$wage76
exp <- RetSchool$exp76
edu <- RetSchool$grade76

#linear regression

lm_levels <- lm(wage76 ~ exp76 + grade76, data = RetSchool)
lm_levels

# interpretation: one additional year of schooling / experience brings Beta more
# dollars an hour
# Opportunity cost of schooling for one more year vs. working that year is: beta2 - beta1

Op_cost <- 
  
# getting the sense of the data

plot(RetSchool$grade76, RetSchool$wage76, main = "Grade vs Wage in 1976",
     xlab = "Grade in 1976", ylab = "Wage in 1976", pch = 19, cex = 0.1)
lines(x, )

plot(RetSchool$exp76, RetSchool$wage76, main = "Experience vs Wage in 1976",
     xlab = "Grade in 1976", ylab = "Wage in 1976", pch = 19, cex = 0.1)

# taking the logarithm of the data ??? should I
# adding one to avoid -infty problems 

RetSchool$log_wage76 <- log(RetSchool$wage76+1)
RetSchool$log_exp76 <- log(RetSchool$exp76+1)
RetSchool$log_grade76 <- log(RetSchool$grade76+1)

lm_log <- lm(log_wage76 ~ log_exp76 + log_grade76, data = RetSchool)
lm_log

plot(RetSchool$log_exp76, RetSchool$wage76,xlab = "Grade in 1976", ylab = "Wage in 1976", 
     pch = 19, cex = 0.3)
