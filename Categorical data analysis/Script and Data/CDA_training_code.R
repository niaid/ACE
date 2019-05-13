#####################################
###    Categorial Data Analysis   ###
###       Date: May 6, 2019       ###
###        Version: 2             ###
###       Author: J. Gu           ###
#####################################

####################################################################################################
## Topic:  Contingency table, Odds ratio, Relative risk, 
## Independence test (Pearson, LR test, Cochran-Mantel-Haenszel, McNemar's, Ordinal data test),
####################################################################################################

## install packages
install.packages("conting")

########################### Contingency table ################################
##############################################################################

## When you have a complete contingency table, you can manually type in:
C <- matrix(byrow=T, nrow=2, ncol=2, c(14,9,1,21))
dimnames(C) <- list(ThalUse=c("yes","no"), Healing=c("yes","no"))

ctable <- as.table(C)
ctable

## When have data frame where each row represents one case, create table from "long" to "wide":
cases <- data.frame(
  Sex = c("M", "M", "F", "F", "F"), 
  EyeColor = c("brown", "blue", "brown", "brown", "brown")
)
cases

ctable <- table(cases)
ctable

## Create table from aggregate table:
table <- data.frame(
  smoker = c("Yes","No","Yes","No"),
  lung_cancer = c("Cases","Cases","Control","Control"),
  count = c(688,21,650,59)
)
table

ctable <- xtabs(count ~ smoker+lung_cancer, data=table)
ctable

# Reference : http://www.cookbook-r.com/Manipulating_data/Converting_between_data_frames_and_contingency_tables/

#################################### Odds Ratio ###################################
###################################################################################

odds_ratio <- function(n11,n12,n21,n22) {
  
  OR <- (n11*n22)/(n12*n21)
  se.logOR <- sqrt(sum(1/n11, 1/n12, 1/n21,1/n22))
  
  L <- log(OR) - 1.96*se.logOR
  U <- log(OR) + 1.96*se.logOR
  
  exp(L)
  exp(U)
  
  paste0("The odds ratio is ",
         round(OR,digits = 3),
         ". The 95% confidence interval for odds ratio is ",
         sprintf("[ %4.3f, %4.3f ]", exp(L), exp(U)),".")
  
}

odds_ratio(688,650,21,59)


################################## Relative Risk ##################################
###################################################################################

# rows for exposure and columns for disease status

relative_risk <- function(n11,n12,n21,n22) {
  
  p1 <- n11/sum(n11,n12)
  p2 <- n21/sum(n21,n22)
  
  se.logRR <- sqrt( (1-p1)/n11 + ((1-p2)/n21) )
  L <- log(p1/p2) - 1.96 * (se.logRR)
  U <- log(p1/p2) + 1.96 * (se.logRR)
  
  p1/p2
  exp(L) ; exp(U) 
  
  paste0("The relative risk is ",
         round(p1/p2,digits = 3),
         ". The 95% confidence interval for relative risk is ",
         sprintf("[ %4.3f, %4.3f ]", exp(L), exp(U)),".")
  
}

relative_risk(688,650,21,59)


######################################### Test of Independence ####################################
######################### Pearson Chi-square test ##############################
################################################################################

table <- data.frame(
  smoker = c("Yes","No","Yes","No"),
  lung_cancer = c("Cases","Cases","Control","Control"),
  count = c(688,21,650,59)
)
table

ctable <- xtabs(count ~ smoker + lung_cancer, data=table)
chisq.test(ctable, correct = FALSE)
#chisq.test(table(table) ) #for the Yates-corrected chi-square statistic 


######################### Likelihood ratio test ##############################
##############################################################################

G2 <- function(table) {
  X2 <- chisq.test(table)
  e <- X2$expected
  o <- X2$observed
  
  G2 <- 2*sum(o*log(o/e))
  df <- (nrow(table)-1)*(ncol(table)-1)
  G2.p_value <- 1-pchisq(G2, df)
  return(data.frame(G2=G2, df=df,p_value=G2.p_value))

# print result  
  paste0("The G2 statistic is ",
         round(G2,digits = 3),
         ", degree of freedom is ",
         df,
         ", p-value is ",
         round(G2.p_value,digits = 6), ", at significant level 0.05.")
  
}

G2(ctable)

######################### Fisher exact test ##################################
##############################################################################

tea <- matrix(c(3,1,1,3),ncol=2,byrow=TRUE)
dimnames(tea) <- list(PouringFirst=c("Milk","Tea"), GuessPouredFirst=c("Milk","Tea"))
tea
fisher.test(tea,alternative="greater")

######################### Cochran-Mantel-Haenszel test #######################
##############################################################################

Input <- ("
          Age       CVD      obese    Count
          <50       CVD      Obese     10
          <50       Non-CVD  Obese     90
          <50       CVD      Non-Obese 35
          <50       Non-CVD  Non-Obese 465
          >=50      CVD      Obese     36
          >=50      Non-CVD  Obese     164
          >=50      CVD      Non-Obese 25
          >=50      Non-CVD  Non-Obese 175
          ")

Data <- read.table(textConnection(Input),header=TRUE)


Table <- xtabs(Count ~ obese  + CVD + Age , 
               data=Data)
ftable(Table)
mantelhaen.test(Table)

########################### McNemar's test #####################################
################################################################################

voting <-
  matrix(c(175, 54, 16, 188),
         nrow = 2, byrow = F,
         dimnames = list("2004 Election" = c("Democrat", "Repulican"),
                         "2008 Election" = c("Democrat", "Repulican")))
voting

mcnemar.test(voting) # A closer approximation to the chi-squared distributin uses a continuity correction
mcnemar.test(voting, correct = F) # without continuity correction


########################### Ordinal data ########################################
#################################################################################

# data AOH is storage in counting package; library to use it
library(conting)

# Way 1: importing data from package
data(AOH)

# Way 2: importing data from project folder; 
AOH <- read.table("./data/AOH.txt", header = T) # you can change it to your path

# see basic info about this data set
names(AOH); dim(AOH); head(AOH)

# assign scores to ordinal variables
assign <- function(data) {
  data$hyp <- ifelse(data$hyp=="no",1,2)
  data$obe <- ifelse(data$obe=="low",1,ifelse(data$obe=="average",2,3))
  #data$alc <- recode(data$alc, "0" = 0, "1-2" = 1.5, "3-5"=4, "6+"=7)
  data$alc <- ifelse(data$alc=="0",0, data$alc)
  data$alc <- ifelse(data$alc=="1-2",1.5, data$alc)
  data$alc <- ifelse(data$alc=="3-5",4, data$alc)
  data$alc <- ifelse(data$alc=="6+",7, data$alc)
  return(data)
}

AOH_adj <- assign(AOH)

AOH_adj_1 <- AOH_adj[ rep(1:nrow(AOH_adj),AOH_adj$y), 2:4]

table(AOH_adj_1$obe, AOH_adj_1$alc)

# add column and row margins
addmargins(table(AOH_adj_1$obe, AOH_adj_1$alc))
addmargins(table(AOH_adj_1$obe, AOH_adj_1$hyp))
addmargins(table(AOH_adj_1$alc, AOH_adj_1$hyp))

# calculate pearson correlation
n <- nrow(AOH_adj_1)
r_obe_alc <- cor(AOH_adj_1$obe, AOH_adj_1$alc)
r_hyp_alc <- cor(AOH_adj_1$hyp, AOH_adj_1$alc)
r_obe_hyp <- cor(AOH_adj_1$obe, AOH_adj_1$hyp)

# calculate M2 for linear trend statistics
M2_obe_alc <- (n-1)*r_obe_alc^2
M2_hyp_alc <- (n-1)*r_hyp_alc^2
M2_obe_hyp <- (n-1)*r_obe_hyp^2

# calculate p-value
1-pchisq(M2_obe_alc,1)
1-pchisq(M2_hyp_alc,1)
1-pchisq(M2_obe_hyp,1)

# perform pearson chi-square test for result comparison
chisq.test(table(AOH_adj_1$obe, AOH_adj_1$alc))
chisq.test(table(AOH_adj_1$hyp, AOH_adj_1$alc))
chisq.test(table(AOH_adj_1$obe, AOH_adj_1$hyp))

# perform likelihood ratio test for result comparison
G2(table(AOH_adj_1$obe, AOH_adj_1$alc))
G2(table(AOH_adj_1$hyp, AOH_adj_1$alc))
G2(table(AOH_adj_1$obe, AOH_adj_1$hyp))


##### END #####
