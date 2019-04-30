####### Workshop with Statistical Testing
####### Presenter:  Qinlu (Claire) Wang M.S.
####### Email: qinlu.wang@nih.gov

######################################################################################
############################### 0. Data Preprocessing ################################
######################################################################################

##### Data Transformation #####
## 1st Method: Click from the "Files"
transformation_data <- data.frame(transformation)
transformation_data
## 2nd Method: Read with read.csv
transformation_data <- read.csv(file.choose()) # Choose the file "transformation.csv" #
transformation_data
# Or set the path and read the file name 
#setwd("/Users/Desktop/desktop/Seminar/data")
#transformation_data <- read.csv("transformation.csv")
#transformation_data

library(rcompanion)
plotNormalHistogram(transformation_data$Turbidity)
qqnorm(transformation_data$Turbidity, ylab="Sample Quantiles for Turbidity")
qqline(transformation_data$Turbidity, col="red")

## First, Square Root Transformation ##
transformation_data$square_root <- sqrt(transformation_data$Turbidity)
plotNormalHistogram(transformation_data$square_root)

## Second, Cube root transformation ##
transformation_data$cube_root = sign(transformation_data$Turbidity) * abs(transformation_data$Turbidity)^(1/3)
plotNormalHistogram(transformation_data$cube_root)

## Third, log transformation ##
transformation_data$log = log(transformation_data$Turbidity)
plotNormalHistogram(transformation_data$log)

## Fourth, square transformation ##
transformation_data$square = (transformation_data$Turbidity) ^ 2
plotNormalHistogram(transformation_data$square)


######################################################################################
################################## 1. Correlations ###################################
######################################################################################

# Link: http://rcompanion.org/handbook/I_10.html

install.packages("psych")
install.packages("PerformanceAnalytics")
install.packages("ggplot2")
install.packages("rcompanion")

######## Step 1: Read Data ########
## 1st Method: Click from the "Files"
correlation_data <- data.frame(correlation_data)
correlation_data
## 2nd Method: Read with read.csv
correlation_data <- read.csv(file.choose()) # Choose the file "correlation.csv" #
correlation_data
# Or set the path and read the file name 
#setwd("/Users/Desktop/desktop/Seminar/data")
#correlation_data <- read.csv("correlation.csv")
#correlation_data

######## Step 2: Visualize correlated variables ########
pairs(data=correlation_data, ~ Ozone + Solar.R + Wind + Temp)

######## Step 3: Correlation ########
library(psych)
# correlation matrix, p-value matrix and confidence interval
print(corr.test(correlation_data, 
                use    = "pairwise",
                method = "pearson",
                adjust = "none"), short=FALSE)

######## Step 4 (Optional): Put all results in one plot ########
library(PerformanceAnalytics)
chart.Correlation(correlation_data,
                  method="pearson",
                  histogram=TRUE,
                  pch=16)
# In the output, the numbers represent the correlation coefficients. 
# The stars represent the p-value of the correlation:
#   *   p < 0.05
#  **   p < 0.01
# ***   p < 0.001

######################################################################################
################################ 1. Two-Sample T-Test ################################
######################################################################################

# Link: http://www.sthda.com/english/wiki/normality-test-in-r

######## Step 1: Read Data ########
## 1st Method: Click from the "Files"
two_sample_t_test <- data.frame(two_sample_t_test)
two_sample_t_test
## 2nd Method: Read with read.csv
two_sample_t_test <- read.csv(file.choose()) # Choose the file "two-sample t-test.csv" #
two_sample_t_test
# Or set the path and read the file name 
#setwd("/Users/Desktop/desktop/Seminar/data")
#two_sample_t_test <- read.csv("two-sample t-test.csv")
#two_sample_t_test

######## Step 2: Check Normality ########
# Install package "dplyr" for data manipulation #
install.packages("dplyr")
library(dplyr)

### QQ-plot ###
color = ifelse(two_sample_t_test$Gender == "Male", "Blue", "Red")
qqnorm(two_sample_t_test$Y, col = color)
qqline(two_sample_t_test$Y)
legend(-1.6, 75, legend=c("Male", "Female"), col=c("Blue", "Red"), pch = c(1, 1))
# -> As all the points fall approximately along this reference line, we can assume normality.

### Normalty Test ###
# Visual inspection, described in the previous section, it’s possible to use a significance test comparing the sample distribution to a normal one in order to ascertain 
# whether data show or not a serious deviation from normality.
# There are several methods for normality test such as Kolmogorov-Smirnov (K-S) normality test and Shapiro-Wilk’s test.
# **The null hypothesis of these tests is that “sample distribution is normal”. If the test is significant, the distribution is non-normal.

# Shapiro-Wilk’s method is widely recommended for normality test and it provides better power than K-S. It is based on 
# the correlation between the data and the corresponding normal scores.
# Note that, normality test is sensitive to sample size. Small samples most often pass normality tests. 
# Therefore, it’s important to combine visual inspection and significance test in order to take the right decision.
shapiro.test(two_sample_t_test[two_sample_t_test$Gender == "Male", "Y"]) 
shapiro.test(two_sample_t_test[two_sample_t_test$Gender == "Female", "Y"]) 

# -> From the output, the p-value > 0.05 implying that the distribution of the data are not significantly 
#    different from normal distribution. In other words, we can assume the normality.

######## Step 3: Check same SDs ########
# We’ll use F-test to test for homogeneity in variances. This can be performed with the function var.test() as follow:
res.ftest <- var.test(Y ~ Gender, data = two_sample_t_test)
res.ftest

######## Step 4: Compute t-test ########
res <- t.test(Y ~ Gender, data = two_sample_t_test, var.equal = TRUE)
res
# -> Since the p-value of the test is larger than the significance level alpha = 0.05, we cannot conclude that a 
# significant difference exists between men and women.

######################################################################################
############################ 2. Unmatched One-Way ANOVA ##############################
######################################################################################

# Link: http://www.sthda.com/english/wiki/one-way-anova-test-in-r

######## Step 1: Read Data ########
## 1st Method: Click from the "Files"
unmatched_one_way_anova <- data.frame(unmatched_one_way_anova)
unmatched_one_way_anova
## 2nd Method: Read with read.csv
unmatched_one_way_anova <- read.csv(file.choose()) # Choose the file "unmatched one-way anova.csv" #
unmatched_one_way_anova
# Or set the path and read the file name 
#setwd("/Users/Desktop/desktop/Seminar/data")
#unmatched_one_way_anova <- read.csv("unmatched one-way anova.csv")
#unmatched_one_way_anova

######## Step 2: Check Normality ########
shapiro.test(unmatched_one_way_anova[unmatched_one_way_anova$Group == "Control", "Y"]) 
shapiro.test(unmatched_one_way_anova[unmatched_one_way_anova$Group == "Treated", "Y"]) 
shapiro.test(unmatched_one_way_anova[unmatched_one_way_anova$Group == "Treated+Antagonist", "Y"]) 

# Box-plots
install.packages("ggpubr")
library(ggpubr)
ggboxplot(unmatched_one_way_anova, x = "Group", y = "Y", 
          color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Control", "Treated", "Treated+Antagonist"),
          ylab = "Y", xlab = "Groups")
ggline(unmatched_one_way_anova, x = "Group", y = "Y", 
       add = c("mean_se", "jitter"), 
       order = c("Control", "Treated", "Treated+Antagonist"),
       ylab = "Y", xlab = "Groups")

######## Step 3: One-Way ANOVA ########
# Compute the analysis of variance
res.aov <- aov(Y ~ Group, data = unmatched_one_way_anova)
summary(res.aov)

# Non-parametric method: Kruskal-Wallis Test
kruskal.test(Y ~ Group, data = unmatched_one_way_anova)
# -> As the p-value is less than the significance level 0.05, we can conclude that there are 
#    significant differences between the treatment groups.

######################################################################################
######################## 3. Repeated measures one-way ANOVA ##########################
######################################################################################

# Link: http://rcompanion.org/handbook/I_09.html

######## Step 1: Read Data ########
## 1st Method: Click from the "Files"
repeated_one_way_anova <- data.frame(repeated_measures_one_way_anova)
repeated_one_way_anova
## 2nd Method: Read with read.csv
repeated_one_way_anova <- read.csv(file.choose()) # Choose the file "repeated measures one-way anova.csv" #
repeated_one_way_anova
# Or set the path and read the file name 
#setwd("/Users/Desktop/desktop/Seminar/data")
#repeated_one_way_anova <- read.csv("repeated measures one-way anova.csv")
#repeated_one_way_anova

######## Step 2: Check Normality ########
shapiro.test(repeated_one_way_anova$Control) 
shapiro.test(repeated_one_way_anova$T1) 
shapiro.test(repeated_one_way_anova$T2) 
shapiro.test(repeated_one_way_anova$T3) 

######## Step 3: Repeated measures one-way ANOVA ########
multmodel=lm(cbind(Control,T1,T2,T3) ~ 1, data = repeated_one_way_anova) # we column bind the 4 response columns together, and use only the intercept as our predictor variable.
Trials=factor(c("Control","T1","T2","T3"), ordered=F) # create a factor called “Trials” with labels control, T1, T2 and T3.
Trials

install.packages("car")
library(car)
model1=Anova(multmodel,idata=data.frame(Trials),idesign=~Trials,type="III") # We name the model 'multmodel', 
#then the ‘idata’ parameter is for the repeated-measures part of the data, the ‘idesign’ is where
#you specify the repeated part of the design. This is a one-way within design, so we have only the ‘Trials’ variable.
summary(model1,multivariate=F) 
# omit the ‘multivariate=F’ part, and you will get the MANOVA results as well.

######################################################################################
############################ 4. Ordinary Two-way ANOVA ###############################
######################################################################################

######## Step 1: Read Data ########
## 1st Method: Click from the "Files"
two_way_anova <- data.frame(two_way_anova)
two_way_anova
## 2nd Method: Read with read.csv
two_way_anova <- read.csv(file.choose()) # Choose the file "two-way anova.csv" #
two_way_anova
# Or set the path and read the file name 
#setwd("/Users/Desktop/desktop/Seminar/data")
#two_way_anova <- read.csv("two-way anova.csv")
#two_way_anova

#Generate frequency tables:
table(two_way_anova$Treatment, two_way_anova$Cell)

#Visualization
install.packages("ggpubr")
library("ggpubr")
ggboxplot(two_way_anova, x = "Cell", y = "Y", color = "Treatment", palette = c("#00AFBB", "#E7B800"))

######## Step 2: Check Normality ########
shapiro.test(two_way_anova[two_way_anova$Treatment == "Serum Starved" & two_way_anova$Cell == "Wild-Type cells", "Y"]) 
shapiro.test(two_way_anova[two_way_anova$Treatment == "Serum Starved" & two_way_anova$Cell == "GPP5 cell line", "Y"]) 
shapiro.test(two_way_anova[two_way_anova$Treatment == "Normal Culture" & two_way_anova$Cell == "Wild-Type cells", "Y"]) 
shapiro.test(two_way_anova[two_way_anova$Treatment == "Normal Culture" & two_way_anova$Cell == "GPP5 cell line", "Y"]) 

######## Step 3: Ordinary Two-way ANOVA ########
res.aov2 <- aov(Y ~ Treatment + Cell + Treatment:Cell, data = two_way_anova)
summary(res.aov2)
#-> It can be seen that the two main effects (supp and dose) are statistically significant, as well as their interaction.

## Run one-way ANOVA/two-sample t-test for each cell groups
# Wild-Type cells
res_wild_type <- t.test(Y ~ Treatment, data = two_way_anova[two_way_anova$Cell == "Wild-Type cells",], var.equal = TRUE)
res_wild_type
# GPP5 cell line
res_gpp5_cell<- t.test(Y ~ Treatment, data = two_way_anova[two_way_anova$Cell == "GPP5 cell line",], var.equal = TRUE)
res_gpp5_cell

######## Step 4: Check Assumptions ########
##Homogeneity of variances
plot(res.aov2, 1)
#-> Points 32 and 23 are detected as outliers, which can severely affect normality and homogeneity of variance. 
#   It can be useful to remove outliers to meet the test assumptions.

##Levene’s test to check the homogeneity of variances.
library(car)
leveneTest(Y ~ Treatment*Cell, data = two_way_anova)
#-> From the output above we can see that the p-value is not less than the significance level of 0.05. This means 
#   that there is no evidence to suggest that the variance across groups is statistically significantly different. 
#   Therefore, we can assume the homogeneity of variances in the different treatment groups.

##Normality Assumption
plot(res.aov2, 2)
#-> As all the points fall approximately along this reference line, we can assume normality.
# Extract the residuals
aov_residuals <- residuals(object = res.aov2)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals)

######## Step 5: Multiple Comparisons ########
TukeyHSD(res.aov2, which = "dose")
TukeyHSD(res.aov2)