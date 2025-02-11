# MIstol 2019: fruit set and total fruit production. 
# We carried out 3 treatments on S.mistol, here we use the average of these fruit set values

library(lme4)
library(Matrix)
library(ggplot2)
library(hrbrthemes)
library(multcomp)
library(mvtnorm)
library(survival)
library(TH.data)
library(MASS)
library(dplyr)
install.packages("DHARMa")
library(DHARMa)
library(glmmTMB)
library(ggeffects)

mis19<- read.table(file.choose(), header = TRUE)
head(mis19) # mistol2019_paper

# Check te interaction between diameter and land use condition
# Create the base plot with individual points
p3 <- ggplot(mis19, aes(x = uso, y = diam, fill = uso)) + 
  geom_boxplot(alpha = 1.2, outlier.shape = NA) +  # Add boxplot and hide outliers
  geom_jitter(width = 0.2, aes(color = uso), size = 2, alpha = 0.6) +  # Add individual points with jitter
  labs(
    title = "Diámetro Mistol",
    x = "Gradiente de uso",
    y = "Diámetro"
  ) +
  theme_classic() +
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 21,  # Use a filled point shape
    fill = "red", 
    size = 3
  )

# Reorder the x-axis items
c1 <- p3 + scale_x_discrete(limits = c("CSF", "SF", "CS", "OF"))

# Display the plot
c1

# It seam that tree diameter is bigger in CSF (conserve secundary forest) 
# We test this using a model

mis19 %>%
  group_by(parcela) %>%
  summarise(mean_diam = mean(diam), sd_diam = sd(diam))
# it looks like there is variability in both the mean diameter (mean_diam) 
#and the standard deviation (sd_diam) across parcela

mis19 %>%
  group_by(parcela) %>%
  summarise(n = n())
# it seems that each parcela has a reasonable number of observations 
#(mostly 5, with one having 4 and another having 6)

mod1<- lm(diam~uso,data=mis19) 
summary(mod1)

mod0<- lm(diam~1,data=mis19)
anova(mod0,mod1) # el gradiente de uso es sig. 

#I cannot use diameter as a covariate because there is no independence. 
#However, I can use it in an interaction

Cont.1<-glht(mod1,linfct=mcp(uso ="Tukey")) 
summary(Cont.1)

#Linear Hypotheses:
#Estimate Std. Error z value Pr(>|z|)   
#CSF - CS == 0    7.852     12.097   0.649  0.91585   
#OF - CS == 0     4.867     11.886   0.409  0.97683   
#SF - CS == 0   -37.318     12.335  -3.025  0.01373 * 
#OF - CSF == 0   -2.986     12.097  -0.247  0.99472   
#SF - CSF == 0  -45.170     12.538  -3.603  0.00180 **
#SF - OF == 0   -42.185     12.335  -3.420  0.00332 **

diagnostic_plots <- function(model) {
  par(mfrow = c(2, 2))  
  plot(model)           
  par(mfrow = c(1, 1))  
}

#model adjustment
diagnostic_plots(mod1)

#Residuals vs Fitted: detect systematic patterns in residuals
#There seems to be more disperssion in the higher values of the fitted values, 
#suggesting heteroscedasticity (non-constant variance of the residuals).
#Normal Q-Q:verify whether the residuals follow a normal distribution 
#There appears to be some deviation at the extreme values, suggesting that the residuals may not be completely normalThere appears to be some deviation at the extreme values, 
#suggesting that the residuals may not be completely normal
#Scale-Location: evaluate homoscedasticity
# There some heteroscedasticity
#Residuals vs. Factor Levels: assess whether the residuals have a homogeneous distribution at different levels of the categorical factors.
#No clear trend is observed, which is good, but there are some points with extreme values, which could be influencing the model.
library(performance)
# Check for multicollinearity
check_collinearity(mod1) # NULL

# Check for influential observations
check_outliers(mod1) # no outliers

# Check model assumptions
check_heteroscedasticity(mod1)
# Error variance appears to be homoscedastic (p = 0.550).
check_normality(mod1)
# residuals appear as normally distributed (p = 0.491).

###########################################################
#Check the interaction between conespecificos and land use condition
p3 <- ggplot(mis19, aes(x = uso, y = cones, fill = uso)) + 
  geom_boxplot(alpha = 1.2, outlier.shape = NA) +  # Add boxplot and hide outliers
  geom_jitter(width = 0.2, aes(color = uso), size = 2, alpha = 0.6) +  # Add individual points with jitter
  labs(
    title = "Diámetro Mistol",
    x = "Gradiente de uso",
    y = "Conespecificos"
  ) +
  theme_classic() +
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 21,  # Use a filled point shape
    fill = "red", 
    size = 3
  )
# Reorder the x-axis items
c1 <- p3 + scale_x_discrete(limits = c("CSF", "SF", "CS", "OF"))
# Display the plot
c1
#############
## We check with a linear model 
c1<-glmer(cones~uso +(1|parcela), data=mis19,family="poisson", na.action=na.fail)
summary(c1)

check_overdispersion(c1) # No overdispersion detected.

c2<-glm(cones ~ uso, family = poisson, data=mis19)
summary(c2)

anova(c1,c2,test="Chisq") # I can remove the random 


Cont.c2<-glht(c2,linfct=mcp(uso ="Tukey")) 
summary(Cont.c2)

#Estimate Std. Error z value Pr(>|z|)  
#CSF - CS == 0   0.3185     0.2078   1.532   0.4165  
#OF - CS == 0    0.5158     0.1998   2.582   0.0479 *
#SF - CS == 0    0.5596     0.1982   2.823   0.0242 *
#OF - CSF == 0   0.1974     0.1820   1.085   0.6977  
#SF - CSF == 0   0.2412     0.1802   1.338   0.5369  
#SF - OF == 0    0.0438     0.1709   0.256   0.9941  
# I can not use the conespecifics as a covariables, only interacting witg "uso"
##############################################################
##Fuit set (binomial distribution)
mis19$y<-cbind(mis19$fs,mis19$nfs)
# Select random effect

q2_glmer <- glmer(y ~ uso + (1 | parcela), data = mis19, family = binomial)
summary(q2_glmer)

check_overdispersion(q2_glmer)
#dispersion ratio = 0.557
#p-value = 0.016
#If the underdispersion is strong, sometimes the random effect is not necessary. 

q2_glm <- glm(y ~ uso, data = mis19, family = binomial)
anova(q2_glmer,q2_glm, test = "Chisq") # Random effect is no need

summary(q2_glm)
anova(q2_glm, test = "Chisq") # p= 0.0988 . land use is partially significant

mis19$uso <- as.factor(mis19$uso)
Cont.q2<-glht(q2_glm,linfct=mcp(uso ="Tukey")) 
summary(Cont.q2)

#Linear Hypotheses:
#Estimate Std. Error z value Pr(>|z|)
#CSF - CS == 0  0.45237    0.32226   1.404    0.495
#OF - CS == 0   0.53069    0.31763   1.671    0.337
#SF - CS == 0  -0.13555    0.36780  -0.369    0.983
#OF - CSF == 0  0.07832    0.27999   0.280    0.992
#SF - CSF == 0 -0.58792    0.33583  -1.751    0.295
#SF - OF == 0  -0.66624    0.33139  -2.010    0.182

#############################################
###########total fruit production: possion distribution #################
########################################
# Create the base plot with individual points
p3 <- ggplot(mis19, aes(x = uso, y = frutostot, fill = uso)) + 
  geom_boxplot(alpha = 1.2, outlier.shape = NA) +  # Add boxplot and hide outliers
  geom_jitter(width = 0.2, aes(color = uso), size = 2, alpha = 0.6) +  # Add individual points with jitter
  labs(
    title = "Frutos totales",
    x = "Gradiente de uso",
    y = "Diámetro"
  ) +
  theme_classic() +
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 21,  # Use a filled point shape
    fill = "red", 
    size = 3
  )

# Reorder the x-axis items
c1 <- p3 + scale_x_discrete(limits = c("CSF", "SF", "CS", "OF"))

# Display the plot
c1

#Select random effect 
ft1<-glmer(frutostot~uso +(1|parcela), data=mis19,family=poisson, na.action=na.fail)
summary(ft1)

check_overdispersion(ft1)
#dispersion ratio =   436.052
#Pearson's Chi-Squared = 23982.846
#p-value =   < 0.001
# We have a problem here


#### Negative Binomial Model ###############
library(glmmTMB)
ft1_nb <- glmmTMB(frutostot ~ uso + (1|parcela), data = mis19, family = nbinom2, na.action = na.fail)
summary(ft1_nb)

check_overdispersion(ft1_nb)
#dispersion ratio = 0.935
#p-value = 0.856

ft0_nb <- glmmTMB(frutostot ~ uso, data = mis19, family = nbinom2, na.action = na.fail)

check_overdispersion(ft0_nb)
#dispersion ratio = 1.014
#p-value = 0.808
anova(ft1_nb,ft0_nb) # p= 0.006087 ** The random effect is important 

# Select fixed effects

ft00_nb <- glmmTMB(frutostot ~ 1 + (1|parcela) , data = mis19, family = nbinom2, na.action = na.fail)

anova(ft1_nb,ft00_nb) #  Land use is significant  p=0.01414

##### Final model ######

ft1_nb <- glmmTMB(frutostot ~ uso + (1|parcela), data = mis19, family = nbinom2, na.action = na.fail)
summary(ft1_nb)


Cont.ft4<-glht(ft1_nb,linfct=mcp(uso="Tukey")) 
summary(Cont.ft4)
#Estimate Std. Error z value Pr(>|z|)   
#CSF - CS == 0  -0.7412     0.4238  -1.749  0.29843   
#OF - CS == 0    0.1060     0.4241   0.250  0.99452   
#SF - CS == 0   -1.4140     0.4241  -3.334  0.00442 **
#OF - CSF == 0   0.8472     0.4232   2.002  0.18721 
#SF - CSF == 0  -0.6728     0.4235  -1.589  0.38514   
#SF - OF == 0   -1.5200     0.4236  -3.588  0.00207 **

# Interaction with the diameter
ft1_nb_diam <- glmmTMB(frutostot ~ uso*diam + (1|parcela), data = mis19, family = nbinom2, na.action = na.fail)
summary(ft1_nb_diam)

check_overdispersion(ft1_nb_diam)

#dispersion ratio = 0.807
#p-value = 0.952

anova(ft1_nb_diam,ft1_nb, test = "Chisq") # The interaction is significant 
#p=  0.03969 *

### total fruit production vs diam
ggplot(mis19, aes(x=diam, y=frutostot, color=uso)) + 
  geom_point(size=4, alpha=0.7) +          # Slight transparency for better visualization
  geom_smooth(method="lm", se=FALSE, linetype="dashed", size=1) +  # Add trend lines
  labs(
    title = "Frutos Totales vs Diámetro",
    subtitle = "Comparación de frutos totales en función del diámetro por uso",
    x = "Diámetro (cm)",
    y = "Frutos Totales",
    color = "Uso"
  ) +                                      # Add titles and axis labels
  theme_ipsum(base_size = 15) +            # Increase base font size for better readability
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    plot.subtitle = element_text(hjust = 0.5) # Center the plot subtitle
  )


# Extract fitted values from the model
mis19$fitted <- predict(ft1_nb, type = "response")

# Plot observed vs. fitted values
ggplot(mis19, aes(x = fitted, y = frutostot)) +
  geom_point(alpha = 0.6, color = "blue") +  # Observed values
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Reference line
  labs(
    title = "Observed vs. Fitted Values",
    x = "Fitted Values",
    y = "Observed Values"
  ) +
  theme_minimal()

# Get predicted values for 'uso'
#This plot shows the predicted values of frutostot across the levels of uso, 
#while accounting for the random effect of parcela.
pred_uso <- ggpredict(ft1_nb, terms = "uso")

# Plot predicted values
plot(pred_uso) +
  labs(
    title = "Predicted frutostot by uso",
    x = "Uso",
    y = "Predicted frutostot"
  ) +
  theme_minimal()

# Get predicted values with confidence intervals
pred_data <- as.data.frame(pred_uso)

# Plot with ggplot2
ggplot(pred_data, aes(x = x, y = predicted)) +
  geom_point(size = 3, color = "blue") +  # Predicted values
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, color = "red") +  # Confidence intervals
  labs(
    title = "Predicted frutostot by uso",
    x = "Uso",
    y = "Predicted frutostot"
  ) +
  theme_minimal()

# Interaction with the conespecifics
ft1_nb_cones <- glmmTMB(frutostot ~ uso*cones + (1|parcela), data = mis19, family = nbinom2, na.action = na.fail)
summary(ft1_nb_cones)

check_overdispersion(ft1_nb_cones) #No overdispersion detected

anova(ft1_nb_cones,ft1_nb, test = "Chisq") # The interaction is not significant 

### total fruit production vs cones
ggplot(mis19, aes(x=cones, y=frutostot, color=uso)) + 
  geom_point(size=4, alpha=0.7) +          # Slight transparency for better visualization
  geom_smooth(method="lm", se=FALSE, linetype="dashed", size=1) +  # Add trend lines
  labs(
    title = "Frutos Totales vs Conespecifics",
    subtitle = "Comparación de frutos totales en función del numero de conespecíficos por uso",
    x = "Conespecíficos",
    y = "Frutos Totales",
    color = "Uso"
  ) +                                      # Add titles and axis labels
  theme_ipsum(base_size = 15) +            # Increase base font size for better readability
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    plot.subtitle = element_text(hjust = 0.5) # Center the plot subtitle
  )
