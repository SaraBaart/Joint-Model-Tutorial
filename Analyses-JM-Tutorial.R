###################################
##   Load necessary libraries    ##
###################################
library(JMbayes)
library(nlme)
library(latticeExtra)
library(splines)


###############################
##     Data preparation      ##
###############################

# The data set needs to be in the long format

# Variables in the data set:
# id    = Patient identifier
# y     = Longitudinal outcome
# y2    = Second longitudinal outcome
# time  = Time of the longitudinal measurement
# time2 = Variable "time" shifted one measurement. 
#         Necessary for interval censored model
# Time  = Time of the event
# event = Event identifier
# group = Binary covariate

# For the interval censored model, time cannot start on exactly 0
data$time <- ifelse(data$time == 0, data$time+0.01, 
                    data$time)

# Make the data in the short format 
data.id <-  data[!rev(duplicated(rev(data$id))),] 

# Get the event rate for the  data set
table(data.id$event)[2] / nrow(data.id)

# Plot the data
xyplot(y ~ time | event, group=id, data = data, type = "l")



##############################################
##     Fit the time-dependent Cox model     ##
##############################################
##       	 Use y as biomarker   	        ##
##############################################


# TD-Cox model
TD.Cox <- (coxph((Surv(time, time2 , event) ~ y + group + 
                    cluster(id)), data = data))

summary(TD.Cox)


########################################
##      Fit the Basic Joint Model     ##
########################################
##        Use y as biomarker          ##
########################################

# Fit the survival model
Surv <- coxph(Surv(Time, event) ~ group, 
              data = data.id, x = TRUE, model = TRUE)

# Fit the mixed model
multMixedFit <- mvglmer(list(y ~ ns(time, knots = c(2,10)) + group  +
                               (ns(time, knots = c(2,10)) | id)), 
                        data = data, families = list(gaussian))

# Fit the joint model
JM1 <- mvJointModelBayes(multMixedFit, Surv, timeVar = "time")

# Inspect the traceplots 
plot(JM1)
# Obtain the results
summary(JM1)
# Obtain the HRs
exp(summary(JM1)$Survival)[,c(1,4,5)]



########################################
##      Fit the Basic Joint Model     ##
########################################
##          Use y as biomarker        ##
##      Fit the model with IC data    ##
########################################

# Fit the survival model with IC data
SurvInt <- survreg(Surv(time, time2, event,  type = "interval") ~ group, 
                   data = data.id, x = TRUE, model = TRUE)
summary(SurvInt)

# Fit the mixed model
multMixedFit <- mvglmer(list(y ~ ns(time, knots = c(2,10)) + group  +
                               (ns(time, knots = c(2,10)) | id)), 
                        data = data, families = list(gaussian))

# Fit the joint model
JM1.IC <- mvJointModelBayes(multMixedFit, SurvInt, timeVar = "time")

# Inspect the traceplots 
plot(JM1.IC)
# Obtain the results
summary(JM1.IC)
# Obtain the HRs
exp(summary(JM1.IC)$Survival)[,c(1,4,5)]



############################################
##            Fit the JM model            ##
##  Use slope as additional association   ##
############################################
##          Use y as biomarker            ##
############################################

Forms <- list("y" = "value",
              "y" = list(fixed = ~ 0 + dns(time, knots = c(2,10)), 
                         indFixed = c(2:4) , 
                         random = ~ 0 + dns(time, knots = c(2,10)),
                         indRandom = 2:4, name = "slope"))
#Fit the second joint model
JM2 <- update(JM1, Formulas = Forms)

# Inspect the traceplots
plot(JM2)
# Obtain the results
summary(JM2)
# Obtain the HRs
exp(summary(JM2)$Survival)[,c(1,4,5)]



############################################
##            Fit the JM model            ##
##           Multimarker Model            ##
############################################
##              Use y and y2              ##
############################################

# Fit the mixed model for two markers
multMixedFit2 <- mvglmer(list(y ~ ns(time, knots = c(2,10)) + group + 
                                (ns(time, knots = c(2,10)) | id),
                              y2 ~ ns(time, knots = c(2,10)) + group + 
                                (ns(time, knots = c(2,10)) | id)), 
                         data = data, families = list(gaussian, gaussian))

# Fit the joint model
JM3 <- mvJointModelBayes(multMixedFit2, Surv, timeVar = "time")

# Inspect the traceplots
plot(JM3)
# Obtain the results
summary(JM3)
# Obtain the HRs
exp(summary(JM3)$Survival)[,c(1,4,5)]




###########################################
##      Fit the dynamic predictions      ##
###########################################

# Make a data set for a specific patient A
NDA <- data[data$id == 28,] 

# Estimate survival probabilities for different time points
survPredsA <- vector("list", nrow(NDA))
for (i in 1:nrow(NDA)){
  survPredsA[[i]] <- survfitJM(JM1, newdata = NDA[1:i,], idVar = "id")
}

# Plot the graph on four different time points
par(mfrow=c(2,2))
for (i in c(1,3,5,6)) {
  plot(survPredsA[[i]], estimator = "mean", conf.int = TRUE, 
       fill.area = TRUE, col.area = "lightgrey", col.abline = "black", 
       pch = 21, add.last.time.axis.tick = FALSE, include.y = TRUE, 
       main = NULL, ylab = "", xlim = c(0,25))
}
