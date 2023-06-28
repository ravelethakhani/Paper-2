library(MASS)
library(e1071)
library(quantreg)
library(qgam)
library(devtools)
library(gefcom2017)
library(caret)
library(relaimpo)
library(glmnet)
library(evgam)
library(extremefit)
library(scoringRules)
library(fitdistrplus)

######################Hourly GHI#############################################
attach(Hourlydata)
head(Hourlydata)
#tail(data)

####Summary statistics for hours
summary(GHI)

####Standard deviation
sd(GHI)

####skewness
skewness(GHI)

####kurtosis
kurtosis(GHI)

#####ts and density plots for GHI
win.graph()
par(mfrow=c(1,2))
yGHI=ts(GHI)
plot(yGHI,main="(a) Plot of GHI",xlab = "Observation number",
     ylab="Hourly GHI (W/sqre m)",col="blue")
vGHI=density(GHI)
plot(vGHI,main="(b) Density of plot",xlab = "Hourly GHI (W/sqre m)",
     col="blue")

######Q-Q plot
qqnorm(GHI, pch = 1, frame = FALSE, main="Normal Q-Q plot")
qqline(GHI, col = "blue", lwd = 2)

######Bloxplot
boxplot(GHI, horizontal = TRUE, 
        main="(d) Box plot",xlab="Hourly GHI (W/sqre m)", 
        col = "blue")

##########################################################################################################
#####nonlinear trend
###GHI
win.graph()

y <- ts(GHI)
plot(y, ylab="Hourly GHI (W/sqre m)",col="blue",
     xlab="Observation number")
z = (smooth.spline(time(y), y))
z
lines(smooth.spline(time(y), y, spar = 0.2034927),col="red",lwd=2)
nlifit = fitted((smooth.spline(time(y), y, spar=0.2034927)))
nlifit <- round(nlifit,0)
write.table(nlifit,"~/nlifit.txt",sep="\t")

###GHI residuals 
r2=residuals((smooth.spline(time(y), y, spar = 0.2034927)))
r2 <- round(r2,0)
win.graph()
par(mfrow=c(1,2))
plot(r2,col="blue",ylab="Residuals (W/sq m)", 
     xlab="Observation number", lwd = 2)
plot(density(r2), col="blue", main="", xlab="Residuals (W/sq m)")

#r2
dpdfits = fitted((smooth.spline(time(y), y, spar = 0.2034927)))
excess <- round(GHI-dpdfits,0)
#excess

write.table(excess,"~/ResidualsHourly.txt",sep="\t")
#-----------------------------------------------------
# Percentiles
#------------------------------------------
win.graph()
plot(y,ylab="Hourly GHI (W/sqre m)",col="black",
     xlab="Observation number",ylim=c(0,2000))
# 95% u = 483
# 97% u = 523
# 99% u = 588
# 99.9% u = 720
# 99.99% u = 743

lines(smooth.spline(time(y), y+483, spar = 0.2034927),col="red",lwd=2)
lines(smooth.spline(time(y), y+523, spar = 0.2034927),col="yellow",lwd=2)
lines(smooth.spline(time(y), y+588, spar = 0.2034927),col="brown",lwd=2)
lines(smooth.spline(time(y), y+720, spar = 0.2034927),col="green",lwd=2)
lines(smooth.spline(time(y), y+743, spar = 0.2034927),col="blue",lwd=2)
legend("topright",col=c("black","red", "yellow", "brown", "green", "blue"),
       lty=1:3,lwd=2,legend=c("GHI", "95%", "97%", "99%", "99.9%", "99,99%"))

dpdfit95 = fitted((smooth.spline(time(y), y+483, spar = 0.2034927)))
dpdfit95 <- round(dpdfit95,0)

dpdfit97 = fitted((smooth.spline(time(y), y+523, spar = 0.2034927)))
dpdfit97 <- round(dpdfit97,0)

dpdfit99 = fitted((smooth.spline(time(y), y+588, spar = 0.2034927)))
dpdfit99 <- round(dpdfit99,0)

dpdfit999 = fitted((smooth.spline(time(y), y+720, spar = 0.2034927)))
dpdfit999 <- round(dpdfit999,0)

dpdfit9999 = fitted((smooth.spline(time(y), y+743, spar = 0.2034927)))
dpdfit9999 <- round(dpdfit9999,0)

write.table(dpdfit9999,"~/dpdfit9999.txt",sep="\t")

########################################################################################################
##### Calculates the pinball loss score for a given quantile.
library(gefcom2017)
####95%, 97%, 99%, 99,9% and 99,99%
pinball_loss <- function(tau, y, q) {
        pl_df <- data.frame(tau = tau,
                            y = y,
                            q = q)
        
        pl_df <- pl_df %>%
                mutate(L = ifelse(y>=q,
                                  tau/100 * (y-q),
                                  (1-tau/100) * (q-y)))
        
        return(pl_df)
}

tau= 99.99
y= GHI
q= dpdfit9999   # Flstm, Fsvr, Fffnn, Fconvex, Fqra
z = pinball_loss(tau, y, q)
z
write.table(z,"~/pinballFconvex.txt",sep="\t")
qloss =z$L
a=ts(qloss)
plot(a)
mean <- round(mean(qloss),4)
mean
qloss <- round(qloss,4)
write.table(qloss,"~/pinballBM95.txt",sep="\t")

#########################################################################################################
#### EXTREMEFIT OF Hourly GHI#
#########################################################################################################
attach(Hourlydata)
head(Hourlydata)

y=GHI
min(y)
max(y)

################Extreme conditional quantiles#########
library(extremefit)

## Estimate extreme quantile of DPED depending on index 

win.graph()
index <- 1:length(y) # index is 1,2,...,4083
y1 <- ts(GHI)
plot(y1 ,ylim=c(0,6600),lwd=2,lty=1, ylab = "Hourly GHI (W/sqre m)", 
     xlab = "Observation number")
TgridIndex <- seq(min(index), max(index),length = 4083) # here is on which points do you want to perform the estimation

#new.Tgrid <- Tgrid
HHIndex <- hill.ts(y, index, TgridIndex,h = 100, #I set the parameter h at 50 but you can change it freely and play with it
                   kernel = TruncGauss.kernel, CritVal = 3.4)

####### 0.95, 0.97, 0.99, 0.999, 0.9999-quantile
QuantIndex95<- as.numeric(predict(HHIndex,
                                newdata = 0.95, type = "quantile")$y)
QuantIndex97<- as.numeric(predict(HHIndex,
                                newdata = 0.97, type = "quantile")$y)
QuantIndex99<- as.numeric(predict(HHIndex,
                                newdata = 0.99, type = "quantile")$y)
QuantIndex999<- as.numeric(predict(HHIndex,
                                 newdata = 0.999, type = "quantile")$y)
QuantIndex9999<- as.numeric(predict(HHIndex,
                                 newdata = 0.9999, type = "quantile")$y)

######## Rounding
QuantIndex95 <- round(QuantIndex95,0)
QuantIndex97 <- round(QuantIndex97,0)
QuantIndex99 <- round(QuantIndex99,0)
QuantIndex999 <- round(QuantIndex999,0)
QuantIndex9999 <- round(QuantIndex9999,0)

######## Plotting lines
lines(TgridIndex,QuantIndex95,col="red",lwd=2)
lines(TgridIndex,QuantIndex97,col="yellow",lwd=2)
lines(TgridIndex,QuantIndex99,col="brown",lwd=2)
lines(TgridIndex,QuantIndex999,col="green",lwd=2)
lines(TgridIndex,QuantIndex9999,col="blue",lwd=2)
legend("topright",col=c("black","red", "yellow", "brown", "green", "blue"),
       lty=1:3,lwd=2,legend=c("GHI", "95%", "97", "99%", "99.9%", "99.99%"))

write.table(QuantIndex9999,"~/EM9999.txt",sep="\t")

########################################################################################################
##### Calculates the pinball loss score for a given quantile.
####95%, 97%, 99%, 99.9% and 99.99%
pinball_loss <- function(tau, y, q) {
        pl_df <- data.frame(tau = tau,
                            y = y,
                            q = q)
        
        pl_df <- pl_df %>%
                mutate(L = ifelse(y>=q,
                                  tau/100 * (y-q),
                                  (1-tau/100) * (q-y)))
        
        return(pl_df)
}

tau= 99.99
y= GHI
q= QuantIndex9999   # Flstm, Fsvr, Fffnn, Fconvex, Fqra
z = pinball_loss(tau, y, q)
z
write.table(z,"~/pinballFconvex.txt",sep="\t")
qloss =z$L
a=ts(qloss)
plot(a)
mean <- round(mean(qloss),4)
mean
qloss <- round(qloss,4)
write.table(qloss,"~/pinballEM95.txt",sep="\t")

##########################################################################################
########### Data with t

attach(Hourlydata)
head(Hourlydata)

y=GHI
min(y)
max(y)
################Extreme conditional quantiles######### evgam package
library(evgam)

win.graph()
index <- 1:length(y) # index is 1,2,...,4083
plot(y1 ,ylim=c(0,1800),lwd=2,lty=1, ylab = "Hourly GHI (W/sqre m)", 
     xlab = "Observation number")

#################### 0.95-quantile
zeta <- 0.05
Hourlydata$cyc <- as.integer(Hourlydata$t)
FC_fmla_ald <- list(GHI ~ s(cyc, bs = "cc", k = 100), ~ s(cyc, bs = "cc"))
FC_ald <- evgam(FC_fmla_ald, Hourlydata, family = "ald", ald.args = list(tau = 1 - zeta))
Hourlydata$threshold95 <- predict(FC_ald)$location

#################### 0.97-quantile
zeta <- 0.03
Hourlydata$cyc <- as.integer(Hourlydata$t)
FC_fmla_ald <- list(GHI ~ s(cyc, bs = "cc", k = 110), ~ s(cyc, bs = "cc"))
FC_ald <- evgam(FC_fmla_ald, Hourlydata, family = "ald", ald.args = list(tau = 1 - zeta))
Hourlydata$threshold97 <- predict(FC_ald)$location

#################### 0.99-quantile
zeta <- 0.01
Hourlydata$cyc <- as.integer(Hourlydata$t)
FC_fmla_ald <- list(GHI ~ s(cyc, bs = "cc", k = 120), ~ s(cyc, bs = "cc"))
FC_ald <- evgam(FC_fmla_ald, Hourlydata, family = "ald", ald.args = list(tau = 1 - zeta))
Hourlydata$threshold99 <- predict(FC_ald)$location

#################### 0.999-quantile
zeta <- 0.001
Hourlydata$cyc <- as.integer(Hourlydata$t)
FC_fmla_ald <- list(GHI ~ s(cyc, bs = "cc", k = 135), ~ s(cyc, bs = "cc"))
FC_ald <- evgam(FC_fmla_ald, Hourlydata, family = "ald", ald.args = list(tau = 1 - zeta))
Hourlydata$threshold999 <- predict(FC_ald)$location

#################### 0.9999-quantile
zeta <- 0.0001
Hourlydata$cyc <- as.integer(Hourlydata$t)
FC_fmla_ald <- list(GHI ~ s(cyc, bs = "cc", k = 90), ~ s(cyc, bs = "cc"))
FC_ald <- evgam(FC_fmla_ald, Hourlydata, family = "ald", ald.args = list(tau = 1 - zeta))
Hourlydata$threshold9999 <- predict(FC_ald)$location

######## Rounding
GAEV95 <- round(Hourlydata$threshold95,0)
GAEV97 <- round(Hourlydata$threshold97,0)
GAEV99 <- round(Hourlydata$threshold99,0)
GAEV999 <- round(Hourlydata$threshold999,0)
GAEV9999 <- round(Hourlydata$threshold9999,0)

######## Plotting lines
lines(index,GAEV95,col="red",lwd=2)
lines(index,GAEV97,col="yellow",lwd=2)
lines(index,GAEV99,col="brown",lwd=2)
lines(index,GAEV999,col="green",lwd=2)
lines(index,GAEV9999,col="blue",lwd=2)
legend("topright",col=c("black","red", "yellow", "brown", "green", "blue"),
       lty=1:3,lwd=2,legend=c("GHI", "95%", "97", "99%", "99.9%", "99.99%"))

write.table(GAEV97,"~/GAEV97.txt",sep="\t")
########################################################################################################
##### Calculates the pinball loss score for a given quantile.
####95%, 97%, 99%, 99.9% and 99.99%
pinball_loss <- function(tau, y, q) {
        pl_df <- data.frame(tau = tau,
                            y = y,
                            q = q)
        
        pl_df <- pl_df %>%
                mutate(L = ifelse(y>=q,
                                  tau/100 * (y-q),
                                  (1-tau/100) * (q-y)))
        
        return(pl_df)
}

tau= 97
y= GHI
q= GAEV97   # Flstm, Fsvr, Fffnn, Fconvex, Fqra
z = pinball_loss(tau, y, q)
z
write.table(z,"~/pinballFconvex.txt",sep="\t")
qloss =z$L
a=ts(qloss)
plot(a)
mean <- round(mean(qloss),4)
mean
qloss <- round(qloss,4)
write.table(qloss,"~/pinballEM95.txt",sep="\t")

########################################################################################################
###########ADDITIVE QUANTILE REGRESSION 
# Calibrate learning rate on a grid
###### Model 1
win.graph()
index <- 1:length(y) # index is 1,2,...,4083
plot(y1 ,ylim=c(0,1800),lwd=2,lty=1, ylab = "Hourly GHI (W/sqre m)", xlab = "Observation number")

#################### 0.95-quantile 
set.seed(5235)
tun <- tuneLearnFast(form=GHI~s(t, bs="ad"), err = 0.05, qu = 0.95, data = Hourlydata)
tun
fit95 <-qgam(GHI~s(t, bs="ad"), err = 0.05, qu = 0.95, lsig = tun$lsig, data = Hourlydata)
summary(fit95, se="boot") #se =" ker" " nid" "boot"

#################### 0.97-quantile 
set.seed(5235)
tun <- tuneLearnFast(form=GHI~s(t, bs="ad"), err = 0.05, qu = 0.97, data = Hourlydata)
tun
fit97 <-qgam(GHI~s(t, bs="ad"), err = 0.05, qu = 0.97, lsig = tun$lsig, data = Hourlydata)
summary(fit97, se="boot") #se =" ker" " nid" "boot"

#################### 0.99-quantile 
set.seed(5235)
tun <- tuneLearnFast(form=GHI~s(t, bs="ad"), err = 0.05, qu = 0.99, data = Hourlydata)
tun
fit99 <-qgam(GHI~s(t, bs="ad"), err = 0.05, qu = 0.99, lsig = tun$lsig, data = Hourlydata)
summary(fit99, se="boot") #se =" ker" " nid" "boot"

#################### 0.999-quantile 
set.seed(5235)
tun <- tuneLearnFast(form=GHI~s(t, bs="ad"), err = 0.05, qu = 0.999, data = Hourlydata)
tun
fit999 <-qgam(GHI~s(t, bs="ad"), err = 0.05, qu = 0.999, lsig = tun$lsig, data = Hourlydata)
summary(fit999, se="boot") #se =" ker" " nid" "boot"

#################### 0.9999-quantile 
set.seed(5235)
tun <- tuneLearnFast(form=GHI~s(t, bs="ad"), err = 0.05, qu = 0.9999, data = Hourlydata)
tun
fit9999 <-qgam(GHI~s(t, bs="ad"), err = 0.05, qu = 0.9999, lsig = tun$lsig, data = Hourlydata)
summary(fit9999, se="boot") #se =" ker" " nid" "boot"

######## Rounding
AQR95 <- round(fit95$fit,0)
AQR97 <- round(fit97$fit,0)
AQR99 <- round(fit99$fit,0)
AQR999 <- round(fit999$fit,0)
AQR9999 <- round(fit9999$fit,0)

######## Plotting lines
lines(AQR95,col="red",lwd=2)
lines(AQR97,col="yellow",lwd=2)
lines(AQR99,col="brown",lwd=2)
lines(AQR999,col="green",lwd=2)
lines(AQR9999,col="blue",lwd=2)
legend("topright",col=c("black","red", "yellow", "brown", "green", "blue"),
       lty=1:3,lwd=2,legend=c("GHI", "95%", "97", "99%", "99.9%", "99.99%"))

write.table(AQR9999,"~/AQR9999.txt",sep="\t")
########################################################################################################
##### Calculates the pinball loss score for a given quantile.
####95%, 97%, 99%, 99.9% and 99.99%
pinball_loss <- function(tau, y, q) {
        pl_df <- data.frame(tau = tau,
                            y = y,
                            q = q)
        
        pl_df <- pl_df %>%
                mutate(L = ifelse(y>=q,
                                  tau/100 * (y-q),
                                  (1-tau/100) * (q-y)))
        
        return(pl_df)
}

tau= 99.99
y= GHI
q= AQR9999   # Flstm, Fsvr, Fffnn, Fconvex, Fqra
z = pinball_loss(tau, y, q)
z
write.table(z,"~/pinballFconvex.txt",sep="\t")
qloss =z$L
a=ts(qloss)
plot(a)
mean <- round(mean(qloss),4)
mean
qloss <- round(qloss,4)
write.table(qloss,"~/pinballEM95.txt",sep="\t")

################################
###### Model 1
win.graph()
index <- 1:length(y) # index is 1,2,...,4083
plot(y1 ,ylim=c(0,2400),lwd=2,lty=1, ylab = "Hourly GHI (W/sqre m)", xlab = "Observation number")

#################### 0.95-quantile 
set.seed(5235)
tun <- tuneLearnFast(form=GHI~s(t, bs="ad")+s(Temp, bs="ad"), err = 0.05, qu = 0.95, data = Hourlydata)
tun
fit95s <-qgam(GHI~s(t, bs="ad")+s(Temp, bs="ad"), err = 0.05, qu = 0.95, lsig = tun$lsig, data = Hourlydata)
summary(fit95s, se="boot") #se =" ker" " nid" "boot"

#################### 0.97-quantile 
set.seed(5235)
tun <- tuneLearnFast(form=GHI~s(t, bs="ad")+s(Temp, bs="ad"), err = 0.05, qu = 0.97, data = Hourlydata)
tun
fit97s <-qgam(GHI~s(t, bs="ad")+s(Temp, bs="ad"), err = 0.05, qu = 0.97, lsig = tun$lsig, data = Hourlydata)
summary(fit97s, se="boot") #se =" ker" " nid" "boot"

#################### 0.99-quantile 
set.seed(5235)
tun <- tuneLearnFast(form=GHI~s(t, bs="ad")+s(Temp, bs="ad"), err = 0.05, qu = 0.99, data = Hourlydata)
tun
fit99s <-qgam(GHI~s(t, bs="ad")+s(Temp, bs="ad"), err = 0.05, qu = 0.99, lsig = tun$lsig, data = Hourlydata)
summary(fit99s, se="boot") #se =" ker" " nid" "boot"

#################### 0.999-quantile 
set.seed(5235)
tun <- tuneLearnFast(form=GHI~s(t, bs="ad")+s(Temp, bs="ad"), err = 0.05, qu = 0.999, data = Hourlydata)
tun
fit999s <-qgam(GHI~s(t, bs="ad")+s(Temp, bs="ad"), err = 0.05, qu = 0.999, lsig = tun$lsig, data = Hourlydata)
summary(fit999s, se="boot") #se =" ker" " nid" "boot"

#################### 0.9999-quantile 
set.seed(5235)
tun <- tuneLearnFast(form=GHI~s(t, bs="ad")+s(Temp, bs="ad"), err = 0.05, qu = 0.9999, data = Hourlydata)
tun
fit9999s <-qgam(GHI~s(t, bs="ad")+s(Temp, bs="ad"), err = 0.05, qu = 0.9999, lsig = tun$lsig, data = Hourlydata)
summary(fit9999s, se="boot") #se =" ker" " nid" "boot"

######## Rounding
AQR95s <- round(fit95s$fit,0)
AQR97s <- round(fit97s$fit,0)
AQR99s <- round(fit99s$fit,0)
AQR999s <- round(fit999s$fit,0)
AQR9999s <- round(fit9999s$fit,0)

######## Plotting lines
lines(AQR95s,col="red",lwd=2)
lines(AQR97s,col="yellow",lwd=2)
lines(AQR99s,col="brown",lwd=2)
lines(AQR999s,col="green",lwd=2)
lines(AQR9999s,col="blue",lwd=2)
legend("topright",col=c("black","red", "yellow", "brown", "green", "blue"),
       lty=1:3,lwd=2,legend=c("GHI", "95%", "97", "99%", "99.9%", "99.99%"))

write.table(AQR9999s,"~/AQR9999s.txt",sep="\t")
########################################################################################################
##### Calculates the pinball loss score for a given quantile.
####95%, 97%, 99%, 99.9% and 99.99%
pinball_loss <- function(tau, y, q) {
        pl_df <- data.frame(tau = tau,
                            y = y,
                            q = q)
        
        pl_df <- pl_df %>%
                mutate(L = ifelse(y>=q,
                                  tau/100 * (y-q),
                                  (1-tau/100) * (q-y)))
        
        return(pl_df)
}

tau= 99.9
y= GHI
q= AQR999s   # Flstm, Fsvr, Fffnn, Fconvex, Fqra
z = pinball_loss(tau, y, q)
z
write.table(z,"~/pinballFconvex.txt",sep="\t")
qloss =z$L
a=ts(qloss)
plot(a)
mean <- round(mean(qloss),4)
mean
qloss <- round(qloss,4)
write.table(qloss,"~/pinballEM95.txt",sep="\t")

#######################################################################################################
###########Average model 
attach(Fdata)
head(Fdata)

y=GHI
min(y)
max(y)
#####################

win.graph()
index <- 1:length(y) # index is 1,2,...,4083
y1 <- ts(GHI)
plot(y1 ,ylim=c(0,2300),lwd=2,lty=1, ylab = "Hourly GHI (W/sqre m)", xlab = "Observation number")

######## Plotting lines
lines(AV95,col="red",lwd=2)
lines(AV97,col="yellow",lwd=2)
lines(AV99,col="brown",lwd=2)
lines(AV999,col="green",lwd=2)
lines(AV9999,col="blue",lwd=2)
legend("topleft",col=c("black","red", "yellow", "brown", "green", "blue"),
       lty=1:3,lwd=2,legend=c("GHI", "95%", "97", "99%", "99.9%", "99.99%"))

########################################################################################################
##### Calculates the pinball loss score for a given quantile.
####95%, 97%, 99%, 99.9% and 99.99%
pinball_loss <- function(tau, y, q) {
        pl_df <- data.frame(tau = tau,
                            y = y,
                            q = q)
        
        pl_df <- pl_df %>%
                mutate(L = ifelse(y>=q,
                                  tau/100 * (y-q),
                                  (1-tau/100) * (q-y)))
        
        return(pl_df)
}

tau= 99.99
y= GHI
q= AV9999   # Flstm, Fsvr, Fffnn, Fconvex, Fqra
z = pinball_loss(tau, y, q)
z
write.table(z,"~/pinballFconvex.txt",sep="\t")
qloss =z$L
a=ts(qloss)
plot(a)
mean <- round(mean(qloss),4)
mean
qloss <- round(qloss,4)
write.table(qloss,"~/pinballEM95.txt",sep="\t")

#####################Median model

win.graph()
index <- 1:length(y) # index is 1,2,...,4083
y1 <- ts(GHI)
plot(y1 ,ylim=c(0,1900),lwd=2,lty=1, ylab = "Hourly GHI (W/sqre m)", xlab = "Observation number")

######## Plotting lines
lines(MED95,col="red",lwd=2)
lines(MED97,col="yellow",lwd=2)
lines(MED99,col="brown",lwd=2)
lines(MED999,col="green",lwd=2)
lines(MED9999,col="blue",lwd=2)
legend("topleft",col=c("black","red", "yellow", "brown", "green", "blue"),
       lty=1:3,lwd=2,legend=c("GHI", "95%", "97", "99%", "99.9%", "99.99%"))

########################################################################################################
##### Calculates the pinball loss score for a given quantile.
####95%, 97%, 99%, 99.9% and 99.99%
pinball_loss <- function(tau, y, q) {
        pl_df <- data.frame(tau = tau,
                            y = y,
                            q = q)
        
        pl_df <- pl_df %>%
                mutate(L = ifelse(y>=q,
                                  tau/100 * (y-q),
                                  (1-tau/100) * (q-y)))
        
        return(pl_df)
}

tau= 99.99
y= GHI
q= MED9999   # Flstm, Fsvr, Fffnn, Fconvex, Fqra
z = pinball_loss(tau, y, q)
z
write.table(z,"~/pinballFconvex.txt",sep="\t")
qloss =z$L
a=ts(qloss)
plot(a)
mean <- round(mean(qloss),4)
mean
qloss <- round(qloss,4)
write.table(qloss,"~/pinballEM95.txt",sep="\t")

########################################################################################################
##### Calculates the crps for a given quantile.
####95%, 97%, 99%, 99.9% and 99.99%
a<- crps_gev(MED9999, shape=0.6723277 , location = -27571.71, scale =297.2173438)
mean <- round(mean(a),4)
mean

#######################################################################################################
################# Best model plotting

win.graph()
index <- 1:length(y) # index is 1,2,...,4083
y1 <- ts(GHI)

######## Plotting lines
plot(y1 ,ylim=c(0,1400),lwd=2,lty=1, ylab = "Hourly GHI (W/sqre m)", xlab = "Observation number")
lines(AQR95s,col="red",lwd=1)
legend("topleft",col=c("black","red"), lty=1:3,lwd=2,legend=c("GHI", "95%"))

plot(y1 ,ylim=c(0,1400),lwd=2,lty=1, ylab = "Hourly GHI (W/sqre m)", xlab = "Observation number")
lines(AQR97s,col="red",lwd=1)
legend("topleft",col=c("black","red"), lty=1:3,lwd=2,legend=c("GHI", "97%"))

plot(y1 ,ylim=c(0,1400),lwd=2,lty=1, ylab = "Hourly GHI (W/sqre m)", xlab = "Observation number")
lines(AQR99s,col="red",lwd=1)
legend("topleft",col=c("black","red"), lty=1:3,lwd=2,legend=c("GHI", "99%"))

plot(y1 ,ylim=c(0,1650),lwd=2,lty=1, ylab = "Hourly GHI (W/sqre m)", xlab = "Observation number")
lines(AQR999s,col="red",lwd=1)
legend("topleft",col=c("black","red"), lty=1:3,lwd=2,legend=c("GHI", "99.9%"))

plot(y1 ,ylim=c(0,1500),lwd=2,lty=1, ylab = "Hourly GHI (W/sqre m)", xlab = "Observation number")
lines(GA9999,col="red",lwd=1)
legend("topright",col=c("black","red"), lty=1:3,lwd=2,legend=c("GHI", "99.99%"))


#######################################################################################################
attach(pindata)
head(pindata)



#####ts and density plots for pinball losses for best model
win.graph()
par(mfrow=c(1,2))
a1=ts(EV999)
plot(a1, main="(a) Plot of Pinball losses",xlab = "Observation number",ylab="Pinball losses (W/sqre m)"
     ,col="blue")
vGHI=density(EV999)
plot(vGHI, main="(b) Density of Pinball losses",xlab = "Pinball losses (W/sqre m)",col="blue")

######Bloxplot
win.graph()
modl <- c("Extremefit", "Evgam", "Benchmark")
boxplot(EX999, EV999, NL999, names= modl, horizontal = FALSE, 
        main="",ylab="Pinball losses (W/sqre m)", col = "blue")
