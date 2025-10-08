#加载包
library(zoo)
library(ggplot2)
library(dlnm)
library(mgcv)
library(splines)
library(plot3D)
library(gridExtra)


file_path_1 = "C:/Users/Huang/Desktop/iowa2013-2024_1.csv" 
data1 <- read.csv(file_path_1)
beta_iowa = data1[c(15:322),26]
climate_IA = read.csv("C:/Users/Huang/Desktop/iowa2013-2020daily.csv")
temp_iowa = as.numeric(climate_IA[c(276:2425),10])
hum_iowa = as.numeric(climate_IA[c(276:2425),8])
time_iowa = as.numeric(climate_IA[c(276:2425),1])
weeks = data1[c(15:322),27]
start_date <- as.Date("2013-9-30")
week_dates <- start_date + (weeks - 1) * 7
daily_dates <- seq(from = start_date, to = max(week_dates), by = "day")


daily_values1 <- approx(week_dates, beta_iowa, xout = daily_dates)$y


daily_data <- data.frame(Date = daily_dates, Value = daily_values1)

beta_iowa = as.numeric(daily_values1)
cb_temp_ia = crossbasis(temp_iowa, lag=c(0,7), argvar=list(fun="ns",df=3), arglag=list(fun="bs",df=4))
cb_hum_ia = crossbasis(hum_iowa, lag=c(0,7), argvar=list(fun="bs",df=6), arglag=list(fun="ns",df=4))

model_IA = glm(beta_iowa ~ cb_temp_ia  +cb_hum_ia + ns(daily_dates,6*2),  family = quasipoisson(),climate_IA)


summary(model_IA)

pred.temp_ia = crosspred(cb_temp_ia, model_IA, cen=round(median(temp_iowa)), bylag=0.2)


crall11 <- crossreduce(cb_temp_ia,model_IA,cen=round(median(temp_iowa)),type="overall")


file_path_3 <- "C:/Users/Huang/Desktop/Illinois2013-2024_1.csv"  
data3 <- read.csv(file_path_3)
beta_il = data3[c(40:347),6]
climate_IL = read.csv("C:/Users/Huang/Desktop/illinois2013-2020daily.csv")
temp_il = as.numeric(climate_IL[c(276:2425),10])
hum_il = as.numeric(climate_IL[c(276:2425),8])
time_il = as.numeric(climate_IL[c(276:2425),1])
start_date <- as.Date("2013-9-30")
week_dates <- start_date + (weeks - 1) * 7



daily_dates <- seq(from = start_date, to = max(week_dates), by = "day")


daily_values3 <- approx(week_dates, beta_il, xout = daily_dates)$y


daily_data <- data.frame(Date = daily_dates, Value = daily_values3)
beta_il = as.numeric(daily_values3)

cb_temp_il = crossbasis(temp_il, lag=c(0,7), argvar=list(fun="ns",df=3), arglag=list(fun="bs",df=4))
cb_hum_il = crossbasis(hum_il, lag=c(0,7), argvar=list(fun="bs",df=6), arglag=list(fun="ns",df=4))
model_IL = glm(beta_il ~ cb_temp_il + cb_hum_il + ns(daily_dates,6*2),  family = quasipoisson(), climate_IL)
summary(model_IL)
pred.temp_il = crosspred(cb_temp_il, model_IL, cen=round(median(temp_il)), bylag=0.2)




crall31 <- crossreduce(cb_temp_il,model_IL,cen=round(median(temp_il)),type="overall")


file_path_4 <- "C:/Users/Huang/Desktop/IN2013-2019_1.csv" 
data4 <- read.csv(file_path_4)
beta_in = data4[c(3:310),6]
climate_IN = read.csv("C:/Users/Huang/Desktop/IN.csv")
temp_in = as.numeric(climate_IN[c(276:2425),10])
hum_in = as.numeric(climate_IN[c(276:2425),8])
time_in = as.numeric(climate_IN[c(276:2425),1])
start_date <- as.Date("2013-9-30")
week_dates <- start_date + (weeks - 1) * 7



daily_dates <- seq(from = start_date, to = max(week_dates), by = "day")


daily_values4 <- approx(week_dates, beta_in, xout = daily_dates)$y


daily_data <- data.frame(Date = daily_dates, Value = daily_values4)
beta_in = as.numeric(daily_values4)

cb_temp_in = crossbasis(temp_in, lag=c(0,7), argvar=list(fun="ns",df=3), arglag=list(fun="bs",df=4))
cb_hum_in = crossbasis(hum_in, lag=c(0,7), argvar=list(fun="bs",df=6), arglag=list(fun="ns",df=4))
model_IN = glm(beta_in ~ cb_temp_in + cb_hum_in + ns(daily_dates,6*2),  family = quasipoisson(), climate_IN)


summary(model_IN)

pred.temp_in = crosspred(cb_temp_in, model_IN, cen=round(median(temp_in)), bylag=0.2)




crall41 <- crossreduce(cb_temp_in,model_IN,cen=round(median(temp_in)),type="overall")








png("C:/Users/Huang/Desktop/IA3d.png", width = 3000, height = 2500, res=350) 
#par(mfrow=c(1,3), mar = c(6, 7, 5, 5))par(mgp = c(3, 2, 0))
par(mar = c(2, 2, 2, 2),mgp = c(5, 2, 0))

plot(pred.temp_ia,ticktype='detailed',border='#3366FF',xlab="(°C)",ylab="day",zlab="RR",col='#99FFCC',shade = 0.1 ,cex.lab=1.2,cex.axis=1.2,lwd=1,theta = 20, phi = 25,ltheta = -35)
dev.off()
png("C:/Users/Huang/Desktop/IL3d.png", width = 3000, height = 2500, res=350)
par(mar = c(2, 2, 2, 2))
plot(pred.temp_il,ticktype='detailed',border='#3366FF',xlab="(°C)",ylab="day",zlab="RR",col='#99FFCC',shade = 0.1,cex.lab=1.2,cex.axis=1.2,lwd=1,theta = 20, phi = 25,ltheta = -35)

dev.off()
png("C:/Users/Huang/Desktop/IN3d.png", width = 3000, height = 2500, res=350)
par(mar = c(2, 2, 2, 2))
plot(pred.temp_in,ticktype='detailed',border='#3366FF',xlab="(°C)",ylab="day",zlab="RR",col='#99FFCC',shade = 0.1,cex.lab=1.2,cex.axis=1.2,lwd=1,theta = 20, phi = 25,ltheta = -35)

dev.off()





png("C:/Users/Amax/Desktop/IA.png", width = 1600, height = 1600, res=300) 
#par(mfrow=c(1,3), mar = c(6, 7, 5, 5))
plot(crall11,xlab="(°C)",ylab="RR",col=2,lwd=2,cex.lab=1.2,cex.axis=1.2,mar=c(1,2,0,1))

dev.off()
png("C:/Users/Amax/Desktop/IL.png", width = 1600, height = 1600, res=300) 
plot(crall31,xlab=" (°C)",ylab="RR",col=2,lwd=2,cex.lab=1.2,cex.axis=1.2,mar=c(1,2,0,1))

dev.off()
png("C:/Users/Huang/Desktop/IN.png", width = 1600, height = 1600, res=300) 
plot(crall41,xlab="(°C)",ylab="RR",col=2,lwd=2,cex.lab=1.2,cex.axis=1.2,mar=c(1,2,0,1))
dev.off()

yall_temp = matrix(NA,3,3)
yall_temp[1,] = coef(crall11)
yall_temp[2,] = coef(crall41)
yall_temp[3,] = coef(crall31)
sall_temp = vector("list",3)
sall_temp[[1]] = vcov(crall11)
sall_temp[[2]] = vcov(crall41)
sall_temp[[3]] = vcov(crall31)
library(mvmeta)

method <- "reml"

mvall_temp <- mvmeta(yall_temp~1,sall_temp,method=method)
mvall_tempgdp = mvmeta(yall_temp~gdp,sall_temp,method=method)
mvall_temppop = mvmeta(yall_temp~pop,sall_temp,method=method)

summary(mvall_temp)
argvar1 = attr(cb_temp_ia, "argvar")

bvar_temp=do.call("onebasis", c(list(x=seq(-22,27,0.5)), argvar1))
cpall_temp = crosspred(bvar_temp, coef=coef(mvall_temp), vcov=vcov(mvall_temp), 
                       model.link="log", by=0.5,bylag=0.2, cen=13)


png("meta.png", width=2000, height=1500, res=300)
plot(cpall_temp, xlab="(°C)", ylab="RR", lwd=2, col='red')
for(i in 1:3){
  cp = crosspred(bvar_temp, coef=yall_temp[i,], vcov=sall_temp[[i]] ,cen=13,model.link="log")
  lines(cp, lty=i, lwd=2, col="grey")
}

dev.off()







