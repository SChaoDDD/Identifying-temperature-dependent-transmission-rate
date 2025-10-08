rm(list = ls())
gc()

library("readxl")
library(readr)
library(rEDM)
library(ggplot2)
setwd("C:/Users/XLYC/Desktop/Rprom/liuganccm")
data01<-read_excel("ILcase.xlsx")
data02<-read_excel('ILtemp.xlsx')
x<-as.matrix(data01[2:308,1])
y<-as.matrix(data02[1:307,1])


data<-data.frame(time = 1:length(x), x = x, y = y)
lags <- -5:5
time_lag <- matrix(NA,length(lags),3)
colnames(time_lag) <- c("tp", "Temp:Cases", "Cases:Temp")
time_lag[,1] <- lags
for (i in 1:length(lags)){
  ccm_result <- CCM(dataFrame = data,
                    E = 4,
                    Tp = lags[i],
                    tau = -1, 
                    columns = "temp",
                    target = "ILcase",
                    random = TRUE,
                    sample = 10,
                    seed = 0,
                    libSizes = "290",
                    showPlot = FALSE)
  time_lag[i,2:3] <- as.vector(as.matrix(ccm_result[1,2:3]))
}
time_lag <- as.data.frame(time_lag)
ggplot(time_lag, aes(x = tp)) +
  geom_line(aes(y = `Cases:Temp`, color = "Cases:Temp"), linewidth = 1.5) +
  geom_line(aes(y = `Temp:Cases`, color = "Temp:Cases"), linewidth = 1.5) +
  geom_point(aes(y = `Cases:Temp`, color = "Cases:Temp"), size = 3, shape = 21, fill = "white", stroke = 1.5) +  # 添加数据点
  geom_point(aes(y = `Temp:Cases`, color = "Temp:Cases"), size = 3, shape = 21, fill = "white", stroke = 1.5) +
  theme_bw() +
  labs(y = "Cross map skill",x = "Time lag", color = "Iowa") +
  scale_color_manual(values = c("Cases:Temp" = "red", "Temp:Cases" = "blue")) +
  theme(legend.position= c(0.85, 0.6),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 24), 
        axis.text.y = element_text(size = 24),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.title.align = 0.5
  )



# ccm_result <- CCM(dataFrame = data,
#                   E = 4,
#                   Tp = -1,
#                   tau = -1, 
#                   columns = "IARH",
#                   target = "beta",
#                   random = TRUE,
#                   sample = 10,
#                   seed = 0,
#                   libSizes = "6,295,1",
#                   showPlot = FALSE)
# fanxiang<-as.array(ccm_result$`IARH:beta`)
# zhengxiang<-as.array(ccm_result$`beta:IARH`)
# librarysize<-as.array(ccm_result$LibSize)
# data <- data.frame(x = librarysize, y1 = zhengxiang, y2 = fanxiang)
# 
# #绘制图形
# ggplot(data, aes(x = x)) +
#   geom_line(aes(y = y1, color = "Cases:RH"), linewidth = 1) +
#   geom_line(aes(y = y2, color = "RH:Cases"), linewidth = 1) +
#   #geom_point(aes(y = y1, color = "Cases:Temp"), size = 3, shape = 21, fill = "white", stroke = 1.5)+
#   labs(x = "Length of library", y = "Cross map skill", color = "Iowa") +
#   scale_color_manual(values = c("Cases:RH" = "red", "RH:Cases" = "blue")) +
#   theme_minimal()+
#   theme_bw() +
#   theme(legend.position= c(0.85, 0.55),
#         axis.title.x = element_text(size = 24),
#         axis.title.y = element_text(size = 24),
#         axis.text.x = element_text(size = 24), 
#         axis.text.y = element_text(size = 24),
#         legend.text = element_text(size = 18),
#         legend.title = element_text(size = 16),
#         legend.title.align = 0.5
#   )

