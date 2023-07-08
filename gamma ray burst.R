#Data collection

library(cluster)
library(factoextra)
set.seed(100)
astro_data = read.csv(file.choose())
astro_data = astro_data[,-c(14:17)]
#a = head(astro_data)
#is.numeric(astro_data)
#is.numeric(astro_data$TrigNo)
#is.numeric(astro_data$F2)
#a
#tail(astro_data)
attach(astro_data)

#data imputation and missing value treatment

library(mice)
library(missForest)
astr.mis = prodNA(astro_data, noNA = 0.1)
#summary(astr.mis)
#md.pattern(astr.mis)
library(VIM)
mice_plot = aggr(astro_data, col=c('navyblue','yellow'),
                 numbers=TRUE, sortVars=TRUE,
                 labels=names(astro_data), cex.axis=.7,
                 gap=3, ylab=c("Missing data","Pattern"))
imputed_Data = mice(astr.mis, m=3, maxit = 1000, method = 'pmm', seed = 500)
complete_data = complete(imputed_Data , 3)
head(complete_data)
astro_data_log = log(complete_data[,c(2:13)])
head(astro_data_log)
colnames(astro_data_log) = c("logF1","logF2","logF3","logF4","logF64","logF256",
                             "logF1024","logFT","logH32","logH321","logT50","logT90")
head(astro_data_log)
comb_astro_data = data.frame(complete_data , astro_data_log)
head(comb_astro_data)
final_astro_data = comb_astro_data[,-c(2:8 , 10:13)]
head(final_astro_data)
attach(final_astro_data)


#Distortion curve
#
set.seed(100)
wss = (nrow(final_astro_data)-1)*sum(apply(final_astro_data , 2, var))
for (i in 2:13)
{
  wss[i] =
    sum(kmeans(final_astro_data , centers = i)$withinss)
}
plot(1:13 , wss , type = "l" , 
     xlab = "number of clusters" ,
     ylab = "within group of sum of  square")
#
#
#
fviz_nbclust(final_astro_data, kmeans, method = "wss")
#
#
#
#
set.seed(100)
km_res = kmeans(final_astro_data, centers = 2, nstart = 20)

sil = silhouette(km_res$cluster, dist(astro_data))
fviz_silhouette(sil)
#
#
#
library(fpc)
plotcluster(final_astro_data , km_res$cluster)

fviz_cluster(km_res,data=final_astro_data)
#
#
clus = cbind(final_astro_data , km_res$cluster)
clus
summary(clus)
summary(km_res)
#
#
#
astro_cluster = data.frame(final_astro_data, cluster = as.factor(km_res$cluster))
head(astro_cluster)
#
# LDA
set.seed(100)
library(MASS)
group = c(rep (1, 907), rep (2, 1055))
head(group) 
discr = lda (final_astro_data[,c(8,10:14)],group)
tab1 = table(predict(discr)$class, group,dnn = c('Actual Group','Predicted Group'))
tab1

nrow(final_astro_data)
# EDA and Distribution Analysis

set.seed(100)

x1=logT50
x2=logT90

hist(x1)
hist(x2)

shapiro.test(x1)
shapiro.test(x2)

ks.test(x1,x2)

cor(x1,x2)
plot(x1,x2,pch=19)


x3=logF64
x4=logF256
x5=logF1024

d345=data.frame(x3,x4,x5)

hist(x3)
hist(x4)
hist(x5)

shapiro.test(x3)
shapiro.test(x4)
shapiro.test(x5)

ks.test(x3,x4)
ks.test(x3,x5)
ks.test(x4,x5)

cor(d345)

x6 = logH32 
x7 = logH321
x8 = logF1
x9 = logF2
x10 = logF3
x11 = logF4
cor(x6 , x7)
shapiro.test(x6)
shapiro.test(x7)
shapiro.test(x8)
shapiro.test(x9)
shapiro.test(x10)
shapiro.test(x11)
hist(x6)
hist(x7)
hist(x8)
hist(x9)
hist(x10)
hist(x11)
ks.test(x6 , x7)
ks.test(x8 , x9)
ks.test(x8 , x10)
ks.test(x8 , x11)
ks.test(x9 , x10)
ks.test(x9 , x11)
ks.test(x10 , x11)
df891011 = data.frame(logF1 , logF2 , logF3 ,logF4)
cor(df891011)

x12=logFT

hist(x12)
shapiro.test(x12)

cor(x1,x12)
cor(x2,x12)

library('fpc')
library('MVN')

result=MVN::mvn(data=final_astro_data,mvnTest="hz",univariateTest="AD",univariatePlot="histogram",multivariatePlot="qq",multivariateOutlierMethod="adj",showOutliers=TRUE,showNewData=TRUE)
result
summary(results)
rcorr(as.matrix(final_astro_data[c(9:13 , 6)]))