#Data collection
library(readxl)
library(readr)
library(cluster)
library(factoextra)
set.seed(100)
astro_data = read_csv("astro_data.csv")
astro_data = data.frame(astro_data)
astro_data = astro_data[,-c(14:17)]
head(astro_data)
#a = head(astro_data)
#is.numeric(astro_data)
#is.numeric(astro_data$TrigNo)
#is.numeric(astro_data$F2)
#a
#tail(astro_data)
attach(astro_data)
sum(is.na(astro_data))
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
library(mice)
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
#fviz_nbclust(final_astro_data, kmeans, method = "wss")
#
library(clv)
Dunn <- function(data,clust) 
  clv.Dunn( cls.scatt.data(data,clust),
            intracls = c("complete","average","centroid"), 
            intercls = c("single", "complete", "average","centroid", "aveToCent", "hausdorff")
  )
Davies.Bouldin <- function(data,clust) 
  clv.Davies.Bouldin( cls.scatt.data(data,clust),
                      intracls = c("complete","average","centroid"),
                      intercls = c("single", "complete", "average","centroid", "aveToCent", "hausdorff")
  )
#
#
fviz_nbclust(final_astro_data, kmeans, method = "wss")
#
#
set.seed(100)
km_res = kmeans(final_astro_data, centers = 3, nstart = 20)

sil = silhouette(km_res$cluster, dist(astro_data))
fviz_silhouette(sil)
Dunn(final_astro_data , km_res$cluster)
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
km_res$size
#
#
#
astro_cluster = data.frame(final_astro_data, cluster = as.factor(km_res$cluster))
head(astro_cluster)
#
# LDA
set.seed(100)
library(MASS)
group = c(rep (1, 711), rep (2, 901), rep(3 , 350))
head(group) 
discr = lda(final_astro_data[,c(8,10:14)],group)
tab1 = table(predict(discr)$class, group,dnn = c('Actual Group','Predicted Group'))
tab1

nrow(final_astro_data)

Data = read_csv("GRB.csv")
Data = data.frame(Data)
head(Data)
dim(Data)
dim(final_astro_data)
merged_dataset <- merge(Data, final_astro_data, by = "TrigNo", all = TRUE)
head(merged_dataset) 
final_data = na.omit(merged_dataset[,-4])
head(final_data)
dim(final_data)
library(circular)
library(CircStats)
watson.test(final_data$GLON , alpha = 0.05 , dist = "uniform")
watson.test(final_data$GLAT , alpha = 0.05 , dist = "uniform")
watson.test(final_data$RAJ2000 , alpha = 0.05 , dist = "vonmises")
watson.test(final_data$DEJ2000 , alpha = 0.05 , dist = "vonmises")
watson.test(final_data$Angle , alpha = 0.05 , dist = "vonmises")
head(final_data)
library(ggplot2)
library(sp)
# Assuming 'df' contains the equatorial coordinates
L_Cluster_Data = final_data[c(13:25)]
head(L_Cluster_Data)
Spherical_Clustering = cbind(final_data$RAJ2000 , final_data$DEJ2000 , final_data$GLON,
                             final_data$GLAT , final_data$Angle)
head(Spherical_Clustering)
library(plotly)
# Define the galactic coordinates of your object
object_glon <- final_data$GLON
object_glat <- final_data$GLAT

# Convert the galactic coordinates to Cartesian coordinates
object_x <- cos(object_glat * pi / 180) * cos(object_glon * pi / 180)
object_y <- cos(object_glat * pi / 180) * sin(object_glon * pi / 180)
object_z <- sin(object_glat * pi / 180)

# Create the universe sphere plot
universe_plot <- plot_ly(type = "scatter3d", mode = "markers")

# Add the universe sphere
theta <- seq(0, 2 * pi, length.out = 100)
phi <- seq(0, pi, length.out = 50)
grid_theta <- rep(theta, length(phi))
grid_phi <- rep(phi, each = length(theta))
x <- cos(grid_phi) * cos(grid_theta)
y <- cos(grid_phi) * sin(grid_theta)
z <- sin(grid_phi)
universe_plot <- add_trace(universe_plot, x = x, y = y, z = z, mode = "lines", line = list(color = "white"))

# Add the object
universe_plot <- add_trace(universe_plot, x = object_x, y = object_y, z = object_z, mode = "markers",
                           marker = list(size = 5, color = "red"))

# Customize the plot layout
universe_plot <- layout(universe_plot, scene = list(aspectmode = "data",
                                                    xaxis = list(title = "X"),
                                                    yaxis = list(title = "Y"),
                                                    zaxis = list(title = "Z")),
                        title = "Object in Galactic Coordinate System")

# Display the plot
universe_plot
#
#
grb_plot <- ggplot(final_data) +
  coord_map("mollweide") +  # Set the coordinate system to Mollweide projection
  #geom_path(data = map_data("world"), aes(x = long, y = lat, group = group), color = "gray") +  # Add map outline
  geom_point(aes(x = final_data$GLON, y = final_data$GLAT), color = "blue",fill = "black", size = 1) +  # Plot the GRBs as points
  theme_void()  # Remove unnecessary plot elements
grb_plot
#
#
set.seed(100)
km_res = kmeans(L_Cluster_Data, centers = 2, nstart = 20)
km_res3 = kmeans(L_Cluster_Data, centers = 3, nstart = 20)
km_res4 = kmeans(L_Cluster_Data, centers = 4, nstart = 20)

Dunn(final_astro_data , km_res$cluster)
Dunn(final_astro_data , km_res3$cluster)
Dunn(final_astro_data , km_res4$cluster)

library(skmeans)
SphCl3 = skmeans(Spherical_Clustering , 3)
