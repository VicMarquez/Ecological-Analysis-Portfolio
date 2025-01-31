# Ecological diversity analysis :
# Libraries and Setup
library(vegan)
library(lattice)
library(permute)
library (vegan3d)
library (dplyr)
library (labdsv)
library(MASS)
library (ggplot2)
install.packages('jsonlite', repos='http://cran.rstudio.com/')


###### DIVIDE THE COMPOSITION BY STRATA ######

#########tree strata ##########
tree<- read.table(file.choose(), header = TRUE)
head(tree) # data set: veg_tree

###Analysis of Similarity (ANOSIM)
# for abundance data, bray distance is used, the first cell of the database must be EMPTY. 
tree.dist<-vegdist(tree, "bray")
summary(tree.dist)
#I load the sites info.
sitios<- read.table(file.choose(), header = TRUE)
attach(sitios)

dispersion<-betadisper(tree.dist, group=Condicion) # Condicion= fixed factor (C and S)
permutest(dispersion) 
#If the permutation test (permutest) or Tukey HSD test (mod.HSD) does not show significant 
#differences (e.g., p > 0.05), it suggests that group dispersions are similar. 
#This means you can proceed with further analyses (e.g., ANOSIM, PERMANOVA)
#without worrying about heterogeneity of variances.

mod.HSD <- TukeyHSD(dispersion) # non significant 

#The plot (plot(mod.HSD)) helps visually interpret the pairwise comparisons, 
#making it easier to understand which groups (if any) differ in dispersion.

plot(mod.HSD)

# Performs ANOSIM to test for significant differences in species composition between groups (e.g., Condicion).
data.ano.t <- anosim(tree.dist, Condicion)
summary(data.ano.t)

#Dissimilarity ranks between and within classes:
#  0%  25%  50%   75% 100%  N
#Between  3 22.5 38.5 52.25   66 36
#P        1  7.0 17.0 36.50   55 15
#S        6 21.5 33.0 54.00   64 15

plot(data.ano.t)
#R=0.24 P=0.06

#The ANOSIM statistic compares the mean of ranked dissimilarities between groups to 
#the mean of ranked dissimilarities within groups. An R value close to "1.0" suggests
#dissimilarity between groups while an R value close to "0" suggests an even distribution 
#of high and low ranks within and between groups. R values below "0" suggest that dissimilarities are greater within groups
#than between groups. See Clarke and Gorley (2001) for a guide to interpreting ANOSIM R values.

# test the effect of Condicion on species composition

veg.adonis <- adonis2(tree.dist ~ Condicion,  contr.unordered = "contr.sum", data = tree)
veg.adonis

#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#Condicion  1   0.23940 0.239399  2.8093 0.21932   0.07 .
#Residuals 10   0.85215 0.085215         0.78068         
#Total     11   1.09155                  1.00000
      
##############################################################################
#Non-metric multidimensional scaling NMDS
# Load necessary libraries
library(ggplot2)
ord1 <- metaMDS(tree, distance = "bray")
plot(ord1, type = "t")


data.scores <- as.data.frame(scores(ord1))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(ord1)  # create a column of site names, from the rownames of data.scores
head(data.scores)  #look at the data
species.scores <- as.data.frame(scores(ord1, "species"))  
species.scores$species <- rownames(species.scores)  
head(species.scores)  #look at the data
data.scores <- as.data.frame(scores(ord1))  
data.scores$site <- rownames(data.scores)
data.scores$Condicion <- Condicion
head(data.scores)  


ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, shape=Condicion, colour=Condicion),size=3) + 
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=SitioA),size=6,vjust=0) +  
  scale_colour_manual(values=c("S" = "red", "P" = "blue")) +
  coord_equal() +
  theme_bw()

#make it nicer

ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Condicion,colour=Condicion),size=4) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Condicion),size=8,vjust=0,hjust=0) +  # add the site labels
  scale_colour_manual(values=c("S" = "red", "P" = "blue")) +
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), #
        axis.ticks = element_blank(),  
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())


c.s <- data.scores[data.scores$Condicion == "S", ][chull(data.scores[data.scores$Condicion == 
                                                                    "S", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
c.p <- data.scores[data.scores$Condicion== "P", ][chull(data.scores[data.scores$Condicion == 
                                                                   "P", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
hull.data <- rbind(c.s, c.p)  

hull.data


ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Condicion,group=Condicion),alpha=0.30) + # add the convex hulls
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=SitioA),size=4) + # add the point markers
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())




library(ggplot2)
p1 <-ggplot(data.scores, aes(NMDS1, NMDS2))+
  geom_point(aes(NMDS1, NMDS2, color=Condicion),position=position_jitter(.1))+##separates overlapping points
  stat_ellipse(aes(fill=Condicion), alpha=.2,type='t',size =1, geom="polygon")+ ##changes shading on ellipses
  theme_minimal()

p1



#el stress indice permite saber si las diferencias y similitudes de los datos
#�estan bien expresadas en el graficos. En este caso si. hay que fijrse q tan dispersos 
#estan los puntos de la linea

stressplot(ord1)

#Clarke 1993 suggests the following guidelines for acceptable stress values: 
#<0.05 = excellent, <0.10 = good, <0.20 = usable, >0.20 = not acceptable. The plot shows the border of the 0.20 stress value limit. Solutions with higher 
#stress values should be interpreted with caution and those with stress above 0.30 are highly suspect.

metaMDS(
  data,
  distance = "bray",
  k = 2,
  trymax = 20,
  autotransform = TRUE) #  0.09775674  el valor de stress es bueno

#################################
require(graphics)
install.packages("ggdendro")
library("ggplot2")
library(ggdendro)
head(data)


a <- hclust(data.dist, "ave") #########mepa q es este 
a <- hclust(dist(data), "ave")


a <- hclust(data.dist, "ave")
plot(a)
plot(a, hang = -1) # aca se ve que lo sitios estan re mezclados 
# la condicion no explica ni la variable del suelo ni la composicion de la vegetacion 

ggdendrogram(a)

# Build dendrogram object from hclust results
dend <- as.dendrogram(a)
# Extract the data (for rectangular lines)
# Type can be "rectangle" or "triangle"
dend_data <- dendro_data(dend, type = "rectangle")
# What contains dend_data
names(dend_data)
head(dend_data$segments)
head(dend_data$labels)

###############ESTRATO ARBUSTIVO####################
shrub<- read.table(file.choose(), header = TRUE)
head(shrub) # base de datos veg_shrub


#Analysis of Similarities
shrub.dist<-vegdist(shrub, "bray") # para datos de abundancia se usa bray, la primer celda de la base de datos debe estar VACIA 
shrub.dist
summary(shrub.dist)

sitios<- read.table(file.choose(), header = TRUE)#cargo la info de los sitios. Archivo"Sitios.txt" o"sitios+alt+DAS"
attach(sitios)

dispersion<-betadisper(shrub.dist, group=Condicion)
permutest(dispersion) 
mod.HSD <- TukeyHSD(dispersion) # no da significativa puedo seguir 
plot(mod.HSD)

#La matriz.de.distancias es la que obtienes con vegdist
#El grupo en tu caso sería la condic


data.ano.s <- anosim(shrub.dist, Condicion)

summary(data.ano.s)

#Dissimilarity ranks between and within classes:
#  0%   25%  50%   75% 100%  N
#Between  1 17.75 34.5 53.75   66 36
#P        7 21.00 31.0 36.50   44 15
#S        2 11.00 45.0 53.00   60 15

plot(data.ano.s)
#R=0.12. P=0.158

#The ANOSIM statistic compares the mean of ranked dissimilarities between groups to 
#the mean of ranked dissimilarities within groups. An R value close to "1.0" suggests
#dissimilarity between groups while an R value close to "0" suggests an even distribution 
#of high and low ranks within and between groups. R values below "0" suggest that dissimilarities are greater within groups
#than between groups. See Clarke and Gorley (2001) for a guide to interpreting ANOSIM R values.

veg.adonis <- adonis(shrub.dist ~ Condicion,  contr.unordered = "contr.sum", data = shrub)
veg.adonis

#Terms added sequentially (first to last)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Condicion  1   0.28621 0.28621  1.1156 0.10036  0.15
#Residuals 10   2.56565 0.25656         0.89964       
#Total     11   2.85186                 1.00000      




##############################################################################
#Gr?ficos de ordenaci?n: Non-metric multidimensional scaling NMDS
ord1 <- metaMDS(shrub, distance = "bray")
ord1
plot(ord1, type = "t")
data.scores <- as.data.frame(scores(ord1))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(ord1)  # create a column of site names, from the rownames of data.scores
head(data.scores)  #look at the data
species.scores <- as.data.frame(scores(ord1, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data
data.scores <- as.data.frame(scores(ord1))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)# create a column of site names, from the rownames of data.scores
data.scores$Condicion <- Condicion  #  add the grp variable created earlier
head(data.scores)  
library(ggplot2)
ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, shape=Condicion, colour=Condicion),size=3) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=SitioA),size=6,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("S" = "red", "P" = "blue")) +
  coord_equal() +
  theme_bw()

#make it nicer

ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Condicion,colour=Condicion),size=4) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Condicion),size=8,vjust=0,hjust=0) +  # add the site labels
  scale_colour_manual(values=c("S" = "red", "P" = "blue")) +
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())


c.s <- data.scores[data.scores$Condicion == "S", ][chull(data.scores[data.scores$Condicion == 
                                                                       "S", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
c.p <- data.scores[data.scores$Condicion== "P", ][chull(data.scores[data.scores$Condicion == 
                                                                      "P", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
hull.data <- rbind(c.s, c.p)  

hull.data


ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Condicion,group=Condicion),alpha=0.30) + # add the convex hulls
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=SitioA),size=4) + # add the point markers
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())


########NMDS ELIPTICO
library(ggplot2)
p1 <-ggplot(data.scores, aes(NMDS1, NMDS2))+
  geom_point(aes(NMDS1, NMDS2, color=Condicion),position=position_jitter(.1))+##separates overlapping points
  stat_ellipse(aes(fill=Condicion), alpha=.2,type='t',size =1, geom="polygon")+ ##changes shading on ellipses
  theme_minimal()

p1


#el stress indice permite saber si las diferencias y similitudes de los datos
#�estan bien expresadas en el graficos. En este caso si. hay que fijrse q tan dispersos 
#estan los puntos de la linea

stressplot(ord1)

#Clarke 1993 suggests the following guidelines for acceptable stress values: 
#<0.05 = excellent, <0.10 = good, <0.20 = usable, >0.20 = not acceptable. The plot shows the border of the 0.20 stress value limit. Solutions with higher 
#stress values should be interpreted with caution and those with stress above 0.30 are highly suspect.

metaMDS(
  data,
  distance = "bray",
  k = 2,
  trymax = 20,
  autotransform = TRUE) #  0.09775674  el valor de stress es bueno

#################################
require(graphics)
install.packages("ggdendro")
library("ggplot2")
library(ggdendro)
head(data)

a <- hclust(data.dist, "ave")
plot(a)
plot(a, hang = -1) # aca se ve que lo sitios estan re mezclados 
# la condicion no explica ni la variable del suelo ni la composicion de la vegetacion 

ggdendrogram(a)

# Build dendrogram object from hclust results
dend <- as.dendrogram(a)
# Extract the data (for rectangular lines)
# Type can be "rectangle" or "triangle"
dend_data <- dendro_data(dend, type = "rectangle")
# What contains dend_data
names(dend_data)
head(dend_data$segments)
head(dend_data$labels)
#########################################################################
#######################################################################
##########################################################################
#############COMPOSICION DE HERBACEAS############
herb<- read.table(file.choose(), header = TRUE)
head(herb) # base de datos veg_herb


#Analysis of Similarities
data.dist.h<-vegdist(herb, "bray") # para datos de abundancia se usa bray, la primer celda de la base de datos debe estar VACIA 
data.dist.h
summary(data.dist.h)

sitios<- read.table(file.choose(), header = TRUE)#cargo la info de los sitios. Archivo"Sitios.txt" o"sitios+alt+DAS"
attach(sitios)

dispersion<-betadisper(data.dist.h, group=Condicion)
permutest(dispersion) 
mod.HSD <- TukeyHSD(dispersion) # no da significativa puedo seguir 
plot(mod.HSD)

#La matriz.de.distancias es la que obtienes con vegdist
#El grupo en tu caso sería la condic
data.ano <- anosim(data.dist.h, Condicion)
summary(data.ano)

#Dissimilarity ranks between and within classes:
#  0%   25%  50%  75% 100%  N
#Between  1 18.75 36.5 49.5   65 36
#P        5 14.50 25.0 32.5   47 15
#S        2 28.00 50.0 61.0   66 15

plot(data.ano)
#R=0.019. P=0.36

#The ANOSIM statistic compares the mean of ranked dissimilarities between groups to 
#the mean of ranked dissimilarities within groups. An R value close to "1.0" suggests
#dissimilarity between groups while an R value close to "0" suggests an even distribution 
#of high and low ranks within and between groups. R values below "0" suggest that dissimilarities are greater within groups
#than between groups. See Clarke and Gorley (2001) for a guide to interpreting ANOSIM R values.

veg.adonis <- adonis(data.dist.h ~ Condicion,  contr.unordered = "contr.sum", data = herb)
veg.adonis

#          Df SumsOfSqs MeanSqs F.Model   R2 Pr(>F)
#Condicion  1   0.28621 0.28621  1.1156 0.10036   0.36
#Residuals 10   2.56565 0.25656         0.89964       
#Total     11   2.85186                 1.00000  


##############################################################################
#Gr?ficos de ordenaci?n: Non-metric multidimensional scaling NMDS
ord1 <- metaMDS(herb, distance = "bray")
ord1
plot(ord1, type = "t")
data.scores <- as.data.frame(scores(ord1))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(ord1)  # create a column of site names, from the rownames of data.scores
head(data.scores)  #look at the data
species.scores <- as.data.frame(scores(ord1, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data
data.scores <- as.data.frame(scores(ord1))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)# create a column of site names, from the rownames of data.scores
data.scores$Condicion <- Condicion  #  add the grp variable created earlier
head(data.scores)

########NMDS ELIPTICO
library(ggplot2)
p1 <-ggplot(data.scores, aes(NMDS1, NMDS2))+
  geom_point(aes(NMDS1, NMDS2, color=Condicion),position=position_jitter(.1))+##separates overlapping points
  stat_ellipse(aes(fill=Condicion), alpha=.2,type='t',size =1, geom="polygon")+ ##changes shading on ellipses
  theme_minimal()

p1



library(ggplot2)
ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, shape=Condicion, colour=Condicion),size=3) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=SitioA),size=6,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("S" = "red", "P" = "blue")) +
  coord_equal() +
  theme_bw()

#make it nicer

ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Condicion,colour=Condicion),size=4) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Condicion),size=8,vjust=0,hjust=0) +  # add the site labels
  scale_colour_manual(values=c("S" = "red", "P" = "blue")) +
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())


c.s <- data.scores[data.scores$Condicion == "S", ][chull(data.scores[data.scores$Condicion == 
                                                                       "S", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
c.p <- data.scores[data.scores$Condicion== "P", ][chull(data.scores[data.scores$Condicion == 
                                                                      "P", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
hull.data <- rbind(c.s, c.p)  

hull.data


ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Condicion,group=Condicion),alpha=0.30) + # add the convex hulls
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=SitioA),size=4) + # add the point markers
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

########NMDS ELIPTICO
library(ggplot2)
p1 <-ggplot(data.scores, aes(NMDS1, NMDS2))+
  geom_point(aes(NMDS1, NMDS2, color=Condicion),position=position_jitter(.1))+##separates overlapping points
  stat_ellipse(aes(fill=Condicion), alpha=.2,type='t',size =1, geom="polygon")+ ##changes shading on ellipses
  theme_minimal()

p1



#el stress indice permite saber si las diferencias y similitudes de los datos
#�estan bien expresadas en el graficos. En este caso si. hay que fijrse q tan dispersos 
#estan los puntos de la linea

stressplot(ord1)

#Clarke 1993 suggests the following guidelines for acceptable stress values: 
#<0.05 = excellent, <0.10 = good, <0.20 = usable, >0.20 = not acceptable. The plot shows the border of the 0.20 stress value limit. Solutions with higher 
#stress values should be interpreted with caution and those with stress above 0.30 are highly suspect.

metaMDS(
  data,
  distance = "bray",
  k = 2,
  trymax = 20,
  autotransform = TRUE) #  0.09775674  el valor de stress es bueno




#################################
require(graphics)
install.packages("ggdendro")
library("ggplot2")
library(ggdendro)
head(data)


a <- hclust(data.dist, "ave") #########mepa q es este 
a <- hclust(dist(data), "ave")


a <- hclust(data.dist, "ave")
plot(a)
plot(a, hang = -1) # aca se ve que lo sitios estan re mezclados 
# la condicion no explica ni la variable del suelo ni la composicion de la vegetacion 

ggdendrogram(a)

# Build dendrogram object from hclust results
dend <- as.dendrogram(a)
# Extract the data (for rectangular lines)
# Type can be "rectangle" or "triangle"
dend_data <- dendro_data(dend, type = "rectangle")
# What contains dend_data
names(dend_data)
head(dend_data$segments)
head(dend_data$labels)













