############################################################################
## The niche of aquatic bdelloid rotifers in the tropics is unclear and different from that of the same species in temperate areas

## Rapeepan Jaturapruek, Diego Fontaneto, Stefano Mammola, Supiyanit Maiphae

############################################################################

## R code to generate the analyses
## Last update: Helsinki, 27 January 2021

############################################################################

# Loading R package --------------------------------------------------------

library("alphahull")
library("Amelia")
library("BAT")
library("car")
library("corrplot")
library("gdata")
library("ggplot2")
library("gridExtra")
library("Hmisc")
library("hypervolume")
library("MuMIn")
library("performance")
library("PerformanceAnalytics")
library("sp")
library("spdep")

# Working directory --------------------------------------------------------

setwd("") #<- Change me

# Loading the database -----------------------------------------------------

db <- read.table(file="Thai_rotifers_db1.csv", sep=',', dec='.', header=T, as.is=F)
str(db)
dim(db)
names(db)
summary(db)

# Analyses -----------------------------------------------------------------

############################################################################
############################################################################
############################################################################

# 1. drivers of occurrence of R. neptunia, R. rotatoria, and R.tardigrada --

# Check the dataset --------------------------------------------------------

# for habitats -------------------------------------------------------------

# number of records/habitat
table(db$HABITAT) # many habitats have too few records...

# merging similar categories
db$hab <- db$HABITAT
levels(db$hab) <- list(Lotic = c("Algae canal SW",
								 "Canal RW",
								 "Canal SW",
								 "River",
								 "Vegetative canal"), 
                       Lentic = c("Lake",
                                  "Pond",
                                  "Reservoir",
                                  "Vegetative pond",
                                  "Rice field"),
                       Swamps = c("Algae swamp",
                                  "Marsh",
                                  "Peat swamp",
                                  "Swamp",
                                  "Vegetative swamp")) 

# number of records/merged habitat
table(db$hab) 

# Lotic Lentic Swamps 
#    42    108    240 

# missing data -------------------------------------------------------------

Amelia::missmap(db)

# remove ELEVATION, with too many missing data
dim(db)
db_reduced <- subset(db, select = -(ELEVATION))
dim(db_reduced)

# distribution of explanatory variables ------------------------------------

names(db_reduced)
Hmisc::hist.data.frame(db_reduced[,c(9:14,26:29)])

# transform in log scale CONDUCTIVITY, TDS, and SALINITY
summary(db_reduced)
db_reduced$CONDUCTIVITY <- log(db_reduced$CONDUCTIVITY)
db_reduced$TDS <- log(db_reduced$TDS)
db_reduced$SALINITY <- log(db_reduced$SALINITY+1)
Hmisc::hist.data.frame(db_reduced[,c(9:14,26:29)])

# Collinearity -------------------------------------------------------------

PerformanceAnalytics::chart.Correlation(db_reduced[,c("PREC",
                              "MEAN_AIR_TEMP",
                              "TEMPERATURE",
                              "SRAD",
                              "TDS",
                              "SALINITY",
                              "DO",
                              "pH",
                              "CONDUCTIVITY")])
# Not collinear:
# SALINITY, DO, pH, TEMPERATURE, MEAN_AIR_TEMP, PREC, SRAD

correlations <- cor(db_reduced[,c("macrura",
                                  "mento",
                                  "megarostris",
                                  "neptunia",
                                  "neptunoida",
                                  "ovata",
                                  "rotatoria",
                                  "tardigrada",
                                  "cf_neptunoida")])
corrplot::corrplot.mixed(correlations, 
							lower = "number", 
							upper = "circle",
							tl.pos = "lt",
							diag = "u")

# No strict association between species

colSums(db_reduced[,c("macrura",
                      "mento",
                      "megarostris",
                      "neptunia",
                      "neptunoida",
                      "ovata",
                      "rotatoria",
                      "tardigrada",
                      "cf_neptunoida")])

# Keep only species with more than 20 occurrences
# mento, neptunia, neptunoida, rotatoria, tardigrada

# Check 0/1 balance --------------------------------------------------------

table(db_reduced$neptunia) 
table(db_reduced$rotatoria)
table(db_reduced$tardigrada)
# all unbalanced towards absences

# Generating a subset database for the analyses ----------------------------

db_analyses <- data.frame(lat = db_reduced$LATITUDE,
                          lon = db_reduced$LONGITUDE,
                          hab =  db_reduced$hab,
                          salinity = db_reduced$SALINITY,
                          DO = db_reduced$DO,
                          Temp = db_reduced$TEMPERATURE,
                          pH = db_reduced$pH,
                          srad = db_reduced$SRAD,
                          air_Temp = db_reduced$MEAN_AIR_TEMP,
                          prec = db_reduced$PREC,
                          men = db_reduced$mento,
                          neptunia = db_reduced$neptunia,
                          noi = db_reduced$neptunoida,
                          rotatoria = db_reduced$rotatoria,
                          tardigrada = db_reduced$tardigrada)

summary(db_analyses)

# removes rows with NA

Amelia::missmap(db_analyses)
dim(db_analyses)
db_analyses_reduced <- db_analyses[complete.cases(db_analyses), ]
Amelia::missmap(db_analyses_reduced)
dim(db_analyses_reduced)

# Run logistic models for the occurrence of each species -------------------

# Rotaria neptunia ---------------------------------------------------------

m_nept <- glm(neptunia ~ hab + 
					 salinity + 
					 DO +
					 Temp +
					 pH +
					 srad +
					 air_Temp +
					 prec +
					 men +
					 noi +
					 rotatoria +
					 tardigrada,
					 data = db_analyses_reduced, 
					 family = binomial(link="cloglog"))

# check model fit
					 
performance::check_model(m_nept)

# check spatial dependence

res_m_nept <- residuals(m_nept)

dev.off()
par(mfrow=c(1,2))
plot(db_analyses_reduced$lat, res_m_nept, 
					main = "Latitude",
					xlab = "latitude",
					ylab = "model residuals")
plot(db_analyses_reduced$lon, res_m_nept, 
					main = "Longitude",
					xlab = "longitude",
					ylab = "model residuals")
dev.off()

coords <- jitter(sp::coordinates(cbind(db_analyses_reduced$lon,
							    db_analyses_reduced$lat)),
							    0.001)
nb <- spdep::tri2nb(coords,row.names=NULL)

spdep::moran.test(res_m_nept, spdep::nb2listw(nb))
# no evidence of spatial dependence of the model

# explore the output of the model

summary(m_nept)
car::Anova(m_nept)

# run multimodel averaging to measure the importance of each variable

options(na.action="na.fail")
m_nept_MuMIn <- MuMIn::model.avg(get.models(dredge(m_nept), 
							cumsum(weight) <= .9999))
summary(m_nept_MuMIn)
sw(m_nept_MuMIn)
options(na.action="na.omit")

# plot significant results

p_n <- ggplot(data=db_analyses_reduced,
				aes(x=as.factor(neptunia), 
					y=prec, 
					fill=neptunia)) +
			geom_violin(alpha=0.4, 
				position=position_dodge(width=.75),
				size=1,
				color="black") +
			geom_boxplot(notch=TRUE,  
				outlier.size=-1, 
				color="black",
				lwd=1.2, 
				alpha=0.7) +
			geom_point(shape=21,
				size=2, 
				position=position_jitterdodge(jitter.width=0.3), 
				color="black",
				alpha=1) +
			theme_bw() +
			theme(legend.position="none")
p_nf <- p_n + xlab("") + 
			ylab("precipitation (mm)") +
			scale_x_discrete(labels=c("absent", "present")) +
			labs(tag="A") +
			theme(axis.text=element_text(size=20), 
				axis.title=element_text(size=20), 
				axis.title.y=element_text(margin=margin(r=20)),
				plot.tag=element_text(size=30))


# Rotaria rotatoria ---------------------------------------------------------

m_rota <- glm(rotatoria ~ hab + 
					 salinity + 
					 DO +
					 Temp +
					 pH +
					 srad +
					 air_Temp +
					 prec +
					 men +
					 noi +
					 neptunia +
					 tardigrada,
					 data = db_analyses_reduced, 
					 family = binomial(link="cloglog"))

# check model fit
					 
performance::check_model(m_rota)

# check spatial dependence

res_m_rota <- residuals(m_rota)

dev.off()
par(mfrow=c(1,2))
plot(db_analyses_reduced$lat, res_m_rota, 
					main = "Latitude",
					xlab = "latitude",
					ylab = "model residuals")
plot(db_analyses_reduced$lon, res_m_rota, 
					main = "Longitude",
					xlab = "longitude",
					ylab = "model residuals")
dev.off()

coords <- jitter(sp::coordinates(cbind(db_analyses_reduced$lon,
							    db_analyses_reduced$lat)),
							    0.001)
nb <- spdep::tri2nb(coords,row.names=NULL)

spdep::moran.test(res_m_rota, spdep::nb2listw(nb))
# no evidence of spatial dependence of the model

# explore the output of the model

summary(m_rota)
car::Anova(m_rota)

# run multimodel averaging to measure the importance of each variable

options(na.action="na.fail")
m_rota_MuMIn <- MuMIn::model.avg(get.models(dredge(m_rota), 
							cumsum(weight) <= .9999))
summary(m_rota_MuMIn)
sw(m_rota_MuMIn)
options(na.action="na.omit")

# plot significant results

p_r <- ggplot(data=db_analyses_reduced,
				aes(x=as.factor(rotatoria), 
					y=pH, 
					fill=rotatoria)) +
			geom_violin(alpha=0.4, 
				position=position_dodge(width=.75),
				size=1,
				color="black") +
			geom_boxplot(notch=TRUE,  
				outlier.size=-1, 
				color="black",
				lwd=1.2, 
				alpha=0.7) +
			geom_point(shape=21,
				size=2, 
				position=position_jitterdodge(jitter.width=0.3), 
				color="black",
				alpha=1) +
			theme_bw() +
			theme(legend.position="none")
p_rf <- p_r + xlab("") + 
			ylab("pH") +
			scale_x_discrete(labels=c("absent", "present")) +
			labs(tag="B") +
			theme(axis.text=element_text(size=20), 
				axis.title=element_text(size=20),
				axis.title.y=element_text(margin=margin(r=20)),
				plot.tag=element_text(size=30, margin=margin(l=40)))


# Rotaria tardigrada -------------------------------------------------------

m_tard <- glm(tardigrada ~ hab + 
					 salinity + 
					 DO +
					 Temp +
					 pH +
					 srad +
					 air_Temp +
					 prec +
					 men +
					 noi +
					 rotatoria +
					 neptunia,
					 data = db_analyses_reduced, 
					 family = binomial(link="cloglog"))

# check model fit
					 
performance::check_model(m_tard)

# check spatial dependence

res_m_tard <- residuals(m_tard)

dev.off()
par(mfrow=c(1,2))
plot(db_analyses_reduced$lat, res_m_tard, 
					main = "Latitude",
					xlab = "latitude",
					ylab = "model residuals")
plot(db_analyses_reduced$lon, res_m_tard, 
					main = "Longitude",
					xlab = "longitude",
					ylab = "model residuals")
dev.off()

coords <- jitter(sp::coordinates(cbind(db_analyses_reduced$lon,
							    db_analyses_reduced$lat)),
							    0.001)
nb <- spdep::tri2nb(coords,row.names=NULL)

spdep::moran.test(res_m_tard, spdep::nb2listw(nb))
# no evidence of spatial dependence of the model

# explore the output of the model

summary(m_tard)
car::Anova(m_tard)

# run multimodel averaging to measure the importance of each variable

options(na.action="na.fail")
m_tard_MuMIn <- MuMIn::model.avg(get.models(dredge(m_tard), 
						cumsum(weight) <= .9999))
summary(m_tard_MuMIn)
sw(m_tard_MuMIn)
options(na.action="na.omit")

# prepare the final figure -------------------------------------------------
pdf("Figure1.pdf",width=14, height=7, paper='special')
gridExtra::grid.arrange(p_nf, p_rf,ncol=2)
dev.off()

############################################################################
############################################################################
############################################################################

# 2. Niche differentiation among species -----------------------------------

# Rearranging the database -------------------------------------------------

for(i in c(17:24)){
  
  if (i == 17) {
    
  db2 <- db_reduced[db_reduced[,i]  == 1,] ; db2 <- droplevels(db2)
  
  db_niche <- data.frame(species = as.factor(rep(colnames(db_reduced)[i], nrow(db2))),
                            salinity = db2$SALINITY,
                            DO = db2$DO,
                            Temp = db2$TEMPERATURE,
                            pH = db2$pH,
                            srad = db2$SRAD,
                            air_Temp = db2$MEAN_AIR_TEMP,
                            prec  = db2$PREC)
  
  } else {
    
    db2 <- db_reduced[db_reduced[,i]  == 1,] ; db2 <- droplevels(db2)
    
    db_niche2 <- data.frame(species = as.factor(rep(colnames(db_reduced)[i], nrow(db2))),
                           salinity = db2$SALINITY,
                           DO = db2$DO,
                           Temp = db2$TEMPERATURE,
                           pH = db2$pH,
                           srad = db2$SRAD,
                           air_Temp = db2$MEAN_AIR_TEMP,
                           prec  = db2$PREC)
    
    db_niche <- rbind(db_niche, db_niche2)
  }
}

# Rescaling variables ------------------------------------------------------

# Scaling all variable to achieve same dimensionality for hypervolume construction
for(i in 2:ncol(db_niche))
  db_niche[,i] <-  scale(db_niche[,i])

## Do we have enough observation?
table(db_niche$species)
db_niche <- droplevels(db_niche[-which(db_niche$species == "cf_neptunoida"), ] )

#Rename
levels(db_niche$species) <- paste(rep("Rotaria",length(levels(db_niche$species))), as.character(levels(db_niche$species)), sep=' ')

# Generating the hypervolumes in a loop ------------------------------------

for(i in 1:nlevels(db_niche$species)){
  
  #Subsetting the i-species
  db_sp <-  db_niche[db_niche$species == levels(db_niche$species)[i],] 
  
  newHv <- hypervolume::hypervolume_gaussian(na.omit(db_sp)[,c(2:8)],  
  			name= levels(db_niche$species)[i],
            kde.bandwidth = 2*estimate_bandwidth(na.omit(db_sp)[,c(2:8)]),
            verbose=FALSE)
  
  if(i == 1){
    hvBoth <- newHv
    cat(paste("Hypervolume #",as.character(i)," out of ",nlevels(db_niche$species)," has been constructed.",sep=''))
  }
  else{
    hvBoth  <- hypervolume_join(hvBoth,newHv)
    cat(paste("\nHypervolume #",as.character(i)," out of ",nlevels(db_niche$species)," has been constructed.",sep=''))
  }
}

# Visualizing all species niches -------------------------------------------
pdf("FigureS1.pdf",width=14, height=9, paper='special')
plot(hvBoth,
     num.points.max.random = 2000,
     show.data=F,
     pairplot=T,
     show.3d=F,
     showdensity=T,
     showrandom=F,
     contour.lwd=1.5,
     cex.names=1.5,
     cex.random=0.3,
     cex.legend=1,
     cex.axis=1,
     col = c("brown4","blue","black","cadetblue4","firebrick4","dodgerblue3","darkmagenta"),
     contour.filled=T,
     contour.filled.alpha=1,
     show.centroid=TRUE,
     cex.centroid=2,
     names= c("Salinity", "Dissolved oxygen", "Temperature","pH","Solar radiation","Air temperature","Precipitation"))
dev.off()

# Selecting the hypervolumes for 3 species target --------------------------

SubSp <- hypervolume_join(hvBoth[[3]],hvBoth[[6]],hvBoth[[7]])

SubSp@HVList[[1]]@Name <- "Rotaria neptunia"
SubSp@HVList[[2]]@Name <- "Rotaria rotatoria"
SubSp@HVList[[3]]@Name <- "Rotaria tardigrada"

# Estimating hypervolume' statistics with BAT ------------------------------

round(BAT::kernel.alpha(SubSp),0) #Volume

BAT::kernel.dispersion(SubSp) #Dispersion

(Beta <- BAT::kernel.beta(SubSp)) #Niche differentiation

(DistanceCentroids <- BAT::kernel.similarity(SubSp)$Distance_centroids) #Centroid distance

# Prepare the final figure -------------------------------------------------

pdf("Figure2.pdf",width=13, height=11, paper='special')
plot(SubSp,
     num.points.max.random = 4000,
     show.data=F,
     pairplot=T,
     show.3d=F,
     showdensity=T,
     showrandom=F,
     contour.lwd=1.5,
     col= c("black","dodgerblue3","darkmagenta"),
     cex.names=1.2,
     cex.random=0.3,
     cex.legend=1.5,
     cex.axis=1,
     contour.filled=T,
     contour.filled.alpha=1,
     show.centroid=TRUE,
     cex.centroid=2,
     names= c("Salinity", "Dissolved oxygen", "Temperature","pH","Solar radiation","Air temperature","Precipitation"))
dev.off()




############################################################################
# repeat analyses
# with the inclusion of two variables for tropic status
# Chla and Phycocyanin_BGA.PC,
# which are missing for several sites
############################################################################

# Loading the database -----------------------------------------------------

db2 <- read.table(file="Thai_rotifers_db2.csv", sep=',', dec='.', header=T, as.is=F)
dim(db2)
names(db2)
summary(db)

# Analyses -----------------------------------------------------------------

############################################################################
############################################################################
############################################################################

# 1. drivers of occurrence of R. neptunia, R. rotatoria, and R.tardigrada --

# Check the dataset --------------------------------------------------------

# for habitats -------------------------------------------------------------

# number of records/habitat
table(db2$HABITAT) # many habitats have too few records...

# merging similar categories
db2$hab <- db2$HABITAT
levels(db2$hab) <- list(Lotic = c("Algae canal SW",
								 "Canal RW",
								 "Canal SW",
								 "River",
								 "Vegetative canal"), 
                       Lentic = c("Lake",
                                  "Pond",
                                  "Reservoir",
                                  "Vegetative pond",
                                  "Rice field"),
                       Swamps = c("Algae swamp",
                                  "Marsh",
                                  "Peat swamp",
                                  "Swamp",
                                  "Vegetative swamp")) 

# number of records/merged habitat
table(db2$hab) 

# Lotic Lentic Swamps 
#    42    108    240 

# missing data -------------------------------------------------------------

Amelia::missmap(db2)

# remove ELEVATION, with too many missing data
dim(db2)
db2_reduced <- subset(db2, select = -(ELEVATION))
dim(db2_reduced)

# removes rows with NA in Chla

Amelia::missmap(db2_reduced)
dim(db2_reduced)
names(db2_reduced)
db2_reduced <- db2_reduced[complete.cases(db2_reduced[,30]), ]
dim(db2_reduced)
Amelia::missmap(db2_reduced)

# removes STATION and the row with no coordinates

db2_reduced <- subset(db2_reduced, select = -(STATION))
dim(db2_reduced)
db2_reduced <- db2_reduced[complete.cases(db2_reduced[,7]), ]
dim(db2_reduced)
Amelia::missmap(db2_reduced)

# distribution of explanatory variables ------------------------------------

names(db2_reduced)
dev.off()
Hmisc::hist.data.frame(db2_reduced[,c(8:13,25:30)])
dev.off()

# transform in log scale CONDUCTIVITY, TDS, SALINITY, and Chla
summary(db2_reduced)
db2_reduced$CONDUCTIVITY <- log(db2_reduced$CONDUCTIVITY)
db2_reduced$TDS <- log(db2_reduced$TDS)
db2_reduced$SALINITY <- log(db2_reduced$SALINITY)
db2_reduced$Chla <- log(db2_reduced$Chla)
Hmisc::hist.data.frame(db2_reduced[,c(8:13,25:30)])

# Collinearity -------------------------------------------------------------

PerformanceAnalytics::chart.Correlation(db2_reduced[,c("PREC",
                              "MEAN_AIR_TEMP",
                              "TEMPERATURE",
                              "SRAD",
                              "TDS",
                              "SALINITY",
                              "DO",
                              "pH",
                              "CONDUCTIVITY",
                              "Chla",
                              "Phycocyanin_BGA.PC")])
# Not collinear:
# SALINITY, pH, TEMPERATURE, MEAN_AIR_TEMP, 
# and Chla, Phycocyanin_BGA.PC

correlations <- cor(db2_reduced[,c("macrura",
                                  "mento",
                                  "megarostris",
                                  "neptunia",
                                  "neptunoida",
                                  "ovata",
                                  "rotatoria",
                                  "tardigrada")])
corrplot::corrplot.mixed(correlations, 
							lower = "number", 
							upper = "circle",
							tl.pos = "lt",
							diag = "u")

# No strict association between species

colSums(db2_reduced[,c("macrura",
                      "mento",
                      "megarostris",
                      "neptunia",
                      "neptunoida",
                      "ovata",
                      "rotatoria",
                      "tardigrada")])

# Keep only species with more than 5 occurrences
# mento, neptunia, neptunoida, rotatoria, tardigrada

# Check 0/1 balance --------------------------------------------------------

table(db2_reduced$neptunia) 
table(db2_reduced$rotatoria)
table(db2_reduced$tardigrada)
# all unbalanced towards absences

# Generating a subset database for the analyses ----------------------------

db2_analyses <- data.frame(lat = db2_reduced$LATITUDE,
                          lon = db2_reduced$LONGITUDE,
                          hab =  db2_reduced$hab,
                          salinity = db2_reduced$SALINITY,
                          Temp = db2_reduced$TEMPERATURE,
                          pH = db2_reduced$pH,
                          air_Temp = db2_reduced$MEAN_AIR_TEMP,
                          chla = db2_reduced$Chla,
                          phyco = db2_reduced$Phycocyanin_BGA.PC,
                          men = db2_reduced$mento,
                          neptunia = db2_reduced$neptunia,
                          noi = db2_reduced$neptunoida,
                          rotatoria = db2_reduced$rotatoria,
                          tardigrada = db2_reduced$tardigrada)

db2_analyses <- gdata::drop.levels(db2_analyses)

summary(db2_analyses)
names(db2_analyses)
dim(db2_analyses)

PerformanceAnalytics::chart.Correlation(db2_analyses[,4:14])
# remove salinity from models, due to high correlation

# Run logistic models for the occurrence of each species -------------------

# Rotaria neptunia ---------------------------------------------------------

m_nept <- glm(neptunia ~ hab + 
					 Temp +
					 pH +
					 air_Temp +
					 chla +
					 phyco +
					 men +
					 noi +
					 rotatoria +
					 tardigrada,
					 data = db2_analyses, 
					 family = binomial(link="cloglog"))

# check model fit
					 
performance::check_model(m_nept)

# check spatial dependence

res_m_nept <- residuals(m_nept)

dev.off()
par(mfrow=c(1,2))
plot(db2_analyses$lat, res_m_nept, 
					main = "Latitude",
					xlab = "latitude",
					ylab = "model residuals")
plot(db2_analyses$lon, res_m_nept, 
					main = "Longitude",
					xlab = "longitude",
					ylab = "model residuals")
dev.off()

coords <- jitter(sp::coordinates(cbind(db2_analyses$lon,
							    db2_analyses$lat)),
							    0.001)
nb <- spdep::tri2nb(coords,row.names=NULL)

spdep::moran.test(res_m_nept, spdep::nb2listw(nb))
# no evidence of spatial dependence of the model

# explore the output of the model

summary(m_nept)
car::Anova(m_nept)

# run multimodel averaging to measure the importance of each variable

options(na.action="na.fail")
m_nept_MuMIn <- MuMIn::model.avg(get.models(dredge(m_nept), 
							cumsum(weight) <= .9999))
sw(m_nept_MuMIn)
options(na.action="na.omit")

# plot significant results

p_n3 <- ggplot(data=db2_analyses,
				aes(x=as.factor(neptunia), 
					y=phyco, 
					fill=neptunia)) +
			geom_violin(alpha=0.4, 
				position=position_dodge(width=.75),
				size=1,
				color="black") +
			geom_boxplot(notch=TRUE,  
				outlier.size=-1, 
				color="black",
				lwd=1.2, 
				alpha=0.7) +
			geom_point(shape=21,
				size=2, 
				position=position_jitterdodge(jitter.width=0.3), 
				color="black",
				alpha=1) +
			theme_bw() +
			theme(legend.position="none")
p_nf3 <- p_n3 + xlab("") + 
			ylab("Phycocyanin") +
			scale_x_discrete(labels=c("absent", "present")) +
			labs(tag="B") +
			theme(axis.text=element_text(size=20), 
				axis.title=element_text(size=20), 
				axis.title.y=element_text(margin=margin(r=20)),
				plot.tag=element_text(size=30))
p_nf3 

p_n2 <- ggplot(data=db2_analyses,
				aes(x=as.factor(neptunia), 
					y=exp(chla), 
					fill=neptunia)) +
			geom_violin(alpha=0.4, 
				position=position_dodge(width=.75),
				size=1,
				color="black") +
			geom_boxplot(notch=TRUE,  
				outlier.size=-1, 
				color="black",
				lwd=1.2, 
				alpha=0.7) +
			geom_point(shape=21,
				size=2, 
				position=position_jitterdodge(jitter.width=0.3), 
				color="black",
				alpha=1) +
			theme_bw() +
			theme(legend.position="none")
p_nf2 <- p_n2 + xlab("") + 
			ylab("Chlorophyll a") +
			scale_x_discrete(labels=c("absent", "present")) +
			labs(tag="A") +
			theme(axis.text=element_text(size=20), 
				axis.title=element_text(size=20), 
				axis.title.y=element_text(margin=margin(r=20)),
				plot.tag=element_text(size=30))
p_nf2 


# Rotaria rotatoria ---------------------------------------------------------

m_rota <- glm(rotatoria ~ hab + 
					 Temp +
					 pH +
					 air_Temp +
					 chla +
					 phyco +
					 men +
					 noi +
					 neptunia +
					 tardigrada,
					 data = db2_analyses, 
					 family = binomial(link="cloglog"))

# check model fit
					 
performance::check_model(m_rota)

# check spatial dependence

res_m_rota <- residuals(m_rota)

dev.off()
par(mfrow=c(1,2))
plot(db2_analyses$lat, res_m_rota, 
					main = "Latitude",
					xlab = "latitude",
					ylab = "model residuals")
plot(db2_analyses$lon, res_m_rota, 
					main = "Longitude",
					xlab = "longitude",
					ylab = "model residuals")
dev.off()

coords <- jitter(sp::coordinates(cbind(db2_analyses$lon,
							    db2_analyses$lat)),
							    0.001)
nb <- spdep::tri2nb(coords,row.names=NULL)

spdep::moran.test(res_m_rota, spdep::nb2listw(nb))
# no evidence of spatial dependence of the model

# explore the output of the model

summary(m_rota)
car::Anova(m_rota)

# run multimodel averaging to measure the importance of each variable

options(na.action="na.fail")
m_rota_MuMIn <- MuMIn::model.avg(get.models(dredge(m_rota), 
							cumsum(weight) <= .9999))
sw(m_rota_MuMIn)
options(na.action="na.omit")


# Rotaria tardigrada -------------------------------------------------------

m_tard <- glm(tardigrada ~ hab + 
					 Temp +
					 pH +
					 air_Temp +
					 chla +
					 phyco +
					 men +
					 noi +
					 neptunia +
					 rotatoria,
					 data = db2_analyses, 
					 family = binomial(link="cloglog"))

# check model fit
					 
performance::check_model(m_tard)

# check spatial dependence

res_m_tard <- residuals(m_tard)

dev.off()
par(mfrow=c(1,2))
plot(db2_analyses$lat, res_m_tard, 
					main = "Latitude",
					xlab = "latitude",
					ylab = "model residuals")
plot(db2_analyses$lon, res_m_tard, 
					main = "Longitude",
					xlab = "longitude",
					ylab = "model residuals")
dev.off()

coords <- jitter(sp::coordinates(cbind(db2_analyses$lon,
							    db2_analyses$lat)),
							    0.001)
nb <- spdep::tri2nb(coords,row.names=NULL)

spdep::moran.test(res_m_tard, spdep::nb2listw(nb))
# no evidence of spatial dependence of the model

# explore the output of the model

summary(m_tard)
car::Anova(m_tard)

# run multimodel averaging to measure the importance of each variable

options(na.action="na.fail")
m_tard_MuMIn <- MuMIn::model.avg(get.models(dredge(m_tard), 
						cumsum(weight) <= .9999))
sw(m_tard_MuMIn)
options(na.action="na.omit")

# prepare the final figure -------------------------------------------------
pdf("Figure try.pdf",width=14, height=7, paper='special')
gridExtra::grid.arrange(p_nf3, p_nf2,ncol=2)
dev.off()

pdf("Figure1 new.pdf",width=14, height=7, paper='special')
gridExtra::grid.arrange(p_nf, p_nf3,ncol=2)
dev.off()


############################################################################
############################################################################
############################################################################

# 2. Niche differentiation among species -----------------------------------

# Rearranging the database ------------------------------------------------

for(i in c(16:22)){
  
  if (i == 16) {
    
  db2 <- db2_reduced[db2_reduced[,i]  == 1,] ; db2 <- droplevels(db2)
  
  db_niche <- data.frame(species = as.factor(rep(colnames(db2_reduced)[i], nrow(db2))),
                            Temp = db2$TEMPERATURE,
                            pH = db2$pH,
                            air_Temp = db2$MEAN_AIR_TEMP,
                            chla  = db2$Chla,
                            phyco = db2$Phycocyanin_BGA.PC)
  
  } else {
    
    db2 <- db2_reduced[db2_reduced[,i]  == 1,] ; db2 <- droplevels(db2)
    
    db_niche2 <- data.frame(species = as.factor(rep(colnames(db2_reduced)[i], nrow(db2))),
                            Temp = db2$TEMPERATURE,
                            pH = db2$pH,
                            air_Temp = db2$MEAN_AIR_TEMP,
                            chla  = db2$Chla,
                            phyco = db2$Phycocyanin_BGA.PC)
    
    db_niche <- rbind(db_niche, db_niche2)
  }
}

# Rescaling variables -----------------------------------------------------

# Scaling all variable to achieve same dimensionality for hypervolume construction
for(i in 2:ncol(db_niche))
  db_niche[,i] <-  scale(db_niche[,i])

## Do we have enough observation?
table(db_niche$species)
db_niche <- droplevels(db_niche[-which(db_niche$species == "ovata"), ] )

#Rename
levels(db_niche$species) <- paste(rep("Rotaria",length(levels(db_niche$species))), as.character(levels(db_niche$species)), sep=' ')

# Generating the hypervolumes in a loop -----------------------------------

for(i in 1:nlevels(db_niche$species)){
  
  #Subsetting the i-species
  db_sp <-  db_niche[db_niche$species == levels(db_niche$species)[i],] 
  
  newHv <- hypervolume::hypervolume_gaussian(na.omit(db_sp)[,c(2:6)],  name= levels(db_niche$species)[i],
                                             kde.bandwidth = 2*estimate_bandwidth(na.omit(db_sp)[,c(2:6)]),
                                             verbose=FALSE)
  
  if(i == 1){
    hvBoth <- newHv
    cat(paste("Hypervolume #",as.character(i)," out of ",nlevels(db_niche$species)," has been constructed.",sep=''))
  }
  else{
    hvBoth  <- hypervolume_join(hvBoth,newHv)
    cat(paste("\nHypervolume #",as.character(i)," out of ",nlevels(db_niche$species)," has been constructed.",sep=''))
  }
}

# Visualizing all species niches ------------------------------------------
pdf("FigureS2.pdf",width=14, height=9, paper='special')
plot(hvBoth,
     num.points.max.random = 2000,
     show.data=F,
     pairplot=T,
     show.3d=F,
     showdensity=T,
     showrandom=F,
     contour.lwd=1.5,
     cex.names=1.5,
     cex.random=0.3,
     cex.legend=1,
     cex.axis=1,
     col = c("brown4","blue","black","cadetblue4","dodgerblue3","darkmagenta"),
     contour.filled=T,
     contour.filled.alpha=1,
     show.centroid=TRUE,
     cex.centroid=2,
     names= c("Temperature","pH","Air temperature","Chlorophyll_a", "Phycocyanin"))
dev.off()

# Selecting the hypervolumes for 3 species target -------------------------

SubSp <- hypervolume_join(hvBoth[[3]],hvBoth[[5]],hvBoth[[6]])

SubSp@HVList[[1]]@Name <- "Rotaria neptunia"
SubSp@HVList[[2]]@Name <- "Rotaria rotatoria"
SubSp@HVList[[3]]@Name <- "Rotaria tardigrada"

# Estimating hypervolume' statistics with BAT -----------------------------

round(BAT::kernel.alpha(SubSp),0) #Volume

BAT::kernel.dispersion(SubSp) #Dispersion

(Beta <- BAT::kernel.beta(SubSp)) #Niche differentiation

(DistanceCentroids <- BAT::kernel.similarity(SubSp)$Distance_centroids) #Centroid distance

# Prepare the final figure -------------------------------------------------

pdf("Figure3.pdf",width=13, height=11, paper='special')
plot(SubSp,
     num.points.max.random = 4000,
     show.data=F,
     pairplot=T,
     show.3d=F,
     showdensity=T,
     showrandom=F,
     contour.lwd=1.5,
     col= c("black","dodgerblue3","darkmagenta"),
     cex.names=1.2,
     cex.random=0.3,
     cex.legend=1.5,
     cex.axis=1,
     contour.filled=T,
     contour.filled.alpha=1,
     show.centroid=TRUE,
     cex.centroid=2,
     names= c("Temperature","pH","Air temperature","Chlorophyll_a", "Phycocyanin"))
dev.off()

