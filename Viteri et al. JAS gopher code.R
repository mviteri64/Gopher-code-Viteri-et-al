# LIBRARIES
library(Bchron)

# CHRONOLOGY__________________________________
#Jasper Ridge
JR <- read.csv("~/Jasper_Ridge_dates.csv") 

color <- c("blue", rep("black",2), rep("blue",4), "black", rep("blue",3), rep("black",3),"blue",rep("black",3)) #black=gopher, blue=non-gopher
plot(JR$X14C_age, JR$depth, xlim=c(max(JR$X14C_age),0), ylim=c(max(JR$depth),0), col=color, pch=16)

#Jasper Ridge gophers
JR_gopher <- JR[c(2:3,8,12:14,16:18),]

#JR non-fossorial species
JR_other <- JR[c(1,4:7,9:11,15),]

#plot
plot(JR_gopher[,4],JR_gopher[,6], xlim=c(max(JR$X14C_age),0), ylim=c(max(JR$depth),0),xlab="time", ylab="depth", col="black")
abline(lm(JR_gopher[,6] ~ JR_gopher[,4]), col = "black")
points(JR_other[,4],JR_other[,6], xlim=c(max(JR$X14C_age),0), ylim=c(max(JR$depth),0),col="blue")
abline(lm(JR_other[,6] ~ JR_other[,4]), col = "blue")

#boxplots
par(mfrow=c(1,2))
boxplot(JR_gopher[,4],  ylim=c(2000,0),xlab="gopher",ylab="Age")
boxplot(JR_other[,4], ylim=c(2000,0),xlab="non-gopher", ylab="Age")

#Add significance using wilcox (Mann-Whitney) test (non-parametric, mult. version of wilcox test)

library(ggplot2)
library(ggsignif)

goph_v_non <- c("non-gopher", rep("gopher",2), rep("non-gopher",4), "gopher", rep("non-gopher",3),rep("gopher",3), "non-gopher", rep("gopher",3))
JR_box <- data.frame(goph_v_non, as.numeric(JR[,4]))
colnames(JR_box) <- c("Species", "Age")


ggplot(JR_box, aes(x=Species, y=Age)) + 
  geom_boxplot() +
  scale_y_reverse() +
  geom_signif(comparisons = list(c("gopher", "non-gopher")), 
              map_signif_level=c("***"=0.01, "**"=0.05, "*"=0.1, " "=2))

wilcox.test(JR_gopher$X14C_age, JR_other$X14C_age, paired=FALSE)

top.of.deposit <- 0
bottom.of.deposit <- 110 #for PL - based on greatest level number (17), but all other units stop by 10 or 11...
bottom.of.deposit <- 230 #for SW
bottom.of.deposit <- 110 #for JR
predict.depths <- seq(top.of.deposit,bottom.of.deposit,10)

# calbrate the 14C dates_____________________ Note: positions = depth
calibrated.dates <- BchronCalibrate(ages=PL[,2],ageSds = PL[,3], calCurves=c(rep("intcal20",nrow(PL))), positions = PL[,4], ids=PL[,1] )
calibrated.dates<-BchronCalibrate(ages=SW[,3],ageSds=SW[,4],calCurves=c(rep("intcal20",nrow(SW))),positions=SW[,5],ids=SW[,1])
calibrated.dates<-BchronCalibrate(ages=JR[,4],ageSds=JR[,5],calCurves=c(rep("intcal20",nrow(JR))),positions=JR[,6],ids=JR[,1])
#for JR divide gopher and non-gopher
calibrated.dates.gopher<-BchronCalibrate(ages=JR_gopher[,4],ageSds=JR_gopher[,5],calCurves=c(rep("intcal20",nrow(JR_gopher))),positions=JR_gopher[,6],ids=JR_gopher[,1])
calibrated.dates.other<-BchronCalibrate(ages=JR_other[,4],ageSds=JR_other[,5],calCurves=c(rep("intcal20",nrow(JR_other))),positions=JR_other[,6],ids=JR_other[,1])

#plot calibrated dates
color <- c("blue", rep("black",2), rep("blue",4), "black", rep("blue",3), rep("black",3),"blue",rep("black",2)) #black=gopher, blue=non-gopher
plot(calibrated.dates,withPositions=T,xlim=c(2000,0),ylim=c(200,0),ylab="Depth (cm)",xlab="Cal yr BP",col=color, dateHeight=8,las=1,main="")

#plot calibrated dates for gophers and for non-gophers seperately
par(mfrow=c(1,2))
plot(calibrated.dates.gopher,withPositions=T,xlim=c(2000,0),ylim=c(300,0),ylab="Depth (cm)",xlab="Cal yr BP",dateHeight=8,las=1,main="")
plot(calibrated.dates.other,withPositions=T,xlim=c(2000,0),ylim=c(300,0),ylab="Depth (cm)",xlab="Cal yr BP",dateHeight=8,las=1,main="")

# get credible intervals for each date_____________________
sample_ages<-sampleAges(calibrated.dates)
calibrated.table<-apply(sample_ages,2,quantile,prob=c(0.025,0.5,0.975))
calibrated.table<-t(calibrated.table)
colnames(calibrated.table)<-c("calibrated 2.5%","calibrated 50%","calibrated 97.5%")

chron.table<-cbind(JR,calibrated.table)

chron.goph <- chron.table[c(2:3,8,12:14,16:18),]
chron.nongoph <- chron.table[c(1,4:7,9:11,15),]

#make gopher/non-gopher 50% date histograms
par(mfrow=c(1,2))

chron.goph <- chron.table[c(2:3,8,12:14,16:18),]
chron.nongoph <- chron.table[c(1,4:7,9:11,15),]

hist(chron.goph[,11], breaks=c(0,200, 400, 600, 800, 1000, 1200, 1400, 1600,1800, 2000))
hist(chron.nongoph[,11],breaks=c(0,200, 400, 600, 800, 1000, 1200, 1400, 1600,1800, 2000))

# Make age-depth model for gophers and non-gophers_____________________
am.gopher <- Bchronology(ages=JR_gopher[,4],ageSds=JR_gopher[,5],positions=JR_gopher[,6],positionThicknesses=c(rep(10,nrow(JR))),calCurves=c(rep("intcal20",nrow(JR_gopher))), jitterPositions=T,predictPositions=predict.depths)
am.other <- Bchronology(ages=JR_other[,4],ageSds=JR_other[,5],positions=JR_other[,6],positionThicknesses=c(rep(10,nrow(JR))),calCurves=c(rep("intcal20",nrow(JR_other))),jitterPositions=T,predictPositions=predict.depths)

# get confidence intervals_____________________
am.quants.gopher<-summary(am.gopher)
am.quants.other<-summary(am.other)

# plot age model______________________________
color <- c("blue", rep("black",2), rep("blue",4), "black", rep("blue",3), rep("black",3),"blue",rep("black",3)) #black=gopher, blue=non-gopher

par(mfrow=c(1,2))
plot(am.quants.gopher[,"50%"],am.quants.gopher$Depth,type="l",ylim=c(100,0),xlim=c(1000,0))
polygon(x=c(am.quants.gopher[,"2.5%"],rev(am.quants.gopher[,"97.5%"])),y=c(am.quants.gopher$Depth,rev(am.quants.gopher$Depth)),col=rgb(0,.7,1,.5),border=NA)
polygon(x=c(am.quants.gopher[,"25%"],rev(am.quants.gopher[,"75%"])),y=c(am.quants.gopher$Depth,rev(am.quants.gopher$Depth)),col=rgb(0,.7,1,.9),border=NA)
lines(am.quants.gopher[,"50%"],am.quants.gopher$Depth)

plot(am.quants.other[,"50%"],am.quants.other$Depth,type="l",ylim=c(100,0),xlim=c(1000,0))
polygon(x=c(am.quants.other[,"2.5%"],rev(am.quants.other[,"97.5%"])),y=c(am.quants.other$Depth,rev(am.quants.other$Depth)),col=rgb(0,.7,1,.5),border=NA)
polygon(x=c(am.quants.other[,"25%"],rev(am.quants.other[,"75%"])),y=c(am.quants.other$Depth,rev(am.quants.other$Depth)),col=rgb(0,.7,1,.9),border=NA)
lines(am.quants.other[,"50%"],am.quants.other$Depth)
points(JR_other[,4], JR_other[,6], pch=19)

#plot on same graph
par(mfrow=c(1,1))
plot(am.quants.gopher[,"50%"],am.quants.gopher$Depth,type="l",ylim=c(100,0),xlim=c(1800,0), xlab="Time", ylab="Depth")
polygon(x=c(am.quants.gopher[,"2.5%"],rev(am.quants.gopher[,"97.5%"])),y=c(am.quants.gopher$Depth,rev(am.quants.gopher$Depth)),col="lightgrey",border=NA)
polygon(x=c(am.quants.gopher[,"25%"],rev(am.quants.gopher[,"75%"])),y=c(am.quants.gopher$Depth,rev(am.quants.gopher$Depth)),col="grey",border=NA)
lines(am.quants.gopher[,"50%"],am.quants.gopher$Depth)
lines(am.quants.other[,"50%"],am.quants.other$Depth,type="l",ylim=c(100,0),xlim=c(1800,0))
polygon(x=c(am.quants.other[,"2.5%"],rev(am.quants.other[,"97.5%"])),y=c(am.quants.other$Depth,rev(am.quants.other$Depth)),col=rgb(0,.7,1,.5),border=NA)
polygon(x=c(am.quants.other[,"25%"],rev(am.quants.other[,"75%"])),y=c(am.quants.other$Depth,rev(am.quants.other$Depth)),col=rgb(0,.7,1,.9),border=NA)
lines(am.quants.other[,"50%"],am.quants.other$Depth)
#points(JR[,4], JR[,6], pch=19,col=color) #plotted on non-calibrated dates, fix?
arrows(x0=chron.table[,12], y0=chron.table[,6], x1=chron.table[,10], y1=chron.table[,6], angle=90, code=3,  length=0.05, col=color, lwd=2)
points(chron.table[,11], chron.table[,6], pch=19, cex=1.2, col=color)
text(chron.table[,11], chron.table[,6], labels=chron.table$unit, pos=1, cex=0.9, font=2)
box()
