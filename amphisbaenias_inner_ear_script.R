#library import
library(Morpho)
library(Rvcg)
library(rgl)
library(geomorph)
library(base)
library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(shapes)
#Subsampling the phylogenic tree to match my sampling
setwd("C:/Users/nande/OneDrive - UT Arlington/Desktop/Work_in_progress/Amphisbaenians_inner ear/")
tree <- read.nexus("tree.nex")
class(tree)
plot(tree)
str(tree)
tree <- drop.tip(tree, tip = ("Latastia_longicaudata"))
#read tps file
landmarks <- readland.tps(file = "./aligned_coordinates.tps", specID = "ID")

#Define landmarks and semilandmarks
c01 <- define.sliders(landmarks, nsliders=14)
c02 <- define.sliders(landmarks, nsliders=13)
c03 <- define.sliders(landmarks, nsliders=13)

curves <- rbind(c01, c02, c03)
write.csv(curves,"C:./curves.csv")

curves <- read.csv("C:./curves.csv", row.names = 1)
##Conducting GPA
aligned_inner_ear <- gpagen(landmarks, curves = curves, PrinAxes = FALSE, print.progress = TRUE)
class(aligned_inner_ear)

#save aligned array
save(aligned_inner_ear, file="C:./aligned_inner_ear_gpagen.RData")
shapes3d(aligned_inner_ear$coords[,,2],type="p", color = 4,rglopen=TRUE)
rgl.snapshot('3dplot.png', fmt = 'png')

# * 6 PCA ---------------------------------------------
families <- read.csv(file="C:./families.csv",row.names = 1, header = TRUE, stringsAsFactors = TRUE)



PCA.w.phylo <- gm.prcomp(aligned_inner_ear$coords, phy= tree)
summary(PCA.w.phylo)
View(PCA.w.phylo)

pc_scores<- as.data.frame(PCA.w.phylo[[10]])
pc_scores
pc_scores1 <- cbind(families, pc_scores)  
plot(PCA.w.phylo, phylo=TRUE, axis1 = 2, axis2 = 3)

#PC1vsPC2
svg(file="PC1vsPC2.svg", height=10, width=10)
plot(PCA.w.phylo, phylo=TRUE, phylo.par = list(tip.labels=FALSE, node.labels=FALSE, node.cex=0), axes=FALSE,  axis1 = 1, axis2 = 2)
axis(1)
axis(2)
points(pc_scores1$Comp1[pc_scores1$Family=="Amphisbaenidae"], pc_scores1$Comp2[pc_scores1$Family=="Amphisbaenidae"], pch=21, cex=1.9, lwd= 1.6,bg="#666e6e", col= "black")
points(pc_scores1$Comp1[pc_scores1$Family=="Trogonophidae"], pc_scores1$Comp2[pc_scores1$Family=="Trogonophidae"], pch=21, cex=1.9, lwd= 1.6, bg="#38ba9a",col= "black" )
points(pc_scores1$Comp1[pc_scores1$Family=="Cadidae"], pc_scores1$Comp2[pc_scores1$Family=="Cadidae"], pch=21, cex=1.9, lwd= 1.6, bg="#ee799a",col= "black"  )
points(pc_scores1$Comp1[pc_scores1$Family=="Bipedidae"], pc_scores1$Comp2[pc_scores1$Family=="Bipedidae"], pch=21, cex=1.9, lwd= 1.6, bg="#f08965", col= "black")
points(pc_scores1$Comp1[pc_scores1$Family=="Blanidae"], pc_scores1$Comp2[pc_scores1$Family=="Blanidae"], pch=21, cex=1.9, lwd= 1.6, bg="#caae39", col= "black")
points(pc_scores1$Comp1[pc_scores1$Family=="Rhineuridae"], pc_scores1$Comp2[pc_scores1$Family=="Rhineuridae"], pch=21, cex=1.9, lwd= 1.6, bg="#74a8da", col= "black")
legend("bottomright", cex=1,
       legend = c("Amphisbaenidae","Bipedidae", "Blanidae", "Cadidae", "Rhineuridae","Trogonophidae"), 
       col="black",
       pch = c(21, 21, 21, 21, 21, 21), bg=NULL,
       pt.bg=c("#666e6e", "#f08965","#caae39", "#ee799a","#74a8da", "#38ba9a"  ),
       pt.cex=2.2, y.intersp = 1,
       xpd=FALSE, bty="n")
dev.off()
#PC2vsPC3
svg(file="PC2vsPC3.svg", height=10, width=10)
plot(PCA.w.phylo, phylo=TRUE, phylo.par = list(tip.labels=FALSE, node.labels=FALSE, node.cex=0), axes=FALSE,  axis1 = 2, axis2 = 3)
axis(1)
axis(2)
points(pc_scores1$Comp2[pc_scores1$Family=="Amphisbaenidae"], pc_scores1$Comp3[pc_scores1$Family=="Amphisbaenidae"], pch=21, cex=1.9, lwd= 1.6,bg="#666e6e", col= "black")
points(pc_scores1$Comp2[pc_scores1$Family=="Trogonophidae"], pc_scores1$Comp3[pc_scores1$Family=="Trogonophidae"], pch=21, cex=1.9, lwd= 1.6, bg="#38ba9a",col= "black" )
points(pc_scores1$Comp2[pc_scores1$Family=="Cadidae"], pc_scores1$Comp3[pc_scores1$Family=="Cadidae"], pch=21, cex=1.9, lwd= 1.6, bg="#ee799a",col= "black"  )
points(pc_scores1$Comp2[pc_scores1$Family=="Bipedidae"], pc_scores1$Comp3[pc_scores1$Family=="Bipedidae"], pch=21, cex=1.9, lwd= 1.6, bg="#f08965", col= "black")
points(pc_scores1$Comp2[pc_scores1$Family=="Blanidae"], pc_scores1$Comp3[pc_scores1$Family=="Blanidae"], pch=21, cex=1.9, lwd= 1.6, bg="#caae39", col= "black")
points(pc_scores1$Comp2[pc_scores1$Family=="Rhineuridae"], pc_scores1$Comp3[pc_scores1$Family=="Rhineuridae"], pch=21, cex=1.9, lwd= 1.6, bg="#74a8da", col= "black")
legend("bottomright", cex=1,
       legend = c("Amphisbaenidae","Bipedidae", "Blanidae", "Cadidae", "Rhineuridae","Trogonophidae"), 
       col="black",
       pch = c(21, 21, 21, 21, 21, 21), bg=NULL,
       pt.bg=c("#666e6e", "#f08965","#caae39", "#ee799a","#74a8da", "#38ba9a"  ),
       pt.cex=2.2, y.intersp = 1,
       xpd=FALSE, bty="n")
dev.off()
shapes3d(PCA.w.phylo$shapes$shapes.comp3$max,type="p", color = 4,rglopen=TRUE)
rgl.snapshot('comp3$max_lateral.png', fmt = 'png')
rgl.snapshot('comp3$max_posterior.png', fmt = 'png')
rgl.snapshot('comp3$max_dorsal.png', fmt = 'png')

#Phylogenetic signal#
physignal(A=aligned_inner_ear$coords, phy=tree, iter = 999, seed = NULL, print.progress = FALSE)

#Allometry
# Test allometry --------------------------------------------

allometry <-procD.lm(aligned_inner_ear$coords~log(aligned_inner_ear$Csize))
summary(allometry)
plotAllometry(allometry, size=aligned_inner_ear$Csize,method="RegScore")

allomtery.phylo <-procD.pgls(aligned_inner_ear$coords~aligned_inner_ear$Csize, phy= tree)
summary(allomtery.phylo)