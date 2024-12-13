##load data

palm1 <- read.table("palm1.txt", header=TRUE, row.names=1, sep="\t")
palm <- palm1[-25,] #outlier in terms of seedling counts, and first transect conducted

comm <- read.table("community3.txt", sep="\t", header=TRUE, row.names=1)
comm.all <- t(comm[,rownames(palm)])

##calculate phylogenetic relatedness to macarthur's palm

library(picante)

phylo <- read.tree("phyloout.txt")
phylo.dist <- cophenetic(phylo)
phylo.dist <- phylo.dist[rownames(phylo.dist)=="Ptychosperma_macarthuri",]
phylo.dist <- phylo.dist[sort(names(phylo.dist))]
phylo.dist <- phylo.dist[names(phylo.dist)!="Ptychosperma_macarthuri"]

palm$phylo.dist <- rep(NA, nrow(palm))
for(i in 1:nrow(palm)){
	palm$phylo.dist[i] <- sum(comm.all[i,]*phylo.dist)/sum(comm.all[i,])
	}

##soil characteristics analysis
palm$pH <- 10^(palm$pH)
palm[,20:23] <- asin(sqrt(palm[,20:23]/100))
palm[,25:26] <- log(palm[,25:26])

soil <- princomp(palm[,19:26], cor=TRUE)

par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(2.5,1,0))

plot(soil$scores[,1],soil$scores[,2], xlab=paste("PCA Axis 1 (", round(soil$sdev[1]^2/8*100, 1), "%)"), ylab=paste("PCA Axis 2 (", round(soil$sdev[2]^2/8*100, 1), "%)"), xlim=c(-4,4), ylim=c(-3,3))
arrows(0,0, soil$loadings[,1]*3.5, soil$loadings[,2]*3.5, length=0.1, col="blue")
text(soil$loadings[,1]*4, soil$loadings[,2]*4, label=rownames(soil$loadings), col="blue")
text(-5.1,3,labels="(a)",xpd=TRUE)

plot(soil$scores[,3],soil$scores[,2], xlab=paste("PCA Axis 3 (", round(soil$sdev[3]^2/8*100, 1), "%)"), ylab=paste("PCA Axis 2 (", round(soil$sdev[2]^2/8*100, 1), "%)"), xlim=c(-4,4), ylim=c(-3,3))
arrows(0,0, soil$loadings[,3]*3.5, soil$loadings[,2]*3.5, length=0.1, col="blue")
text(soil$loadings[,3]*4, soil$loadings[,2]*4, label=rownames(soil$loadings), col="blue")
text(-5.1,3,labels="(b)",xpd=TRUE)

palm$PCA1 <- soil$scores[,1]
palm$PCA2 <- soil$scores[,2]
palm$PCA3 <- soil$scores[,3]

##community composition analysis

library(vegan)

#exclude rare species (occurs in 3 or less transects)

comm.subset <- comm.all[,comm$occurrence>0]
comm.mds <- metaMDS(comm.subset, k=3)
comm.envfit12 <- envfit(comm.mds, cbind(palm[,17:18], palm[,29:31]), choices=c(1,2))
comm.envfit32 <- envfit(comm.mds, cbind(palm[,17:18], palm[,29:31]), choices=c(3,2))

par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(2.5,1,0))

plot(comm.mds$species[,1], comm.mds$species[,2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", type="n")
text(comm.mds$species[,1], comm.mds$species[,2], label=rownames(comm.mds$species), col="red", cex=0.7)
plot(comm.envfit12, arrow.mul=2,col="blue", add=TRUE)
text(-2.6,1.4,labels="(a)",xpd=TRUE)

plot(comm.mds$species[,3], comm.mds$species[,2], xlab="NMDS Axis 3", ylab="NMDS Axis 2", type="n")
text(comm.mds$species[,3], comm.mds$species[,2], label=rownames(comm.mds$species), col="red", cex=0.7)
plot(comm.envfit32, arrow.mul=2, col="blue", add=TRUE)
text(-1.3,1.4,labels="(b)",xpd=TRUE)
text(1,-1,paste("stress = ",round(comm.mds$stress,1)))

palm$NMDS1 <- comm.mds$points[,1]
palm$NMDS2 <- comm.mds$points[,2]
palm$NMDS3 <- comm.mds$points[,3]

##regression tree

library(party)

tree.full <- ctree(seedling~position+ppext.stem+species.richness+phylo.dist+canopy+litter+PCA1+PCA2+PCA3+NMDS1+NMDS2+NMDS3, data=palm, controls=ctree_control(testtype="MonteCarlo",minsplit=2,mincriterion=0.90))
plot(tree.full)

forest.8 <- cforest(seedling~position+ppext.stem+species.richness+phylo.dist+canopy+litter+PCA1+PCA2+PCA3+NMDS1+NMDS2+NMDS3, data=palm, controls=cforest_unbiased(mtry=8, ntree=5000, minsplit=2))
forest.10 <- cforest(seedling~position+ppext.stem+species.richness+phylo.dist+canopy+litter+PCA1+PCA2+PCA3+NMDS1+NMDS2+NMDS3, data=palm, controls=cforest_unbiased(mtry=10, ntree=5000, minsplit=2))
forest.12 <- cforest(seedling~position+ppext.stem+species.richness+phylo.dist+canopy+litter+PCA1+PCA2+PCA3+NMDS1+NMDS2+NMDS3, data=palm, controls=cforest_unbiased(mtry=12, ntree=5000, minsplit=2))

varimp.8 <- varimp(forest.8, conditional=TRUE)
varimp.10 <- varimp(forest.10, conditional=TRUE)
varimp.12 <- varimp(forest.12, conditional=TRUE)

par(mar=c(5,4,2,1), mgp=c(2.5,1,0))
bar<-barplot(varimp.8, axisnames=FALSE, ylim=c(-200,1000), ylab="Conditional Reduction in Variable Importance (Mean Squared Error)")
abline(h=0, lwd=2)
text(bar,-70,adj=1,xpd=TRUE,label=names(varimp.8), srt=50)

#GLMM

library(lme4)

model.full1 <- "seedling~ppext.stem+litter+NMDS1+NMDS2+(1|fragment)"
mod.full1 <- lmer(model.full1, data=palm, family=poisson, REML=FALSE)
c <- sum(residuals(mod.full1)^2)/(length(residuals(mod.full1))-length(mod.full1@fixef)-1)

model.full <- "seedling~ppext.stem+litter+NMDS1+NMDS2+(1|fragment)+(1|rownames(palm))"
mod.full <- lmer(model.full, data=palm, family=poisson, REML=FALSE)

model.null <- "seedling~1+(1|fragment)+(1|rownames(palm))"
model.01 <- "seedling~ppext.stem+(1|fragment)+(1|rownames(palm))"
model.02 <- "seedling~litter+(1|fragment)+(1|rownames(palm))"
model.03 <- "seedling~NMDS1+(1|fragment)+(1|rownames(palm))"
model.04 <- "seedling~NMDS2+(1|fragment)+(1|rownames(palm))"
model.05 <- "seedling~ppext.stem+litter+(1|fragment)+(1|rownames(palm))"
model.06 <- "seedling~ppext.stem+NMDS1+(1|fragment)+(1|rownames(palm))"
model.07 <- "seedling~ppext.stem+NMDS2+(1|fragment)+(1|rownames(palm))"
model.08 <- "seedling~litter+NMDS1+(1|fragment)+(1|rownames(palm))"
model.09 <- "seedling~litter+NMDS2+(1|fragment)+(1|rownames(palm))"
model.10 <- "seedling~NMDS1+NMDS2+(1|fragment)+(1|rownames(palm))"
model.11 <- "seedling~ppext.stem+litter+NMDS1+(1|fragment)+(1|rownames(palm))"
model.12 <- "seedling~ppext.stem+litter+NMDS2+(1|fragment)+(1|rownames(palm))"
model.13 <- "seedling~ppext.stem+NMDS1+NMDS2+(1|fragment)+(1|rownames(palm))"
model.14 <- "seedling~litter+NMDS1+NMDS2+(1|fragment)+(1|rownames(palm))"

model.list<-c(model.null,model.01,model.02,model.03,model.04,model.05,model.06,model.07,model.08,model.09,model.10,model.11,model.12,model.13,model.14,model.full)

k <- rep(NA,length(model.list))
loglik <- rep(NA,length(model.list))
dev <- rep(NA, length(model.list))
agg <- rep(NA, length(model.list))
agg.patch <- rep(NA, length(model.list))
for (i in 1:length(model.list)) {
	mod <- lmer(model.list[i], data=palm, family=poisson)
	k[i] <- length(mod@fixef)
	AIC[i] <- AIC(mod)
	dev[i] <- sum(mod@resid^2)
	vc <- VarCorr(mod)
	agg[i] <- exp(vc[[1]][1]+vc[[2]][1])-1
	agg.patch[i] <- exp(vc[[1]][1])-1
	}
#AIC <- -2*loglik+2*(k+1) #+1 for random effects
AICc <- AIC+2*(k+1)*((k+1)+1)/(length(residuals(mod.full))-(k+1)-1)
dAICc <- AICc-min(AICc)
rL <- exp(-0.5*dAICc)
w <- rL/sum(rL)
chdev <- (dev[1]-dev)/dev[1]

#chagg <- (agg[1]-agg)/agg[1]
#chagg.patch <- (agg.patch[1]-agg.patch)/agg.patch[1]

#QAIC <- -2*loglik/c+2*((k+1)+1) #+1 again for overdispersion
#QAICc <- QAIC+2*((k+1)+1)*(((k+1)+1)+1)/(length(residuals(mod.full))-((k+1)+1)-1)
#dQAICc <- QAICc-min(QAICc)
#QrL <- exp(-0.5*dQAICc)
#Qw <- QrL/sum(QrL)

bestmod <- data.frame(model.list,k,AIC,AICc,dAICc,rL,w,chdev)
write.csv(bestmod, "bestmod.csv", row.names=FALSE, col.names=TRUE, quote=FALSE)

##model averaging

library(AICcmodavg)

new.litter <- seq(min(palm$litter),max(palm$litter),length.out=100)
new.ppext.stem <- seq(min(palm$ppext.stem),max(palm$ppext.stem),length.out=100)

bestmod.AICc <- bestmod$model[bestmod$dAICc<2]
bestmod.AICc.mods <- vector("list", length=length(bestmod.AICc))
for (i in 1:length(bestmod.AICc)){
	bestmod.AICc.mods[[i]]<-lmer(as.character(bestmod.AICc[i]),data=palm,family=poisson)
	}

pred.litter <- modavgpred(bestmod.AICc.mods, modnames=bestmod.AICc, newdata=data.frame(litter=new.litter,ppext.stem=mean(palm$ppext.stem)), type="response", uncond.se="revised")
pred.ppext.stem <- modavgpred(bestmod.AICc.mods, modnames=bestmod.AICc, newdata=data.frame(ppext.stem=new.ppext.stem,litter=mean(palm$litter)), type="response", uncond.se="revised")

par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(2.5,1,0))
plot(seedling~litter, data=palm, ylab="Seedling count", xlab="Leaf litter depth /cm")
lines(pred.litter$mod.avg.pred~new.litter)
lines(pred.litter$mod.avg.pred+1.96*pred.litter$uncond.se~new.litter, lty=2)
lines(pred.litter$mod.avg.pred-1.96*pred.litter$uncond.se~new.litter, lty=2)
text(0,175,labels="(a)", xpd=TRUE)
plot(seedling~ppext.stem, data=palm, ylab="Seedling count", xlab="No. of cultivated stems")
lines(pred.ppext.stem$mod.avg.pred~new.ppext.stem)
lines(pred.ppext.stem$mod.avg.pred+1.96*pred.ppext.stem$uncond.se~new.ppext.stem, lty=2)
lines(pred.ppext.stem$mod.avg.pred-1.96*pred.ppext.stem$uncond.se~new.ppext.stem, lty=2)
text(-40,175,labels="(b)", xpd=TRUE)






##trash code
macarthur <- read.table("macarthur37.txt", header=TRUE, sep="\t")
macarthur <- macarthur[-32,]

mean <- tapply(macarthur$seedling, macarthur$position, mean)
sd <- tapply(macarthur$seedling, macarthur$position, sd)/sqrt(tapply(macarthur$seedling, macarthur$position, length))
bar <- barplot(mean, ylim=c(0,110))
arrows(bar, mean+sd, bar, mean-sd, code=3, angle=90)

par(mfrow=c(3,3), mar=c(4,4,1,1))
plot(seedling~adult.cluster, data=macarthur)
plot(seedling~adult.stem, data=macarthur)
plot(seedling~species.richness, data=macarthur)
plot(seedling~canopy, data=macarthur)
plot(seedling~litter, data=macarthur)
plot(seedling~area, data=macarthur)
plot(seedling~perimeter, data=macarthur)
plot(seedling~shape, data=macarthur)
plot(seedling~width, data=macarthur)

par(mfrow=c(3,3))
plot(seedling+sapling~adult.cluster, data=macarthur)
plot(seedling+sapling~adult.stem, data=macarthur)
plot(seedling+sapling~species.richness, data=macarthur)
plot(seedling+sapling~canopy, data=macarthur)
plot(seedling+sapling~litter, data=macarthur)
plot(seedling+sapling~area, data=macarthur)
plot(seedling+sapling~perimeter, data=macarthur)
plot(seedling+sapling~shape, data=macarthur)
plot(seedling+sapling~width, data=macarthur)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0,1,0,1))
	r.test <- cor.test(x,y)
	txt <- format(c(r.test$estimate, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, font=if(r.test$p.value<0.05) {2} else {1}, cex=cex.cor*abs(r.test$estimate))
	}
panel.hist <- function(x, ...) {
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(usr[1:2], 0, 1.5) )
	h <- hist(x, plot = FALSE)
	breaks <- h$breaks; nB <- length(breaks)
	y <- h$counts; y <- y/max(y)
	rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
	}
pairs(macarthur[3:10], upper.panel=panel.cor, lower.panel=panel.smooth, diag.panel=panel.hist)

models <- read.table("modelset.txt", header=TRUE, sep="\t")

models$formula <- rep("NA", length(models[,1]))
for (i in 1:length(models[,1])) {
	models$formula[i] <- paste(names(models[-1:-3])[models[i,-1:-3]==1], collapse="+")
	}
models$formula <- paste(models$response, models$formula, sep="~")
models$formula[1] <- "seedling~NULL"

models$deviance.po <- rep(NaN, length(models[,1]))
models$aic.po <- rep(NaN, length(models[,1]))
models$dispersion.po <- rep(NaN, length(models[,1]))
for (i in 1:length(models[,1])) {
	model.po <- glm(models$formula[i], data=macarthur, family=poisson)
	models$deviance.po[i] <- (model.po$null-model.po$deviance)/model.po$deviance
	models$aic.po[i] <- model.po$aic
	models$dispersion.po[i] <- model.po$deviance/model.po$df.residual
	}

library(MASS)

models$deviance.nb <- rep(NaN, length(models[,1]))
models$aic.nb <- rep(NaN, length(models[,1]))
models$dispersion.nb <- rep(NaN, length(models[,1]))
models$k <- rep(NaN, length(models[,1]))
for (i in 1:length(models[,1])) {
	model.nb <- glm.nb(models$formula[i], data=macarthur, link="log")
	models$deviance.nb[i] <- (model.nb$null-model.nb$deviance)/model.nb$deviance
	models$aic.nb[i] <- model.nb$aic
	models$dispersion.nb[i] <- model.nb$deviance/model.nb$df.residual
	models$k[i] <- model.nb$df.null-model.nb$df.residual
	}
models$aicc.nb <- models$aic.nb+2*models$k*(models$k+1)/(length(macarthur[,1])-models$k-1)
models$daicc.nb <- models$aicc.nb-min(models$aicc.nb)
models$rLc <- exp(-0.5*models$daicc.nb)
models$waicc <- models$rLc/sum(models$rLc)

write.table(models, "modelsout.txt", sep="\t", quote=FALSE, row.names=FALSE)

model.litter <- glm.nb(seedling~litter, data=macarthur, link="log")
predict <- predict.glm(model.litter, newdata=data.frame(litter=seq(min(macarthur$litter),max(macarthur$litter),0.01)), type="response", se.fit=TRUE)
type <- rep(NA,length(macarthur$adult.stem))
for (i in 1:length(macarthur$adult.stem)) {
	type[i]<-if(macarthur$adult.stem[i]==0) {1} else {48+as.numeric(macarthur$adult.stem[i])}
	}
par(mar=c(4,4,1,1), mgp=c(2.5,1,0))
plot(macarthur$seedling~macarthur$litter, pch=type, xlab="Average leaf litter depth /cm", ylab="No. of Macarthur's Palm seedlings")
lines(seq(min(macarthur$litter),max(macarthur$litter),0.01),predict$fit, lty=1, col="red")
lines(seq(min(macarthur$litter),max(macarthur$litter),0.01),predict$fit+1.96*predict$se.fit, lty=2, col="blue")
lines(seq(min(macarthur$litter),max(macarthur$litter),0.01),predict$fit-1.96*predict$se.fit, lty=2, col="blue")
legend("topright", c("No. of adult stems within transect","No adult stems within transect", "Fitted negative-binomial regression", "95% Confidence interval"), pch=c(35,1,-1,-1), lty=c(-1,-1,1,2), col=c("black", "black", "red", "blue"), bty="n")

######################################

mod.final <- lm(log(seedling+0.1)~leaf+intpp+aspect+extpp+community+extpp:leaf+aspect:community, data=palm)

median(palm$community)->median.community
median(palm$extpp)->median.extpp
median(palm$leaf)->median.leaf

community<-seq(min(palm$community), max(palm$community), 0.1)
extpp <- seq(min(palm$extpp), max(palm$extpp), 1)
leaf<-seq(min(palm$leaf), max(palm$leaf), 0.05)

community.int.pred<-predict(mod.final,newdata=data.frame(community, extpp=rep(median.extpp,length(community)),leaf=rep(median.leaf,length(community)),intpp=rep("Present", length(community)),aspect=rep("Interior",length(community))))
community.edg.pred<-predict(mod.final,newdata=data.frame(community, extpp=rep(median.extpp,length(community)),leaf=rep(median.leaf,length(community)),intpp=rep("Present", length(community)),aspect=rep("Edge",length(community))))
plot(exp(community.edg.pred)~community, type="n", ylab="Predicted no. of seedlings", xlab="Species richness")
lines(exp(community.int.pred)~community, col="blue", lwd=2)
lines(exp(community.edg.pred)~community, col="red", lwd=2)

library(rgl) #you're going to need to install the package rgl

nfacets <- 20+1
grid.x1<-seq(0,210,length=nfacets)
grid.x2<-seq(0.9,8.0,length=nfacets)

fint <- function(a,b){
	predict(mod.final,newdata=data.frame(extpp=a,leaf=b, intpp=rep("Present", length(a)), community=rep(median.community,length(a)), aspect=rep("Interior",length(a))))
	}
fedg <- function(a,b){
	predict(mod.final,newdata=data.frame(extpp=a,leaf=b, intpp=rep("Present", length(a)), community=rep(median.community,length(a)), aspect=rep("Edge",length(a))))
	}

int.pred3d <- outer(grid.x1, grid.x2, fint)
edg.pred3d <- outer(grid.x1, grid.x2, fedg)

open3d()
plot3d(grid.x1, grid.x2, exp(int.pred3d), type="n", xlab="No. of external cultivated stems", ylab="Leaf litter depth /cm", zlab="Predicted no. of seedlings")
surface3d(grid.x1, grid.x2, exp(int.pred3d), col="blue", back="line", front="line", lwd=2, ambient="blue")
surface3d(grid.x1, grid.x2, exp(edg.pred3d), col="red", back="line", front="line", lwd=2, ambient="red")

#the 3d plot is probably sufficient to illustrate the interaction between
#external propagule pressure and leaf litter depth
#but here's some additional 2d plots to complement the visualization
#if you feel the need

extpp.edg.pred <- predict(mod.final,newdata=data.frame(community=rep(median.community,length(extpp)), extpp,leaf=rep(median.leaf,length(extpp)),intpp=rep("Present", length(extpp)),aspect=rep("Edge",length(extpp))))
extpp.int.pred <- predict(mod.final,newdata=data.frame(community=rep(median.community,length(extpp)), extpp,leaf=rep(median.leaf,length(extpp)),intpp=rep("Present", length(extpp)),aspect=rep("Interior",length(extpp))))
plot(exp(extpp.edg.pred)~extpp, type="n", ylim=c(4.5,8.5))
lines(exp(extpp.int.pred)~extpp, col="blue", lwd=2)
lines(exp(extpp.edg.pred)~extpp, col="red", lwd=2)

leaf.edg.pred <- predict(mod.final,newdata=data.frame(community=rep(median.community,length(leaf)), extpp=rep(median.extpp,length(leaf)),leaf,intpp=rep("Present", length(leaf)),aspect=rep("Edge",length(leaf))))
leaf.int.pred <- predict(mod.final,newdata=data.frame(community=rep(median.community,length(leaf)), extpp=rep(median.extpp,length(leaf)),leaf,intpp=rep("Present", length(leaf)),aspect=rep("Interior",length(leaf))))
plot(exp(leaf.int.pred)~leaf, type="n")
lines(exp(leaf.int.pred)~leaf, col="blue", lwd=2)
lines(exp(leaf.edg.pred)~leaf, col="red", lwd=2)

comm.mds <- metaMDS(comm.subsubset, k=3)
comm.envfit <- data.frame(cor(comm.mds$point[,1:3],cbind(palm.subsubset[,17:18], palm.subsubset[,29:31]), method="spearman"))

par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(2,1,0))

plot(comm.mds$species[,1], comm.mds$species[,2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", ylim=c(-1,0.8),xlim=c(-1.3,1.5),type="n")
text(comm.mds$species[,1], comm.mds$species[,2], label=rownames(comm.mds$species), col="red", cex=0.7)
arrows(0,0,as.numeric(comm.envfit[1,])*2,as.numeric(comm.envfit[2,])*2,col="blue", length=0.1)
text(comm.envfit[1,]*2.2,comm.envfit[2,]*2.2, labels=names(comm.envfit),col="blue", cex=0.8)
text(-1.7,0.8,labels="(a)",xpd=TRUE)

plot(comm.mds$species[,3], comm.mds$species[,2], xlab="NMDS Axis 3", ylab="NMDS Axis 2", ylim=c(-1,0.8),xlim=c(-1,1.5), type="n")
text(comm.mds$species[,3], comm.mds$species[,2], label=rownames(comm.mds$species), col="red", cex=0.7)
text(1,-0.8,paste("stress = ",round(comm.mds$stress*100,1)))
arrows(0,0,as.numeric(comm.envfit[3,])*2,as.numeric(comm.envfit[2,])*2,col="blue", length=0.1)
text(comm.envfit[3,]*2.2,comm.envfit[2,]*2.2, labels=names(comm.envfit),col="blue", cex=0.8)
text(-1.4,0.8,labels="(b)",xpd=TRUE)

ordisurf(comm.mds, palm$litter, col="green", add=TRUE)

occ <- specnumber(t(comm))
text(comm.mds$species[occ>10,1],comm.mds$species[occ>10,2], labels=rownames(comm.mds$species[occ>10]), col="blue", cex=0.5)