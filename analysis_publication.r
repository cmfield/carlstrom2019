# This R package and dependencies will need to be installed:
library(devtools)
install_github("cmfield/phylloR")

library(phylloR)

# Import and format all of the data
expt53raw <- read.delim("expt53/otutab_final_classified.txt",stringsAsFactors=F,header=T,row.names=1,comment.char="")
expt61raw <- read.delim("expt61/otutab_final_classified.txt",stringsAsFactors=F,header=T,row.names=1,comment.char="")
expt62raw <- read.delim("expt62/otutab_final_classified.txt",stringsAsFactors=F,header=T,row.names=1,comment.char="")
colnames(expt53raw) <- sub("X","E53S",colnames(expt53raw))
colnames(expt61raw) <- sub("X","E61S",colnames(expt61raw))
colnames(expt62raw) <- sub("X","E62S",colnames(expt62raw))

meta53raw <- read.csv("metadata53.csv",sep=";",header=T,row.names=1)
meta61raw <- read.csv("metadata61.csv",sep=";",header=T,row.names=1)
meta62raw <- read.csv("metadata62.csv",sep=";",header=T,row.names=1)
rownames(meta53raw) <- paste("E53S",rownames(meta53raw),sep="")
rownames(meta61raw) <- paste("E61S",rownames(meta61raw),sep="")
rownames(meta62raw) <- paste("E62S",rownames(meta62raw),sep="")

# The data frame will be row-wise, ie: each row is an experiment, each column a phylogroup
# Remove the inoculum and axenic data, combine unclassified otus
expt53 <- expt53raw[,which(meta53raw$Treatment!="Ax" &  meta53raw$Spray!="Inoc")]
expt53 <- expt53[order(rownames(expt53)),]
expt53 <- rbind(expt53[!grepl("Unclass",rownames(expt53)),],Unclassified=apply(expt53[grepl("Unclass",rownames(expt53)),],2,sum))
meta53 <- meta53raw[which(meta53raw$Treatment!="Ax" &  meta53raw$Spray!="Inoc"),]
meta53 <- cbind(experiment=53,meta53[,1:4])
colnames(meta53) <- c("experiment","initial","spray","repeat","time")
meta53$initial <- factor(meta53$initial,c("ALL","Ax","No A","No B","No G","No P"))
meta53$spray <- factor(meta53$spray,c("U","Mg","A","B","G","P"))

expt61 <- expt61raw[,which(meta61raw$Treatment!="Axenic" & meta61raw$Time!="t0")]
expt61 <- expt61[order(rownames(expt61)),]
expt61 <- rbind(expt61[!grepl("Unclass",rownames(expt61)),],Unclassified=apply(expt61[grepl("Unclass",rownames(expt61)),],2,sum))
expt61 <- expt61[rownames(expt61)!="Leaf281",]
meta61 <- meta61raw[which(meta61raw$Treatment!="Axenic" & meta61raw$Time!="t0"),]
meta61 <- cbind(experiment=61,meta61[,1],"U",meta61[,2:3])
colnames(meta61) <- c("experiment","initial","spray","repeat","time")
levels(meta61$initial) <- sub("-","Leaf",levels(meta61$initial))
meta61$initial <- factor(meta61$initial,c("ALL","Ax",levels(meta61$initial)[!levels(meta61$initial)%in%c("ALL","Ax")]))
meta61$initial <- factor(meta61$initial,c("ALL","Ax",levels(meta61$initial)[grepl("Leaf",levels(meta61$initial))][order(match(levels(meta61$initial)[grepl("Leaf",levels(meta61$initial))],rownames(leafTaxonomy)))]))
meta61$initial <- droplevels(meta61$initial)

expt62 <- expt62raw[,which(meta62raw$Treatment!="Ax" & meta62raw$Time!="t0")]
expt62 <- expt62[order(rownames(expt62)),]
expt62 <- rbind(expt62[!grepl("Unclass",rownames(expt62)),],Unclassified=apply(expt62[grepl("Unclass",rownames(expt62)),],2,sum))
meta62 <- meta62raw[which(meta62raw$Treatment!="Ax" & meta62raw$Time!="t0"),]
meta62 <- cbind(experiment=62,meta62[,1],"U",meta62[,2:3])
colnames(meta62) <- c("experiment","initial","spray","repeat","time")
levels(meta62$initial) <- sub("-","Leaf",levels(meta62$initial))
meta62$initial <- factor(meta62$initial,c("ALL","Ax",levels(meta62$initial)[!levels(meta62$initial)%in%c("ALL","Ax")]))
meta62$initial <- factor(meta62$initial,c("ALL","Ax",levels(meta62$initial)[grepl("Leaf",levels(meta62$initial))][order(match(levels(meta62$initial)[grepl("Leaf",levels(meta62$initial))],rownames(leafTaxonomy)))]))
meta62$initial <- droplevels(meta61$initial)

# Make datasets compatible with analysis functions
#ds51 <- list(counts=expt51,meta=meta51)
ds53 <- list(counts=expt53,meta=meta53)
ds61 <- list(counts=expt61,meta=meta61)
ds62 <- list(counts=expt62,meta=meta62)
ds6c <- list(counts=cbind(expt61,expt62),meta=rbind(meta61,meta62))

# Make control-only datasets and export
control53 <- expt53[,meta53$initial=="ALL"]
control61 <- expt61[,meta61$initial=="ALL"]
control62 <- expt62[,meta62$initial=="ALL"]

write.table(control53,"network/control53.tsv",sep="\t")
write.table(control61,"network/control61.tsv",sep="\t")
write.table(control62,"network/control62.tsv",sep="\t")

cat("Data imported\n")

zcols = c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026")

# Make handy strain sets
strains <- rownames(ds53$counts)
strains <- strains[-length(strains)]
subTax <- leafTaxonomy[strains,]
sNames <- subTax$Name
alpha <- rownames(subTax)[subTax$Class=="Alphaproteobacteria"]
aNames <- subTax[subTax$Class=="Alphaproteobacteria",]$Name
beta <- rownames(subTax)[subTax$Class=="Betaproteobacteria"]
bNames <- subTax[subTax$Class=="Betaproteobacteria",]$Name
gamma <- rownames(subTax)[subTax$Class=="Gammaproteobacteria"]
gNames <- subTax[subTax$Class=="Gammaproteobacteria",]$Name
proteo <- rownames(subTax)[subTax$Phylum=="Proteobacteria"]
pNames <- subTax[subTax$Phylum=="Proteobacteria",]$Name
notalpha <- rownames(subTax)[subTax$Class!="Alphaproteobacteria"]
naNames <- subTax[subTax$Class!="Alphaproteobacteria",]$Name
notbeta <- rownames(subTax)[subTax$Class!="Betaproteobacteria"]
nbNames <- subTax[subTax$Class!="Betaproteobacteria",]$Name
notgamma <- rownames(subTax)[subTax$Class!="Gammaproteobacteria"]
ngNames <- subTax[subTax$Class!="Gammaproteobacteria",]$Name
notproteo <- rownames(subTax)[subTax$Phylum!="Proteobacteria"]
npNames <- subTax[subTax$Phylum!="Proteobacteria",]$Name

if(!exists("nocalc")){
cds53 <- list(conspray=makeCDS(ds53,include=list(initial="ALL"),foi="spray",title="Effect of Mock Spray\n(all strains)",legend=c("Unsprayed","Mock Sprayed")),
    contime=makeCDS(ds53,include=list(initial="ALL"),foi="time",title="Effect of Time\n(all strains)",legend=c("Timepoint 1","Timepoint 2")),
    
    noalpha=makeCDS(ds53,include=list(initial=c("No A","ALL"),spray=c("U","Mg")),foi="initial",title="",legend=c("Control","Absence")),
    alphaback=makeCDS(ds53,include=list(initial="No A",time="t2"),foi="spray",title="",legend=c("Mock Spray","Reintroduction")),
    latealpha=makeCDS(ds53,include=list(initial=c("No A","ALL"),spray=c("U","Mg","A"),time="t2"),exclude=list(initial="No A",spray="Mg"),foi="initial",title="",legend=c("Control","Late Arrival")),
    
    nobeta=makeCDS(ds53,include=list(initial=c("No B","ALL"),spray=c("U","Mg")),foi="initial",title="",legend=c("Control","Absence")),
    betaback=makeCDS(ds53,include=list(initial="No B",time="t2"),foi="spray",title="",legend=c("Mock Spray","Reintroduction")),
    latebeta=makeCDS(ds53,include=list(initial=c("No B","ALL"),spray=c("U","Mg","B"),time="t2"),exclude=list(initial="No B",spray="Mg"),foi="initial",title="",legend=c("Control","Late Arrival")),
    
    nogamma=makeCDS(ds53,include=list(initial=c("No G","ALL"),spray=c("U","Mg")),foi="initial",title="",legend=c("Control","Absence")),
    gammaback=makeCDS(ds53,include=list(initial="No G",time="t2"),foi="spray",title="",legend=c("Mock Spray","Reintroduction")),
    lategamma=makeCDS(ds53,include=list(initial=c("No G","ALL"),spray=c("U","Mg","G"),time="t2"),exclude=list(initial="No G",spray="Mg"),foi="initial",title="",legend=c("Control","Late Arrival")),

    noproteo=makeCDS(ds53,include=list(initial=c("No P","ALL"),spray=c("U","Mg")),foi="initial",title="",legend=c("Control","Absence")),
    proteoback=makeCDS(ds53,include=list(initial="No P",time="t2"),foi="spray",title="",legend=c("Mock Spray","Reintroduction")),
    lateproteo=makeCDS(ds53,include=list(initial=c("No P","ALL"),spray=c("U","Mg","P"),time="t2"),exclude=list(initial="No P",spray="Mg"),foi="initial",title="",legend=c("Control","Late Arrival"))
)

cds6 <- list()
cds6[["batch"]] <- makeCDS(ds6c,include=list(initial="ALL"),foi="experiment",title="R1 vs. R2",legend=c("R1","R2"))
ds633 <- ds6c
ds633$counts <- ds633$counts[,(ds633$meta$experiment=="62" & ds633$meta$initial=="Leaf33") | (ds633$meta$experiment=="61" & ds633$meta$initial=="ALL")]
ds633$meta <- ds633$meta[(ds633$meta$experiment=="62" & ds633$meta$initial=="Leaf33") | (ds633$meta$experiment=="61" & ds633$meta$initial=="ALL"),]
cds6[["dropout"]] <- makeCDS(ds633,foi="experiment",title="R1 Control vs.\nR2 L-33 Sphingomonas Drop-out",legend=c("R1 Control","R2 Drop-out"))
for(x in levels(ds6c$meta$initial)[-1]){
	cds6[[paste(x,1,sep="")]] <- makeCDS(ds61,include=list(initial=c(x,"ALL")),foi="initial",title="")
	cds6[[paste(x,2,sep="")]] <- makeCDS(ds62,include=list(initial=c(x,"ALL")),foi="initial",title="")
    cds6[[paste(x,"c",sep="")]] <- makeCDS(ds6c,include=list(initial=c(x,"ALL")),foi="initial",ftc="experiment",title="")
}
}

pm = 10000
cf = 0.01

bicols <- c("#214478","#A1D662")

adn53 <- list()

cairo_pdf("figures/class_dropout_control_spraytime.pdf",width=16,height=7,family="Arial")
par(mfrow=c(1,2),cex=1.5,cex.lab=1,cex.sub=1.2)
adn53$conspray <- plotPCA(cds53$conspray,soi=strains,perm=pm,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,showArrows=F)
adn53$contime <- plotPCA(cds53$contime,soi=strains,perm=pm,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,showArrows=F)
dev.off()

cairo_pdf("figures/class_dropout_pca.pdf",width=28,height=35,family="Arial")
par(mfrow=c(5,4),cex=1.5,oma=c(0,6,6,0),cex.lab=1,cex.sub=1.2)
adn53$noalpha <- plotPCA(cds53$noalpha,soi=notalpha,perm=pm,cutoff=cf,rowLabs=naNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$nobeta <- plotPCA(cds53$nobeta,soi=notbeta,perm=pm,cutoff=cf,rowLabs=nbNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$nogamma <- plotPCA(cds53$nogamma,soi=notgamma,perm=pm,cutoff=cf,rowLabs=ngNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$noproteo <- plotPCA(cds53$noproteo,soi=notproteo,perm=pm,cutoff=cf,rowLabs=npNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)

adn53$alphaback <- plotPCA(cds53$alphaback,soi=notalpha,perm=pm,cutoff=cf,rowLabs=naNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$betaback <- plotPCA(cds53$betaback,soi=notbeta,perm=pm,cutoff=cf,rowLabs=nbNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$gammaback <- plotPCA(cds53$gammaback,soi=notgamma,perm=pm,cutoff=cf,rowLabs=ngNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$proteoback <- plotPCA(cds53$proteoback,soi=notproteo,perm=pm,cutoff=cf,rowLabs=npNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)

adn53$latealphanot <- plotPCA(cds53$latealpha,soi=notalpha,perm=pm,cutoff=cf,rowLabs=naNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$latebetanot <- plotPCA(cds53$latebeta,soi=notbeta,perm=pm,cutoff=cf,rowLabs=nbNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$lategammanot <- plotPCA(cds53$lategamma,soi=notgamma,perm=pm,cutoff=cf,rowLabs=ngNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$lateproteonot <- plotPCA(cds53$lateproteo,soi=notproteo,perm=pm,cutoff=cf,rowLabs=npNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)

adn53$latealpha <- plotPCA(cds53$latealpha,soi=alpha,perm=pm,cutoff=cf,rowLabs=aNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$latebeta <- plotPCA(cds53$latebeta,soi=beta,perm=pm,cutoff=cf,rowLabs=bNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$lategamma <- plotPCA(cds53$lategamma,soi=gamma,perm=pm,cutoff=cf,rowLabs=gNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$lateproteo <- plotPCA(cds53$lateproteo,soi=proteo,perm=pm,cutoff=cf,rowLabs=pNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)

adn53$latealphaall <- plotPCA(cds53$latealpha,soi=strains,perm=pm,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$latebetaall <- plotPCA(cds53$latebeta,soi=strains,perm=pm,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$lategammaall <- plotPCA(cds53$lategamma,soi=strains,perm=pm,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
adn53$lateproteoall <- plotPCA(cds53$lateproteo,soi=strains,perm=pm,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)

mtext("Drop-Out Condition",side=3,line=3,at=0.515,outer=T,cex=3)
mtext("Alphaproteobacteria",side=3,line=0,at=0.14,outer=T,cex=3)
mtext("Betaproteobacteria",side=3,line=0,at=0.39,outer=T,cex=3)
mtext("Gammaproteobacteria",side=3,line=0,at=0.64,outer=T,cex=3)
mtext("Proteobacteria",side=3,line=0,at=0.89,outer=T,cex=3)

mtext("Group Absent vs. Control (a)",side=2,line=1,at=0.93,outer=T,cex=3)
mtext("Late Arrival vs. Mock (b)",side=2,line=1,at=0.73,outer=T,cex=3)
mtext("Late Arrival vs. Control (c)",side=2,line=1,at=0.53,outer=T,cex=3)
mtext("Late Arrival vs. Control (c)",side=2,line=1,at=0.33,outer=T,cex=3)
mtext("Late Arrival vs. Control (c)",side=2,line=1,at=0.13,outer=T,cex=3)
mtext("Effect on Rest of the Community",side=2,line=4,at=0.73,outer=T,cex=3)
mtext("Effect on Invading Group",side=2,line=4,at=0.33,outer=T,cex=3)
mtext("Effect on Whole Community",side=2,line=4,at=0.13,outer=T,cex=3)

mtext(bquote(underline("                                                                                                                                        ")),side=2,line=3,at=0.73,outer=T,cex=3)

dev.off()

cairo_pdf("figures/class_dropout_control_spraytime_bar.pdf",width=16,height=7,family="Arial")
par(mfrow=c(1,2),cex=1)
plotCommunityChanges(cds53$contime,soi=strains,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,nBars=18)
plotCommunityChanges(cds53$conspray,soi=strains,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,nBars=18)
dev.off()

cairo_pdf("figures/class_dropout_bar.pdf",width=28,height=21,family="Arial")
par(mfrow=c(3,4),cex=1,oma=c(0,6,4,0))
plotCommunityChanges(cds53$noalpha,soi=notalpha,cutoff=cf,rowLabs=naNames,subtitle="",cols=bicols,nBars=18)
plotCommunityChanges(cds53$nobeta,soi=notbeta,cutoff=cf,rowLabs=nbNames,subtitle="",cols=bicols,nBars=18)
plotCommunityChanges(cds53$nogamma,soi=notgamma,cutoff=cf,rowLabs=ngNames,subtitle="",cols=bicols,nBars=18)
plotCommunityChanges(cds53$noproteo,soi=notproteo,cutoff=cf,rowLabs=npNames,subtitle="",cols=bicols,nBars=18)

plotCommunityChanges(cds53$alphaback,soi=notalpha,cutoff=cf,rowLabs=naNames,subtitle="",cols=bicols,nBars=18)
plotCommunityChanges(cds53$betaback,soi=notbeta,cutoff=cf,rowLabs=nbNames,subtitle="",cols=bicols,nBars=18)
plotCommunityChanges(cds53$gammaback,soi=notgamma,cutoff=cf,rowLabs=ngNames,subtitle="",cols=bicols,nBars=18)
plotCommunityChanges(cds53$proteoback,soi=notproteo,cutoff=cf,rowLabs=npNames,subtitle="",cols=bicols,nBars=18)

plotCommunityChanges(cds53$latealpha,soi=strains,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,nBars=18)
plotCommunityChanges(cds53$latebeta,soi=strains,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,nBars=18)
plotCommunityChanges(cds53$lategamma,soi=strains,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,nBars=18)
plotCommunityChanges(cds53$lateproteo,soi=strains,cutoff=cf,rowLabs=sNames,subtitle="",cols=bicols,nBars=18)

mtext("Alphaproteobacteria",side=3,line=0,at=0.14,outer=T,cex=3)
mtext("Betaproteobacteria",side=3,line=0,at=0.39,outer=T,cex=3)
mtext("Gammaproteobacteria",side=3,line=0,at=0.64,outer=T,cex=3)
mtext("Proteobacteria",side=3,line=0,at=0.89,outer=T,cex=3)

mtext("Group Absent vs. Control (a)",side=2,line=1,at=0.853,outer=T,cex=3)
mtext("Late Arrival vs. Mock (b)",side=2,line=1,at=0.52,outer=T,cex=3)
mtext("Late Arrival vs. Control (c)",side=2,line=1,at=0.187,outer=T,cex=3)
mtext("Effect on Rest of the Community",side=2,line=4,at=0.853,outer=T,cex=3)
mtext("Effect on Invading Group",side=2,line=4,at=0.52,outer=T,cex=3)
mtext("Effect on Whole Community",side=2,line=4,at=0.187,outer=T,cex=3)

dev.off()

adn6 <- list()

cairo_pdf("figures/single_dropout_batch_pca.pdf",width=14,height=7,family="Arial")
par(mfrow=c(1,2),cex=1.5,cex.lab=1,cex.sub=1.2)
adn6[["batch"]] <- plotPCA(cds6[["batch"]],soi=strains,perm=pm,cutoff=cf,rowLabs=sNames,cols=bicols,showArrows=F)
adn6[["dropout"]] <- plotPCA(cds6[["dropout"]],soi=strains,perm=pm,cutoff=cf,rowLabs=sNames,cols=bicols,showArrows=F)
dev.off()

cairo_pdf("figures/single_dropout_pca.pdf",width=21,height=175,family="Arial")
par(mfrow=c(25,3),cex=1.5,oma=c(0,8,6,0),cex.lab=1,cex.sub=1.2)
for(x in levels(ds6c$meta$initial)[-1]){
    adn6[[paste(x,1,sep="")]] <- plotPCA(cds6[[paste(x,1,sep="")]],soi=strains[strains!=x],perm=pm,cutoff=cf,rowLabs=sNames[strains!=x],subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
    adn6[[paste(x,2,sep="")]] <- plotPCA(cds6[[paste(x,2,sep="")]],soi=strains[strains!=x],perm=pm,cutoff=cf,rowLabs=sNames[strains!=x],subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
    adn6[[paste(x,"c",sep="")]] <- plotPCA(cds6[[paste(x,"c",sep="")]],soi=strains[strains!=x],perm=pm,cutoff=cf,rowLabs=sNames[strains!=x],subtitle="",cols=bicols,showLegend=F,showArrows=F,showTitle=F)
    mtext(leafTaxonomy[x,]$Name,side=2,line=1,at=1.065-which(levels(ds6c$meta$initial[-1])==x)/25,outer=T,cex=3)
}
mtext("Replicate 1",side=3,line=0,at=0.187,outer=T,cex=3)
mtext("Replicate 2",side=3,line=0,at=0.52,outer=T,cex=3)
mtext("Combined",side=3,line=0,at=0.853,outer=T,cex=3)

mtext("Drop-Out Condition",side=2,line=5,at=0.5,outer=T,cex=4)
mtext("Experiment",side=3,line=4,at=0.52,outer=T,cex=4)
dev.off()

cairo_pdf("figures/single_dropout_highlights_pca.pdf",width=28,height=7,family="Arial")
par(mfrow=c(1,4),cex=1,cex.main=1.5,cex.sub=1.5)
for(x in c("Leaf203","Leaf231","Leaf233","Leaf262")){
    adn6[[paste(x,"c",sep="")]] <- plotPCA(cds6[[paste(x,"c",sep="")]],soi=strains[strains!=x],perm=pm,cutoff=cf,rowLabs=sNames[strains!=x],cols=bicols,showLegend=F,showArrows=F,subtitle=leafTaxonomy[x,]$Name)
}
dev.off()

cairo_pdf("figures/single_dropout_batch_bar.pdf",width=14,height=7,family="Arial")
par(mfrow=c(1,2),cex=1)
plotCommunityChanges(cds6[["batch"]],soi=strains,cutoff=cf,rowLabs=sNames,cols=bicols,nBars=12)
plotCommunityChanges(cds6[["dropout"]],soi=strains,cutoff=cf,rowLabs=sNames,cols=bicols,nBars=12)
dev.off()

cairo_pdf("figures/single_dropout_bar.pdf",width=42,height=175,family="Arial")
par(mfrow=c(25,3),cex=1,oma=c(0,4,4,0))
for(x in levels(ds6c$meta$initial)[-1]){
    plotCommunityChanges(cds6[[paste(x,1,sep="")]],soi=strains[strains!=sub("no","Leaf",x)],cutoff=cf,rowLabs=sNames[strains!=sub("no","Leaf",x)],subtitle="",cols=bicols,nBars=34)
    plotCommunityChanges(cds6[[paste(x,2,sep="")]],soi=strains[strains!=sub("no","Leaf",x)],cutoff=cf,rowLabs=sNames[strains!=sub("no","Leaf",x)],subtitle="",cols=bicols,nBars=34)
    plotCommunityChanges(cds6[[paste(x,"c",sep="")]],soi=strains[strains!=sub("no","Leaf",x)],cutoff=cf,rowLabs=sNames[strains!=sub("no","Leaf",x)],subtitle="",cols=bicols,nBars=34)
    mtext(leafTaxonomy[x,]$Name,side=2,line=1,at=1.065-which(levels(ds6c$meta$initial[-1])==x)/25,outer=T,cex=3)
}
mtext("Replicate 1",side=3,line=0,at=0.187,outer=T,cex=3)
mtext("Replicate 2",side=3,line=0,at=0.52,outer=T,cex=3)
mtext("Combined",side=3,line=0,at=0.853,outer=T,cex=3)
dev.off()

# Summary plot for expt53
library(corrplot)
library(plotrix)
library(circlize)
effects53 <- unlist(lapply(adn53,function(x) x$stats$aov.tab$R2[1]))
pvalues53 <- unlist(lapply(adn53,function(x) x$stats$aov.tab$Pr[1]))

e <- t(matrix(effects53[c(-1,-2)],4))
p <- t(matrix(pvalues53[c(-1,-2)],4))
rownames(e) <- c("Group Absent\nvs. Control (a)","Late Arrival\nvs. Mock (b)","Late Arrival\nvs. Control (c)","Late Arrival\nvs. Control (c)","Late Arrival\nvs. Control (c)")
colnames(e) <- c("Alphaproteobacteria","Betaproteobacteria","Gammaproteobacteria","Proteobacteria")

cairo_pdf("figures/class_dropout_summary.pdf",width=6,height=6,family="Arial")
par(cex=1,xpd=T)
corrplot(e,is.corr=F,cl.lim=c(0,0.7),cl.ratio=0.4,p.mat=p,insig="label_sig",sig=c(0.0001,0.001,0.01,0.05),tl.col=1,pch.cex=1,tl.srt=45,mar=c(0,10,2,2),col=zcols[6:1])

text(-4,4.1,"Effect on Rest\nof the Community",font=2)
text(-4,2.1,"Effect on\nInvading Group",font=2)
text(-4,1.1,"Effect on\nWhole Community",font=2)
lines(c(-2.1,-2.1),c(2.5,5.5))

text(2.5,8.5,"Drop-Out Condition",cex=1,font=2)
text(6,3,"Effect Size",cex=1,srt=-90)
dev.off()

# Summary plot 2 for expt53
library(ape)
library(apextra)
tree <- read.tree("../leafGenomes/16s_full_align.phy_phyml_tree.txt")
tree <- drop.tip(tree,tree$tip.label[!tree$tip.label%in%strains])
tree <- root(tree,node=69,resolve.root=T)
tree <- as.ultrametric(tree)
strainOrder <- tree$tip.label[tipOrder(tree)]
strainCols <- subTax[,'Color',F]
strainCols <- strainCols[strainOrder,]

fc53 <- do.call(cbind,lapply(cds53[names(cds53)[-1:-2]],function(x) x$results[,"log2FoldChange"]))
fc53[is.na(fc53)] <- 0
pv53 <- do.call(cbind,lapply(cds53[names(cds53)[-1:-2]],function(x) x$results[,"padj"]))
pv53[is.na(pv53)] <- 1

rownames(fc53) <- rownames(cds53[[1]]$results)
colnames(fc53) <- rep(c("(a)","(b)","(c)"),4)
fc53 <- fc53[1:62,]
pv53 <- pv53[1:62,]

mask <- fc53==fc53
mask[alpha,1:2] <- FALSE
mask[beta,4:5] <- FALSE
mask[gamma,7:8] <- FALSE
mask[proteo,10:11] <- FALSE

overlay <- matrix(cut(pv53,c(0,0.0001,0.001,0.01,0.05,1),c("****","***","**","*","")),nrow(pv53))

cairo_pdf("figures/class_dropout_details.pdf",width=10,height=15,family="Arial")
par(cex=1)
treatmap(tree,fc53,mask=mask,overlay=overlay,tip.labels=leafTaxonomy[tree$tip.label,]$Name,tip.colors=strainCols,aspect.ratio=0.2,tip.label.width=8,z.cols=zcols)
dev.off()

# Summary plot for expt6
effects6 <- unlist(lapply(adn6,function(x) x$stats$aov.tab$R2[1]))
pvalues6 <- unlist(lapply(adn6,function(x) x$stats$aov.tab$Pr[1]))

e <- matrix(effects6[-c(1,2)],3)
rownames(e) <- c("Replicate 1","Replicate 2","Combined")
shortNames <- levels(ds6c$meta$initial)[-1]
colnames(e) <- leafTaxonomy[shortNames,]$Name
p <- matrix(pvalues6[-c(1,2)],3)

e <- e[,order(match(shortNames,tree$tip.label))]
p <- p[,order(match(shortNames,tree$tip.label))]
shortNames <- shortNames[order(match(shortNames,tree$tip.label))]

subtree <- keep.tip(tree,shortNames)

cairo_pdf("figures/single_dropout_summary.pdf",width=14,height=7,family="Arial")
par(cex=1,xpd=T)
rownames(e) <- rep("",3)
corrplot(e,is.corr=F,cl.lim=c(0,0.25),cl.ratio=0.1,tl.col=c(leafTaxonomy[shortNames,]$Color,"black"),p.mat=p,insig="label_sig",sig=c(0.0001,0.001,0.01,0.05),pch.cex=1,mar=c(0,8,0,2),col=zcols[6:1])
draw.phylo(1,8.5,25,11,subtree,direction="d")
text(0.25,1,"Combined",pos=2)
text(0.25,2,"Replicate 2",pos=2)
text(0.25,3,"Replicate 1",pos=2)
text(12.5,11.5,"Drop-Out Condition",cex=1.5)
text(28,2,"Effect Size",cex=1,srt=-90)
legend(0.5,3.6,legend=names(phylumColors[-4]),fill=phylumColors[-4],xjust=1,yjust=0,cex=0.8)
dev.off()

# Summary plot 2 for expt6 combined
cds6c <- cds6[-c(1,2)][grepl("c",names(cds6)[-c(1,2)])]
fc6c <- do.call(cbind,lapply(cds6c,function(x) x$results$log2FoldChange))
fc6c[is.na(fc6c)] <- 0
pv6c <- do.call(cbind,lapply(cds6c,function(x) x$results$padj))
pv6c[is.na(pv6c)] <- 1

rownames(fc6c) <- rownames(cds6c[[1]]$results)
colnames(fc6c) <- sub("c","",colnames(fc6c)) 
fc6c <- fc6c[1:62,]
pv6c <- pv6c[1:62,]

mask <- fc6c==fc6c
for(col in colnames(fc6c)){
    mask[col,col] <- FALSE
}

overlay <- matrix(cut(pv6c,c(0,0.0001,0.001,0.01,0.05,1),c("****","***","**","*","")),nrow(pv6c))

hc <- as.phylo(hclust(dist(fc6c)))
hc <- reorderTips(hc)

cairo_pdf("figures/single_dropout_details_phylo.pdf",width=15,height=20,family="Arial")
par(cex=1)
treatmap(tree,fc6c,mask=mask,overlay=overlay,tip.labels=leafTaxonomy[tree$tip.label,]$Name,tip.colors=leafTaxonomy[tree$tip.label,]$Color,tip.label.width=10,mat.label.height=24,mat.labels=leafTaxonomy[colnames(fc6c),]$Name,mat.label.color=leafTaxonomy[colnames(fc6c),]$Color,mat.hclust=T,z.cols=zcols)
dev.off()

cairo_pdf("figures/single_dropout_details_hclust.pdf",width=15,height=20,family="Arial")
par(cex=1)
treatmap(hc,fc6c,mask=mask,overlay=overlay,tip.labels=leafTaxonomy[hc$tip.label,]$Name,tip.colors=leafTaxonomy[hc$tip.label,]$Color,tip.label.width=10,mat.label.height=24,mat.labels=leafTaxonomy[colnames(fc6c),]$Name,mat.label.color=leafTaxonomy[colnames(fc6c),]$Color,mat.hclust=T,z.cols=zcols)
dev.off()

# Summary plot with 6-1 and 6-2 separate
cds61 <- cds6[grepl("1$",names(cds6))]
fc61 <- do.call(cbind,lapply(cds61,function(x) x$results$log2FoldChange))
fc61[is.na(fc61)] <- 0
pv61 <- do.call(cbind,lapply(cds61,function(x) x$results$padj))
pv61[is.na(pv61)] <- 1

shortNames1 <- sub("1$","",sub("no","Leaf",names(cds61)))
rownames(fc61) <- rownames(cds61[[1]]$results)
colnames(fc61) <- paste("R1",leafTaxonomy[shortNames1,]$Name,sep="-")
fc61 <- fc61[1:62,]
pv61 <- pv61[1:62,]

cds62 <- cds6[grepl("2$",names(cds6))]
fc62 <- do.call(cbind,lapply(cds62,function(x) x$results$log2FoldChange))
fc62[is.na(fc62)] <- 0
pv62 <- do.call(cbind,lapply(cds62,function(x) x$results$padj))
pv62[is.na(pv62)] <- 1

shortNames2 <- sub("2$","",sub("no","Leaf",names(cds62)))
rownames(fc62) <- rownames(cds62[[1]]$results)
colnames(fc62) <- paste("R2",leafTaxonomy[shortNames2,]$Name,sep="-")
fc62 <- fc62[1:62,]
pv62 <- pv62[1:62,]

fcc <- cbind(fc61,fc62)
pvc <- cbind(pv61,pv62)
shortNames <- c(shortNames1,shortNames2)

mask <- fcc==fcc
for(i in 1:(length(shortNames))){
    mask[shortNames[i],i] <- FALSE
}

overlay <- matrix(cut(pvc,c(0,0.0001,0.001,0.01,0.05,1),c("****","***","**","*","")),nrow(pvc))

cairo_pdf("figures/single_dropout_separate_details_phylo.pdf",width=20,height=20,family="Arial")
par(cex=1)
treatmap(tree,fcc,mask=mask,overlay=overlay,tip.labels=leafTaxonomy[tree$tip.label,]$Name,tip.colors=leafTaxonomy[tree$tip.label,]$Color,tip.label.width=10,mat.label.height=24,mat.labels=colnames(fcc),mat.label.color=leafTaxonomy[shortNames,]$Color,mat.col.order=order(match(shortNames,tree$tip.label),decreasing=T),z.cols=zcols)
dev.off()

# Community plots
control53 <- cds53$conspray$counts
control53 <- control53[rownames(control53)!="Unclassified",]
control61 <- ds61$counts
control61 <- control61[,ds61$meta$initial=="ALL"]
control61 <- control61[rownames(control61)!="Unclassified",]
control62 <- ds62$counts
control62 <- control62[,ds62$meta$initial=="ALL"]
control62 <- control62[rownames(control62)!="Unclassified",]

cairo_pdf("figures/class_dropout_control_community.pdf",width=20,height=10,family="Arial")
par(cex=1)
comm53 <- plotCommunity(control53,type="violinswarm",xlabels=leafTaxonomy[rownames(control53),]$Name,xcols=leafTaxonomy[rownames(control53),]$Color)
dev.off()

cairo_pdf("figures/single_dropout_control_community.pdf",width=20,height=10,family="Arial")
par(cex=1)
comm61 <- plotCommunity(control61,type="violinswarm",xlabels=leafTaxonomy[rownames(control61),]$Name,xcols=leafTaxonomy[rownames(control61),]$Color)
comm62 <- plotCommunity(control62,type="violinswarm",xlabels=leafTaxonomy[rownames(control62),]$Name,xcols=leafTaxonomy[rownames(control61),]$Color)
dev.off()

# PCA of controls
cds5361 <- makeCDS(counts=cbind(control53,control61),meta=rbind(cds53$conspray$meta,ds61$meta[ds61$meta$initial=="ALL",]),foi="experiment",legend=c("Class Droupout","Single Dropout R1"))
cds5362 <- makeCDS(counts=cbind(control53,control62),meta=rbind(cds53$conspray$meta,ds62$meta[ds62$meta$initial=="ALL",]),foi="experiment",legend=c("Class Droupout","Single Dropout R2"))
cds6162 <- makeCDS(counts=cbind(control61,control62),meta=rbind(ds61$meta[ds62$meta$initial=="ALL",],ds62$meta[ds62$meta$initial=="ALL",]),foi="experiment",legend=c("Single Dropout R1","Single Dropout R2"))
cdsControl <- makeCDS(counts=cbind(control53,control61,control62),meta=rbind(cds53$conspray$meta,ds61$meta[ds61$meta$initial=="ALL",],ds62$meta[ds62$meta$initial=="ALL",]),foi="experiment",legend=c("Class Droupout","Single Dropout R1","Single Dropout R2"))

cairo_pdf("figures/controls_comparison_pca.pdf",width=14,height=14,family="Arial")
par(mfrow=c(2,2),cex=1)
plotPCA(cds5361,cols=bicols)
plotPCA(cds6162,cols=bicols)
plotPCA(cds5362,cols=bicols)
plotPCA(cdsControl,cols=c("red",bicols))
dev.off()

# Make tables for SparCC: all controls, expt61+expt62, noalpha, nobeta, nogamma, noproteo, alpharestore, betarestore, gammarestore, proteorestore
sparTabs <- list(control53=t(ds53$counts[,ds53$meta$initial=="ALL"]),
                 control61=t(ds61$counts[,ds61$meta$initial=="ALL"]),
                 control62=t(ds62$counts[,ds62$meta$initial=="ALL"])
)
for(name in names(sparTabs)){
    write.table(t(sparTabs[[name]]),paste("network/",name,".tsv",sep=""),sep="\t")
}

# Pie chart of control community
cairo_pdf("figures/class_dropout_pie.pdf",width=7,height=7,family="Arial")
par(cex=1)
plotCommunityPie(control53,strainTaxa=names(phylumColors)[match(leafTaxonomy[rownames(control53),]$Color,phylumColors)],cols=phylumColors,taxLabels=names(phylumColors),sort=F)
dev.off()

cairo_pdf("figures/single_dropout_pie.pdf",width=7,height=7,family="Arial")
par(cex=1)
plotCommunityPie(control61,strainTaxa=names(phylumColors)[match(leafTaxonomy[rownames(control61),]$Color,phylumColors)],cols=phylumColors,taxLabels=names(phylumColors),sort=F)
plotCommunityPie(control62,strainTaxa=names(phylumColors)[match(leafTaxonomy[rownames(control62),]$Color,phylumColors)],cols=phylumColors,taxLabels=names(phylumColors),sort=F)
dev.off()

# Network analysis based on experiment 6 knockouts: 61, 62 and 6c
cds6c <- cds6[grepl("c$",names(cds6))]
names(cds6c) <- substr(names(cds6c),1,nchar(names(cds6c))-1)
names(cds6c) <- sub("no","Leaf",names(cds6c))

summary6 <- summariseResults(cds6c)
summary6$fcMatrix <- summary6$fcMatrix[strains,]
summary6$pcMatrix <- summary6$pvMatrix[strains,]
net6 <- igraphFromSummary(summary6$fcMatrix,summary6$pvMatrix,cutoff=0.01)
vertex_attr(net6,"shortName") <- sub("-.*","",vertex_attr(net6,"name"))
vertex_attr(net6,"twoLineName") <- sub("-","\n",vertex_attr(net6,"name"))
write.graph(net6,"results/network_001.gml","graphml")

# Including the inocula for time series comparison
inoc53 <- expt53raw[,which(meta53raw$Treatment=="ALL")]
inoc53 <- inoc53[order(rownames(inoc53)),]
inoc53 <- rbind(inoc53[!grepl("Unclass",rownames(inoc53)),],Unclassified=apply(inoc53[grepl("Unclass",rownames(inoc53)),],2,sum))
inocmeta53 <- meta53raw[which(meta53raw$Treatment=="ALL"),]
inocmeta53 <- cbind(experiment=53,inocmeta53[,1:4])
colnames(inocmeta53) <- c("experiment","initial","spray","repeat","time")
inocmeta53$initial <- factor(inocmeta53$initial,c("ALL"))
inocmeta53$spray <- factor(inocmeta53$spray,c("U","Mg","Inoc"))
# Turn t1 into t2
inocmeta53$time[inocmeta53$time=="t1"] <- "t2"

inoc61 <- expt61raw[,which(meta61raw$Treatment=="ALL")]
inoc61 <- inoc61[order(rownames(inoc61)),]
inoc61 <- rbind(inoc61[!grepl("Unclass",rownames(inoc61)),],Unclassified=apply(inoc61[grepl("Unclass",rownames(inoc61)),],2,sum))
inoc61 <- inoc61[rownames(inoc61)!="Leaf281",]
inocmeta61 <- meta61raw[which(meta61raw$Treatment=="ALL"),]
inocmeta61 <- cbind(experiment=61,inocmeta61[,1],"U",inocmeta61[,2:3])
colnames(inocmeta61) <- c("experiment","initial","spray","repeat","time")
inocmeta61$initial <- factor(inocmeta61$initial,c("ALL"))
inocmeta61$spray <- factor(inocmeta61$spray,c("U"))

inoc62 <- expt62raw[,which(meta62raw$Treatment=="ALL")]
inoc62 <- inoc62[order(rownames(inoc62)),]
inoc62 <- rbind(inoc62[!grepl("Unclass",rownames(inoc62)),],Unclassified=apply(inoc62[grepl("Unclass",rownames(inoc62)),],2,sum))
inocmeta62 <- meta62raw[which(meta62raw$Treatment=="ALL"),]
inocmeta62 <- cbind(experiment=62,inocmeta62[,1],"U",inocmeta62[,2:3])
colnames(inocmeta62) <- c("experiment","initial","spray","repeat","time")
inocmeta62$initial <- factor(inocmeta62$initial,c("ALL"))
inocmeta62$spray <- factor(inocmeta62$spray,c("U"))

ids53 <- list(counts=inoc53,meta=inocmeta53)
ids61 <- list(counts=inoc61,meta=inocmeta61)
ids62 <- list(counts=inoc62,meta=inocmeta62)
ids5361 <- list(counts=cbind(inoc53,inoc61),meta=rbind(inocmeta53,inocmeta61))
ids5362 <- list(counts=cbind(inoc53,inoc62),meta=rbind(inocmeta53,inocmeta62))
ids6162 <- list(counts=cbind(inoc61,inoc62),meta=rbind(inocmeta61,inocmeta62))
ids536162 <- list(counts=cbind(inoc53,inoc61,inoc62),meta=rbind(inocmeta53,inocmeta61,inocmeta62))

dds53 <- DESeqDataSetFromMatrix(ds53$counts[rownames(ds53$counts)!="Unclassified",ds53$meta$initial=="ALL"],ds53$meta[ds53$meta$initial=="ALL",],~1)
ddsi <- DESeqDataSetFromMatrix(ids53$counts[rownames(ids53$counts)!="Unclassified",ids53$meta$time=="t0"],ids53$meta[ids53$meta$time=="t0",],~1)
vsti <- assay(varianceStabilizingTransformation(ddsi))
vst53 <- assay(varianceStabilizingTransformation(dds53))

median53 <- apply(vst53,1,median)
mediani <- apply(vsti,1,median)

cairo_pdf("figures/class_dropout_winners_losers.pdf",width=40,height=10,family="Arial")
par(cex=1)
icds53 <- makeCDS(ids53,include=list(time=c("t0","t2")),foi="time",legend=c("Inoculum","Established Community"))
plotCommunityChanges(icds53,cutoff=cf,rowLabs=leafTaxonomy[rownames(ids53$counts),]$Name,subtitle="",cols=bicols,nBars=54)
dev.off()

cairo_pdf("figures/class_dropout_inocula.pdf",width=20,height=10,family="Arial")
par(cex=1)
plotCommunity(ids53$counts[1:62,ids53$meta$spray=="Inoc"],xlabels=leafTaxonomy[rownames(ids53$counts)[1:62],"Name"],xcols=leafTaxonomy[rownames(ids53$counts)[1:62],"Color"],type="points")
dev.off()

write.table(ids53$counts[1:62,ids53$meta$spray=="Inoc"],"results/inocula53.txt")

cairo_pdf("figures/single_dropout_inocula.pdf",width=20,height=10,onefile=T,family="Arial")
par(cex=1)
plotCommunity(ids61$counts[1:62,ids61$meta$time=="t0"],xlabels=leafTaxonomy[rownames(ids61$counts)[1:62],"Name"],xcols=leafTaxonomy[rownames(ids61$counts)[1:62],"Color"],type="points")
plotCommunity(ids62$counts[1:62,ids62$meta$time=="t0"],xlabels=leafTaxonomy[rownames(ids62$counts)[1:62],"Name"],xcols=leafTaxonomy[rownames(ids62$counts)[1:62],"Color"],type="points")
dev.off()

write.table(ids61$counts[1:62,ids61$meta$time=="t0"],"results/inocula61.txt")
write.table(ids62$counts[1:62,ids62$meta$time=="t0"],"results/inocula62.txt")

icds5361 <- makeCDS(ids5361,include=list(time=c("t0")),foi="experiment",title="",legend=c("Class Dropout Inoculum","Single Dropout R1 Inoculum"))
icds5362 <- makeCDS(ids5362,include=list(time=c("t0")),foi="experiment",title="",legend=c("Class Dropout Inoculum","Single Dropout R2 Inoculum"))
icds6162 <- makeCDS(ids6162,include=list(time=c("t0")),foi="experiment",title="",legend=c("Single Dropout R1 Inoculum","Single Dropout R2 Inoculum"))
icds536162 <- makeCDS(ids536162,include=list(time=c("t0")),foi="experiment",title="",legend=c("Class Dropout Inoculum","Single Dropout R1 Inoculum","Single Dropout R2 Inoculum"))

cairo_pdf("figures/inocula_comparison.pdf",width=14,height=14,family="Arial")
par(mfrow=c(2,2),cex=1)
x <- plotPCA(icds5361,cols=bicols)
x <- plotPCA(icds6162,cols=bicols)
x <- plotPCA(icds5362,cols=bicols)
x <- plotPCA(icds536162,cols=c("red",bicols))
dev.off()

# Function to permute a pearson correlation
permute.cor <- function(x,y,n){
    creal <- cor(x,y)
    clist <- c(creal,sapply(1:n,function(i) cor(sample(x,length(x)),sample(y,length(y)))))
    p = sum(creal<clist)/n
    return(p)
}    

# Miscellaneous Correlations
# Strain abundance vs. effect size of single dropout
effects6 <- unlist(lapply(adn6,function(x) x$stats$aov.tab$R2[1]))
pvalues6 <- unlist(lapply(adn6,function(x) x$stats$aov.tab$Pr[1]))

e <- matrix(effects6[-c(1,2)],3)
rownames(e) <- c("Replicate 1","Replicate 2","Combined")
shortNames <- levels(ds6c$meta$initial)[-1]
colnames(e) <- leafTaxonomy[shortNames,]$Name
p <- matrix(pvalues6[-c(1,2)],3)

dds61 <- DESeqDataSetFromMatrix(ds61$counts[rownames(ds61$counts)!="Unclassified",ds61$meta$initial=="ALL"],ds61$meta[ds61$meta$initial=="ALL",],~1)
dds62 <- DESeqDataSetFromMatrix(ds62$counts[rownames(ds62$counts)!="Unclassified",ds62$meta$initial=="ALL"],ds62$meta[ds62$meta$initial=="ALL",],~1)
dds6162 <- DESeqDataSetFromMatrix(cbind(ds61$counts[rownames(ds61$counts)!="Unclassified",ds61$meta$initial=="ALL"],ds62$counts[rownames(ds62$counts)!="Unclassified",ds62$meta$initial=="ALL"]),rbind(ds61$meta[ds61$meta$initial=="ALL",],ds62$meta[ds62$meta$initial=="ALL",]),~1)

vst61 <- assay(varianceStabilizingTransformation(dds61))
vst62 <- assay(varianceStabilizingTransformation(dds62))
vst6162 <- assay(varianceStabilizingTransformation(dds6162))

median61 <- apply(vst61,1,median)
median62 <- apply(vst62,1,median)
median6162 <- apply(vst6162,1,median)

library(calibrate)

cairo_pdf("figures/correlations.pdf",width=7,height=7,onefile=T,family="Arial")
par(cex=1)
plot(100*e[1,],median61[colnames(fc6c)],xlab="Effect Size (%)",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="R1")
plot(100*e[1,],median61[colnames(fc6c)],xlab="Effect Size (%)",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="R1")
textxy(100*e[1,],median61[colnames(fc6c)],sub("eaf","",colnames(fc6c)))

plot(100*e[2,],median62[colnames(fc6c)],xlab="Effect Size (%)",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="R2")
plot(100*e[2,],median62[colnames(fc6c)],xlab="Effect Size (%)",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="R2")
textxy(100*e[2,],median62[colnames(fc6c)],sub("eaf","",colnames(fc6c)))

plot(100*e[3,],median6162[colnames(fc6c)],xlab="Effect Size (%)",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="Combined")
plot(100*e[3,],median6162[colnames(fc6c)],xlab="Effect Size (%)",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="Combined")
textxy(100*e[3,],median6162[colnames(fc6c)],sub("eaf","",colnames(fc6c)))

plot(mediani,median53,xlab="Inoculum Normalized Median Relative Abundance",ylab="Control Normalized Median Relative Abundance",pch=19,col=2)
plot(mediani,median53,xlab="Inoculum Normalized Median Relative Abundance",ylab="Control Normalized Median Relative Abundance",pch=19,col=2)
textxy(mediani,median53,sub("eaf","",names(median53)))

y = igraph::degree(net6)[colnames(e)]
yo = igraph::degree(net6,mode="out")[colnames(e)]
yi = igraph::degree(net6,mode="in")[colnames(e)]
x = 100*e[3,]

#plot(100*e[3,],igraph::degree(net6)[colnames(e)],xlab="Effect Size (%)",ylab="Node Degree",pch=19,col=2,main="Combined",sub=summary(lm(y~x))$r.squared)
#abline(lm(y~x))
#plot(100*e[3,],igraph::degree(net6)[colnames(e)],xlab="Effect Size (%)",ylab="Node Degree",pch=19,col=2,main="Combined")
#textxy(100*e[3,],igraph::degree(net6)[colnames(e)],sub("eaf","",rownames(leafTaxonomy)[match(colnames(e),leafTaxonomy$Name)]))
#abline(lm(y~x))
plot(x,yo,xlab="Effect Size (%)",ylab="Node Out Degree",pch=19,col=2,main="Combined",sub=summary(lm(yo~x))$r.squared)
textxy(x,yo,sub("eaf","",rownames(leafTaxonomy)[match(colnames(e),leafTaxonomy$Name)]))
abline(lm(yo~x))
#plot(100*e[3,],yi,xlab="Effect Size (%)",ylab="Node In Degree",pch=19,col=2,main="Combined",sub=summary(lm(yi~x))$r.squared)
#abline(lm(yi~x))
dev.off()

# Output some data to file
write.table(cbind(NMRA=median6162[colnames(fc6c)],ESize=100*e[3,]),"NMRA-ES.txt")

fc61 <- fcc[,seq(1,ncol(fcc),2)]
fc62 <- fcc[,1+seq(1,ncol(fcc),2)]

colnames(fc61) <- colnames(fc6c)
colnames(fc62) <- colnames(fc6c)

for(x in colnames(fc6c)){
        fc6c[x,x] <- 0
        fc61[x,x] <- 0
        fc62[x,x] <- 0
}
fccsums <- apply(fc6c,1,function(x) sum(abs(x)))
fc61sums <- apply(fc61,1,function(x) sum(abs(x)))
fc62sums <- apply(fc62,1,function(x) sum(abs(x)))

cairo_pdf("figures/correlation2.pdf",width=7,height=7,onefile=T,family="Arial")
par(cex=1)
plot(fc61sums,median61,xlab="Total Fold Changes",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="R1")
plot(fc61sums,median61,xlab="Total Fold Changes",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="R1")
textxy(fc61sums,median61,sub("eaf","",names(median61)))

plot(fc62sums,median62,xlab="Total Fold Changes",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="R1")
plot(fc62sums,median62,xlab="Total Fold Changes",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="R1")
textxy(fc62sums,median62,sub("eaf","",names(median62)))

plot(fccsums,median6162,xlab="Total Fold Changes",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="Combined")
plot(fccsums,median6162,xlab="Total Fold Changes",ylab="Normalized Median Relative Abundance",pch=19,col=2,main="Combined")
textxy(fccsums,median6162,sub("eaf","",names(median6162)))
dev.off()

# Plot the fancy figure
summary <- summariseResults(cds6c)
summary$fcMatrix <- summary$fcMatrix[-nrow(summary$fcMatrix),]
summary$pvMatrix <- summary$pvMatrix[-nrow(summary$pvMatrix),]

cairo_pdf("figures/fancy.pdf",width=14,height=14,family="Arial")
plotBipartiteSummary(summary$fcMatrix,summary$pvMatrix,leftPhylo=as.phylo(hclust(dist(t(summary$fcMatrix)),method="ward.D")),rightPhylo=as.phylo(hclust(dist(summary$fcMatrix),method="ward.D")),leftLabs=leafTaxonomy[colnames(summary$fcMatrix),]$Name,rightLabs=leafTaxonomy[rownames(summary$fcMatrix),]$Name,leftCols=leafTaxonomy[colnames(summary$fcMatrix),]$Color,rightCols=leafTaxonomy[rownames(summary$fcMatrix),]$Color,cutoff=0.01,tip.label.width=0.3)
plotBipartiteSummary(summary$fcMatrix,summary$pvMatrix,leftLabs=leafTaxonomy[colnames(summary$fcMatrix),]$Name,rightLabs=leafTaxonomy[rownames(summary$fcMatrix),]$Name,leftPhylo=keep.tip(tree,colnames(summary$fcMatrix)),rightPhylo=tree,leftCols=leafTaxonomy[colnames(summary$fcMatrix),]$Color,rightCols=leafTaxonomy[rownames(summary$fcMatrix),]$Color,cutoff=0.01,tip.label.width=0.3)
dev.off()

