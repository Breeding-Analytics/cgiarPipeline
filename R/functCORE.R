report_stats<-function(objectS){
	pA <- colMeans(objectS, na.rm = TRUE) / 2 
	He_loci <- 2 * pA * (1 - pA)
	mean_He <- mean(He_loci, na.rm = TRUE)
	sd_He <- sd(He_loci, na.rm = TRUE)
	sd_val <- sd(He_loci, na.rm = TRUE)
	mean_val <- mean(He_loci, na.rm = TRUE)
	CV<-(sd_val / mean_val) * 100
	het_obs <- apply(objectS, 1, function(ind) mean(ind == 1, na.rm = TRUE))  # 1 = heterocigoto
	mean_Ho<-mean(het_obs, na.rm= TRUE)
	sd_Ho<-sd(het_obs, na.rm= TRUE)
	shannon_index <- function(p) {
	p <- p[p > 0]  # eliminar ceros para evitar log(0)
	-sum(p * log(p))
	}
	shannon_values <- numeric(length(pA))
	for (i in seq_along(pA)) {
	freqs <- c(pA[i], 1 - pA[i])  # alelo alt y ref
	shannon_values[i] <- shannon_index(freqs)
	}
	mean_shannon <- mean(shannon_values, na.rm = TRUE)
	sd_shannon <- sd(shannon_values, na.rm = TRUE)
	return(c(mean_He,sd_He,mean_Ho,sd_Ho,CV,mean_shannon,sd_shannon))
}

funchunt<-function(datosgen,dir_fileGen,dir_filePhen,dir_fileDist,size1,w){
library(corehunter)
dist.file<-as.matrix(dist(datosgen))
mds_res <- cmdscale(dist.file, k = 3, eig = TRUE) 
perctCP12=c(round(mds_res$eig[1]/sum(mds_res$eig)*100,2),round(mds_res$eig[2]/sum(mds_res$eig)*100,2),round(mds_res$eig[3]/sum(mds_res$eig)*100,2))

if (dir_fileGen[1]!="none" & dir_filePhen[1]=="none" & dir_fileDist[1]=="none"){
   geno.file=datosgen
   my.data <- genotypes(geno.file,format="biparental")
}

if (dir_fileGen[1]!="none" & dir_filePhen[1]!="none" & dir_fileDist[1]=="none"){
    geno.file=datosgen
    pheno.file=read.autodelim(dir_filePhen$datapath)
    my.data <- coreHunterData(genotypes(geno.file,format="biparental"),phenotypes(pheno.file))
}

if (dir_fileGen[1]!="none" & dir_filePhen[1]!="none" & dir_fileDist[1]!="none"){
   geno.file=datosgen
   pheno.file=read.autodelim(dir_filePhen$datapath)
   dist.file=read.autodelim(dir_fileDist$datapath)
   my.data <- coreHunterData(genotypes(geno.file,format="biparental"),phenotypes(pheno.file),distances(dist.file))
}

if (dir_fileGen[1]=="none" & dir_filePhen[1]!="none" & dir_fileDist[1]=="none"){
   pheno.file=read.autodelim(dir_filePhen$datapath)
   my.data <- phenotypes(pheno.file)
}

if (dir_fileGen[1]=="none" & dir_filePhen[1]=="none" & dir_fileDist[1]!="none"){
   #dist.file=read.autodelim(dir_fileDist$datapath)   
   my.data <- distances(dist.file)
}


values=w
objall=list(objective("EN", "MR", weight=as.numeric(w[1])),
            objective("EN", "CE", weight=as.numeric(w[2])),
            objective("EN", "GD", weight=as.numeric(w[3])),
            objective("EN", "PD", weight=as.numeric(w[4])),
            objective("AN", "MR", weight=as.numeric(w[5])),
            objective("AN", "CE", weight=as.numeric(w[6])),
            objective("AN", "GD", weight=as.numeric(w[7])),
            objective("AN", "PD", weight=as.numeric(w[8])),
            objective("EE", "MR", weight=as.numeric(w[9])),
            objective("EE", "CE", weight=as.numeric(w[10])),
            objective("EE", "GD", weight=as.numeric(w[11])),
            objective("EE", "PD", weight=as.numeric(w[12])),
            objective("SH", weight=as.numeric(w[13])),
            objective("HE", weight=as.numeric(w[14])),
            objective("CV", weight=as.numeric(w[15])))
objuse=list()
use=which(values!=0)
for(i in 1:length(use)){
objuse[[i]]=objall[[use[i]]]
}

core1=sampleCore(my.data,obj=objuse,size=as.numeric(size1),indices=TRUE)

Gen=core1$sel
#GenSelect=rownames(datosGen)[Gen]
mds_res<-mds_res$points
mds_res<-data.frame(rownames(mds_res),mds_res)
mds_res$CoreSet<-"Pop"
mds_res$CoreSet[Gen]<-"Selected"
colnames(mds_res)=c("Gen","Factor1","Factor2","Factor3","CoreSet")

statsF<-data.frame(cbind(c("He","sd_He","Ho","sd_Ho","CV","SH_Index","sd_SH"),rep(c("Pop_Complete","Core_Selected"),each=7),c(report_stats(datosgen),report_stats(datosgen[Gen,]))))
colnames(statsF)<-c("Parameter","Environment","Value")

 return(list(mds_res,statsF,perctCP12))
}
