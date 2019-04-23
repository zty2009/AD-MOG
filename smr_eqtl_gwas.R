

gwas=read.table('gwas.txt',header=T)
mz=mean(gwas$z,na.rm=TRUE)
S=sum((gwas$z-mz)^2,na.rm=TRUE)
gwas$nz=mz+(1-(nrow(gwas)-2)/S)*(gwas$z-mz)
library(data.table)
eqtl=fread("Brain-eMeta_yangjian5e-8.txt",select=c("SNP",'Probe',"b","SE",'Gene'),data.table=F)
eqtl$z=eqtl$b/eqtl$SE
threshold1=0.05/length(unique(eqtl$Probe))
snp1=intersect(eqtl[,1],gwas[,1])

mqtl=fread("Brain-mMeta_yangjian.txt",select=c("SNP",'Probe',"b","SE",'Gene'),data.table=F)
mqtl$z=mqtl$b/mqtl$SE
threshold2=0.05/length(unique(mqtl$Probe))

snp2=intersect(mqtl[,1],gwas[,1])

for (i in 1:length(snp1)){
	index1=which(eqtl[,1]==snp[i])
	index2=which(gwas[,1]==snp[i])
		z_zx=gwas$nz[index2]
		z_zy=eqtl$z[index1]
		gene=eqtl$Gene[index1]
		for (j in 1:length(index1)){
		T=(z_zx^2*z_zy[j]^2)/(z_zx^2+z_zy[j]^2)
		P=pchisq(T,1,lower.tail=F)
		result=data.frame(snp1[i],P,gene[j])
		write.table(result, file = "SMR_eqtl_gwas.txt", append = T, col.names = F, row.names = F,quote=F)
		if (P<=threshold1){
		write.table(result, file = "SMR_eqtl_gwas_selected.txt", append = T, col.names = F, row.names = F,quote=F)
		}
		print(i)
		}
}


for (i in 1:length(snp2)){
	index1=which(mqtl[,1]==snp2[i])
	index2=which(gwas[,1]==snp2[i])
		z_zx=gwas$nz[index2]
		z_zy=mqtl$z[index1]
		gene=mqtl$Gene[index1]
		for (j in 1:length(index1)){
		T=(z_zx^2*z_zy[j]^2)/(z_zx^2+z_zy[j]^2)
		P=pchisq(T,1,lower.tail=F)
		result=data.frame(snp2[i],P,gene[j])
		write.table(result, file = "SMR_mqtl_gwas.txt", append = T, col.names = F, row.names = F,quote=F)
		if (P<=threshold2){
		write.table(result, file = "SMR_mqtl_gwas_selected.txt", append = T, col.names = F, row.names = F,quote=F)
		}
		print(i)
		}
}