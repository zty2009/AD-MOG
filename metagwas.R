library(data.table)
gwas1=fread("IG.txt",select=c("MarkerName","Beta","SE",'Pvalue'),data.table=F)
gwas2=fread("UK",select=c("SNP","BETA","SE",'P'),data.table=F)
names(gwas2)=c("MarkerName","Beta","SE",'Pvalue')
gwas3=fread("co.linear",select=c("SNP","BETA","SE",'P'),data.table=F)
names(gwas3)=c("MarkerName","Beta","SE",'Pvalue')
gwas4=fread("AD.linear",select=c("SNP","BETA","SE",'P'),data.table=F)
names(gwas4)=c("MarkerName","Beta","SE",'Pvalue')
gwas5=fread("bi.linear",select=c("SNP","BETA","SE",'P'),data.table=F)
names(gwas5)=c("MarkerName","Beta","SE",'Pvalue')
data=rbind(gwas1,gwas2,gwas3,gwas4,gwas5)
data=data[data[,4]<0.05,]
name=unique(data[,1])
result1=data.frame()
for (i in 1:length(name)) {
	print(i)
		index=which(data[,1]==name[i])
		se=(1/sum(1/data[index,3]^2))^0.5
		beta=sum(data[index,2]*1/data[index,3]^2)/sum(1/data[index,3]^2)
		Z=beta/se
		P=2*pnorm(abs(Z),lower.tail=F)
		result=c(name[i],Z,P)
		result1=rbind(result1,result)
	}

	
		
