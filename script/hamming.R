args=commandArgs(TRUE)
fe=args[1]
fout=args[2]
perc=args[3]
e=read.table(fe,header=T,sep="\t",stringsAsFactors=F)

hvse=function(e,p){
cutoff=quantile((e$read1+e$read2),0.999)
sumreads=e$read1+e$read2
e2=e[sumreads < cutoff,]
s1=matrix(rep(0), nrow=nrow(e2), ncol=16); s1[e2[,15:30] > p]=1
s2=matrix(rep(0), nrow=nrow(e2), ncol=16); s2[e2[,31:46] > p]=1

hamming.res=rowSums(abs(s1-s2))
boxplot(e2$entropy~hamming.res)

}

h=function(e,p){
cutoff=quantile((e$read1+e$read2),0.999)
sumreads=e$read1+e$read2
e2=e[sumreads < cutoff,]
s1=matrix(rep(0), nrow=nrow(e2), ncol=16); s1[e2[,15:30] > p]=1
s2=matrix(rep(0), nrow=nrow(e2), ncol=16); s2[e2[,31:46] > p]=1

hamming.res=rowSums(abs(s1-s2))
e2$patternDiff=hamming.res
e2
}

htable=h(e,perc)
write.table(htable, file=fout, quote=F, row.names=F, col.names=F, sep="\t")
