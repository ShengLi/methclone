args=commandArgs(TRUE)
fe=args[1]
fs1=args[2]
fs2=args[3]
fout=args[4]
e=read.table(fe,header=T,sep="\t",stringsAsFactors=F)
s1=read.table(fs1,header=F,sep="\t",stringsAsFactors=F)
s2=read.table(fs2,header=F,sep="\t",stringsAsFactors=F)
n=length(strsplit(e$loci[1],":")[[1]])
info=matrix(unlist(strsplit(e$loci,":")), ncol=n,byrow=T)
s.to.exclude=sort(paste(c(s1$V1,s2$V1), c(s1$V2,s2$V2), sep=":"))
m=matrix(rep(0),ncol=ncol(info), nrow=nrow(info))
for(i in 1:ncol(info)) m[which(paste(e$chr, info[,i],sep=":") %in% s.to.exclude),i]=1
id.to.exclude=which(rowSums(m)>0)
if(length(id.to.exclude)!=0){
out=e[-id.to.exclude,]
write.table(out, file=fout, sep="\t",row.names=F, col.names=T, quote=F)
} else {
print("No RRBS covered loci output by methclone overlaped with the supplied vcf files.")
}
