require(data.table)
require(cowplot)
library(tidyverse)
library(scales)
library(Biostrings)
library(edgeR)
library(limma)
library(hues)
library(viridis)
library(impute)
library(ggrepel)

options(stringsAsFactors=F)

dir.create("plots",showWarnings=F)
dir.create("tables",showWarnings=F)

ref_fasta=readDNAStringSet("fasta/reference.fa")
ref_name="Csde1_165_414"
ref_start=165
ref_end=414
primer_f="TGCTTCAAGTTCAGATCAGGCAAGG"
primer_r="AAAGGAGAGCGAGAGGAAATGATCTACC"
data_path="./out/Csde1_165_414"
ext="raw"

#Reference sequence position and nucleotide info
nuc_pos=data.frame(pos=1:length(ref_fasta[[ref_name]])+ref_start-1,nuc_dna=unlist(strsplit(toupper(as.character(ref_fasta[[ref_name]])),"")))
nuc_pos$nuc_pos=paste(nuc_pos$pos,nuc_pos$nuc_dna,sep="")
nuc_pos$nuc_rna=nuc_pos$nuc_dna
nuc_pos$nuc_rna[nuc_pos$nuc_dna=="T"]="U"
nuc_pos$nuc_pos_rna=paste(nuc_pos$pos,nuc_pos$nuc_rna,sep="")
nuc_pos$mask=F
nuc_pos$mask[1:nchar(primer_f)]=T #Primer region, fwd
nuc_pos$mask[(nrow(nuc_pos)-nchar(primer_r)+1):nrow(nuc_pos)]=T #Primer region, rev
rownames(nuc_pos)=nuc_pos$nuc_pos

#Read in count matrices
files=list.files(data_path,pattern=glob2rx(paste("*",ref_name,"*.",ext,sep="")),full.names=T)
counts=list()
for (file in files)
{
  count=as.matrix(data.table::fread(file,sep="\t",header=F))
  counts[[file]]=count
}
counts=as.matrix(t(Reduce('+', counts)))

rownames(counts)=nuc_pos$nuc_pos
colnames(counts)=nuc_pos$nuc_pos

nuc_pos_trimmed=subset(nuc_pos,!mask)$nuc_pos #Masks primer regions

#TMM normalization
expr=DGEList(counts=counts[nuc_pos_trimmed,nuc_pos_trimmed])
expr=calcNormFactors(expr,method="TMM")
v=voom(expr,plot=F,span=0.5,normalize.method="none") #Could use this to handle replicates
ve=t(v$E)

#Check mean, variance ranges
pos_mean_var=cbind(nuc_pos[nuc_pos_trimmed,],data.frame(col_mean=apply(ve,2,mean),
                                                     col_var=apply(ve,2,var),
                                                     row_mean=apply(ve,1,mean),
                                                     row_var=apply(ve,1,var)))

#Row vs. col variances
fig_rowvar_colvar=ggplot(pos_mean_var,aes(x=col_var,y=row_var,colour=nuc_dna))+
  theme_classic()+geom_point(shape=16,stroke=0,alpha=0.5,size=2)+
  xlab("Column variance")+ylab("Row variance")+scale_colour_iwanthue(name=NULL)
cairo_pdf("plots/rowvar_colvar.pdf",width=4,height=3)
print(fig_rowvar_colvar)
dev.off()

#Col mean vs col variances
fig_colmean_colvar=ggplot(pos_mean_var,aes(x=col_mean,y=col_var,colour=nuc_dna))+
  theme_classic()+geom_point(shape=16,stroke=0,alpha=0.5,size=2)+
  xlab("Column mean")+ylab("Column variance")+scale_colour_iwanthue(name=NULL)
cairo_pdf("plots/colmean_colvar.pdf",width=4,height=3)
print(fig_colmean_colvar)
dev.off()

#Filter low signal columns
colmean_cutoff=10.5 #Based on the variane plots
nuc_pos_filter=subset(pos_mean_var,col_mean<colmean_cutoff)$nuc_pos
nuc_pos$filtered=F
nuc_pos$filtered[nuc_pos$mask]=NA
nuc_pos$filtered[nuc_pos$nuc_pos%in%nuc_pos_filter]=T
nuc_pos$used=(!nuc_pos$mask)&(!nuc_pos$filtered)
nuc_pos_used=nuc_pos$nuc_pos[nuc_pos$used]

vef=ve
vef[nuc_pos_filter,]=NA
vef[,nuc_pos_filter]=NA
vef_nona=vef[!apply(is.na(vef),1,all),!apply(is.na(vef),2,all)] #Drop NA col/row

vef_long=reshape2::melt(vef_nona) #Long form for plotting
colnames(vef_long)=c("row_nuc_pos","col_nuc_pos","value")
vef_long$row_nuc_pos=factor(as.character(vef_long$row_nuc_pos),levels=nuc_pos_used)
vef_long$col_nuc_pos=factor(as.character(vef_long$col_nuc_pos),levels=rev(nuc_pos_used))
vef_long$row_nuc=nuc_pos[as.character(vef_long$row_nuc_pos),"nuc_dna"]
vef_long$col_nuc=nuc_pos[as.character(vef_long$col_nuc_pos),"nuc_dna"]
vef_long$row_pos=nuc_pos[as.character(vef_long$row_nuc_pos),"pos"]
vef_long$col_pos=nuc_pos[as.character(vef_long$col_nuc_pos),"pos"]
vef_long[(vef_long$row_pos>=(vef_long$col_pos-4))&(vef_long$row_pos<=(vef_long$col_pos+4)),]$value=NA #Mask +/- 4nt around diagonal

#Imputate NA values
vef_long_casted=reshape2::dcast(vef_long[,c("row_nuc_pos","col_nuc_pos","value")],row_nuc_pos~col_nuc_pos)
rownames(vef_long_casted)=vef_long_casted[,1]
vef_long_casted=as.matrix(vef_long_casted[,-1])
vef_long_casted=vef_long_casted[,rev(colnames(vef_long_casted))]
vef_long_casted_imputed=t(impute.knn(t(vef_long_casted))$data)
vef_long_imputed=reshape2::melt(vef_long_casted_imputed) #Long from for plotting
colnames(vef_long_imputed)=c("row_nuc_pos","col_nuc_pos","value")
vef_long_imputed$row_nuc_pos=factor(as.character(vef_long_imputed$row_nuc_pos),levels=nuc_pos_used)
vef_long_imputed$col_nuc_pos=factor(as.character(vef_long_imputed$col_nuc_pos),levels=rev(nuc_pos_used))
vef_long_imputed$row_nuc=nuc_pos[as.character(vef_long_imputed$row_nuc_pos),"nuc_dna"]
vef_long_imputed$col_nuc=nuc_pos[as.character(vef_long_imputed$col_nuc_pos),"nuc_dna"]
vef_long_imputed$row_pos=nuc_pos[as.character(vef_long_imputed$row_nuc_pos),"pos"]
vef_long_imputed$col_pos=nuc_pos[as.character(vef_long_imputed$col_nuc_pos),"pos"]

#Z-scaling
vefs=apply(vef,2,scale)
rownames(vefs)=colnames(v)
colnames(vefs)=rownames(v)
vefs=vefs-median(vefs,na.rm=T) #Center on median
vefs_nona=vefs[!apply(is.na(vefs),1,all),!apply(is.na(vefs),2,all)] #Drop NA col/row

vefs_long=reshape2::melt(vefs_nona) #Long form for plotting
colnames(vefs_long)=c("row_nuc_pos","col_nuc_pos","value")
vefs_long$row_nuc_pos=factor(as.character(vefs_long$row_nuc_pos),levels=nuc_pos_used)
vefs_long$col_nuc_pos=factor(as.character(vefs_long$col_nuc_pos),levels=rev(nuc_pos_used))
vefs_long$row_nuc=nuc_pos[as.character(vefs_long$row_nuc_pos),"nuc_dna"]
vefs_long$col_nuc=nuc_pos[as.character(vefs_long$col_nuc_pos),"nuc_dna"]
vefs_long$row_pos=nuc_pos[as.character(vefs_long$row_nuc_pos),"pos"]
vefs_long$col_pos=nuc_pos[as.character(vefs_long$col_nuc_pos),"pos"]
vefs_long[(vefs_long$row_pos>=(vefs_long$col_pos-4))&(vefs_long$row_pos<=(vefs_long$col_pos+4)),]$value=NA #Mask +/- 4nt around diagonal

#Plot heatmap of z-scaled accessibility change matrix
#Fig. 5d in manuscript
hm_labelcoords=subset(nuc_pos,used&(pos%in%c(195,215,265,315,330,365,385)))$nuc_pos #axis labels
hm_labeltexts=as.character(subset(nuc_pos,used&(pos%in%c(195,215,265,315,330,365,385)))$pos)
fig_vefs=ggplot(vefs_long,
               aes(x=row_nuc_pos,y=col_nuc_pos,colour=value,fill=value))+
  theme_classic()+coord_fixed()+geom_tile()+
  scale_colour_viridis(limits=c(-2,2),oob=squish,option="magma",na.value="grey80")+
  scale_fill_viridis(limits=c(-2,2),oob=squish,option="magma",na.value="grey80")+
  geom_rect(xmin="215C",xmax="315T",ymin="215C",ymax="315T",fill=NA,colour="white")+
  geom_rect(xmin="334T",xmax="363A",ymin="334T",ymax="363A",fill=NA,colour="white")+
  scale_x_discrete(breaks=hm_labelcoords,labels=hm_labeltexts)+
  scale_y_discrete(breaks=rev(hm_labelcoords),labels=rev(hm_labeltexts))+
  xlab("Nucleotide positions")+ylab("Nucleotide positions")+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5),
        axis.text.y=element_text(hjust=0.5,vjust=0.5))
cairo_pdf("plots/Csde1_Z_Incell.pdf",width=5,height=5)
print(fig_vefs)
dev.off()

#Imputate NA values for z-scaled matrix
vefs_long_casted=reshape2::dcast(vefs_long[,c("row_nuc_pos","col_nuc_pos","value")],row_nuc_pos~col_nuc_pos)
rownames(vefs_long_casted)=vefs_long_casted[,1]
vefs_long_casted=as.matrix(vefs_long_casted[,-1])
vefs_long_casted=vefs_long_casted[,rev(colnames(vefs_long_casted))]
vefs_long_casted_imputed=t(impute.knn(t(vefs_long_casted))$data)
vefs_long_imputed=reshape2::melt(vefs_long_casted_imputed) #Long from for plotting
colnames(vefs_long_imputed)=c("row_nuc_pos","col_nuc_pos","value")
vefs_long_imputed$row_nuc_pos=factor(as.character(vefs_long_imputed$row_nuc_pos),levels=nuc_pos_used)
vefs_long_imputed$col_nuc_pos=factor(as.character(vefs_long_imputed$col_nuc_pos),levels=rev(nuc_pos_used))
vefs_long_imputed$row_nuc=nuc_pos[as.character(vefs_long_imputed$row_nuc_pos),"nuc_dna"]
vefs_long_imputed$col_nuc=nuc_pos[as.character(vefs_long_imputed$col_nuc_pos),"nuc_dna"]
vefs_long_imputed$row_pos=nuc_pos[as.character(vefs_long_imputed$row_nuc_pos),"pos"]
vefs_long_imputed$col_pos=nuc_pos[as.character(vefs_long_imputed$col_nuc_pos),"pos"]

#Variance per position after all masking
pos_mean_var_used=cbind(nuc_pos[nuc_pos_used,],
                        data.frame(col_mean=apply(vefs_long_casted,2,mean,na.rm=T),
                                   col_var=apply(vefs_long_casted,2,var,na.rm=T),
                                   row_mean=apply(vefs_long_casted,1,mean,na.rm=T),
                                   row_var=apply(vefs_long_casted,1,var,na.rm=T)))

#Correlation with "wild-type" mutants
#Fig. 5f in manuscript
wt_rows=head(nuc_pos_used[order(pos_mean_var_used$row_var)],10)
wt_row_mean=apply(vef_long_casted[wt_rows,],2,mean,na.rm=T)
pos_mean_var_used$wt_cor=apply(vef_long_casted,1,cor,wt_row_mean,use="complete",method="spearman")
wt_cor_pval=wilcox.test(subset(pos_mean_var_used,(pos>=215)&(pos<=315))$wt_cor,subset(pos_mean_var_used,!(pos>=215)&(pos<=315))$wt_cor)$p.value
fig_wt_cor=ggplot(pos_mean_var_used,aes(x=nuc_pos,y=wt_cor))+theme_classic()+
  geom_rect(data=data.frame(c(1)),xmin="215C",xmax="315T",ymin=-Inf,ymax=Inf,inherit.aes=F,colour=NA,fill="grey90",alpha=1)+geom_point()+
  xlab("Nucleotide position")+
  ylab("Accessbility correlation with\nlow-variance (\"wild-type\") mutants")+
  scale_x_discrete(breaks=hm_labelcoords,labels=hm_labeltexts)+
  annotate(geom="text",x="365C",y=0.8,label=paste("p=",format(wt_cor_pval,digits=2),sep=""))
cairo_pdf("plots/corr_wt.pdf",width=6,height=3)
print(fig_wt_cor)
dev.off()

#Multidimensional scaling
mds=cmdscale(dist(vefs_long_casted_imputed),k=3) #Default R dist is Euclidean, so this is =PCA
mds=data.frame(x=mds[,1],y=mds[,2],z=mds[,3])
mds=cbind(nuc_pos[nuc_pos_used,],mds)

#Pick clusters; could add something auto
cluster1=rownames(subset(mds,(x<=-3)&((-z)>=0))) #Based on x-z scatter plot
cluster2=rownames(subset(mds,(x>=4)&((-z)<=-0.5)))

mds$clusters=""
mds$clusters[rownames(mds)%in%cluster1]="Cluster 1"
mds$clusters[rownames(mds)%in%cluster2]="Cluster 2"
mds$clusters=factor(mds$clusters,levels=c("Cluster 1","Cluster 2",""))

#Fig. 6a in manuscript
fig_pcoa=ggplot(mds,aes(x=x,y=z))+theme_classic()+
  geom_point(aes(colour=clusters),alpha=1,size=1)+
  geom_text_repel(aes(label=clusters,colour=clusters),force=25,size=3)+
  xlab("Dimension 1")+ylab("Dimension 2")+
  scale_colour_manual(values=c("#aa0000","#0066ff","black"))+
  scale_x_continuous(expand=c(0.1,0.1))+scale_y_continuous(expand=c(0.1,0.1))+
  theme(legend.position="none",axis.text=element_text(size=12))
cairo_pdf("plots/pcoa.pdf",width=5,height=4)
print(fig_pcoa)
dev.off()

#Mean cluster diff. accesibility 
regb_start=215
regb_end=315
regb_len=regb_end-regb_start+1

cluster_avg=data.frame(vefs_long %>% group_by(col_nuc_pos) %>% summarise(
  cluster1=mean(value[row_nuc_pos%in%cluster1],na.rm=T),
  cluster2=mean(value[row_nuc_pos%in%cluster2],na.rm=T)))
rownames(cluster_avg)=cluster_avg$col_nuc_pos
cluster_avg=cbind(nuc_pos[nuc_pos_used,],cluster_avg[nuc_pos_used,])
cluster_avg_regb=subset(cluster_avg,(pos>=regb_start)&(pos<=regb_end))

#1D data for Csde1 5'UTR, from ATP depletion experiments
csde1_atpd=data.frame(data.table::fread("./csde1_1d.tsv",sep="\t",header=T))
rownames(csde1_atpd)=paste(csde1_atpd$pos,csde1_atpd$nuc,sep="")
cluster_avg_regb$atpd=csde1_atpd[cluster_avg_regb$nuc_pos,]$mean_atp_not
cluster_avg_regb$ivf=csde1_atpd[cluster_avg_regb$nuc_pos,]$mean_ivf

#Plot comparison with 1D data
#Fig. 6f in the manuscript
cluster_avg_regb_1d_long=reshape2::melt(cluster_avg_regb[,c("pos","cluster1","cluster2","atpd","ivf")],id.vars="pos")
cluster_avg_regb_1d_labels=as.character(c(215,235,255,290,315))
f1=ggplot(subset(cluster_avg_regb_1d_long,variable=="cluster1"),aes(x=factor(pos),y=value))+
  theme_classic()+geom_hline(yintercept=0,linetype="dashed")+
  geom_bar(stat="identity",fill="#aa0000")+
  ylab("Reactivity difference")+xlab("Nucleotide position")+
  scale_x_discrete(expand=c(0,0),breaks=cluster_avg_regb_1d_labels)
f2=ggplot(subset(cluster_avg_regb_1d_long,variable=="cluster2"),aes(x=factor(pos),y=value))+
  theme_classic()+geom_hline(yintercept=0,linetype="dashed")+
  geom_bar(stat="identity",fill="#0066ff")+
  ylab("Reactivity difference")+xlab("Nucleotide position")+
  scale_x_discrete(expand=c(0,0),breaks=cluster_avg_regb_1d_labels)
f3=ggplot(subset(cluster_avg_regb_1d_long,variable=="atpd"),aes(x=factor(pos),y=value))+
  theme_classic()+geom_hline(yintercept=0,linetype="dashed")+
  geom_bar(stat="identity",fill="#800080")+
  ylab("Reactivity difference")+xlab("Nucleotide position")+
  scale_x_discrete(expand=c(0,0),breaks=cluster_avg_regb_1d_labels)
fig_atpd_clusters_acc=cowplot::plot_grid(f3,f1,f2,ncol=1,align="v")
cairo_pdf("plots/atpd_clusters_acc.pdf",width=6,height=5)
print(fig_atpd_clusters_acc)
dev.off()

#Output constraints for Vienna fold
#use this as --command=$constraint_file with RNAsubopt
cluster_avg_regb_viennacmd_cluster1=data.frame(e="E",pos=cluster_avg_regb$pos-regb_start+1,k=0,j=1,cluster1=cluster_avg_regb$cluster1)
cluster_avg_regb_viennacmd_cluster2=data.frame(e="E",pos=cluster_avg_regb$pos-regb_start+1,k=0,j=1,cluster1=cluster_avg_regb$cluster2)
write.table(cluster_avg_regb_viennacmd_cluster1,file="tables/bonus1.tsv",sep="\t",row.names=F,col.names=F,quote=F)
write.table(cluster_avg_regb_viennacmd_cluster2,file="tables/bonus2.tsv",sep="\t",row.names=F,col.names=F,quote=F)

#Output RDAT for use with REEFFIT
nuc_pos_regb=subset(nuc_pos,(pos>=regb_start)&(pos<=regb_end))
nuc_pos_regb$nuc_rna_comp=unlist(strsplit(as.character(complement(RNAString(paste(nuc_pos_regb$nuc_rna,collapse="")))),""))
nuc_pos_used_regb=subset(nuc_pos,(pos>=regb_start)&(pos<=regb_end)&used)$nuc_pos

ver=ve
ver[nuc_pos_filter,]=NA
ver_nona=ver[!apply(is.na(ver),1,all),!apply(is.na(ver),2,all)] #Drop NA col/row

ver_long=reshape2::melt(ver_nona) #Long form for plotting
colnames(ver_long)=c("row_nuc_pos","col_nuc_pos","value")
ver_long$row_nuc_pos=factor(as.character(ver_long$row_nuc_pos),levels=nuc_pos_used)
ver_long$col_nuc_pos=factor(as.character(ver_long$col_nuc_pos),levels=rev(nuc_pos_trimmed))
ver_long$row_nuc=nuc_pos[as.character(ver_long$row_nuc_pos),"nuc_dna"]
ver_long$col_nuc=nuc_pos[as.character(ver_long$col_nuc_pos),"nuc_dna"]
ver_long$row_pos=nuc_pos[as.character(ver_long$row_nuc_pos),"pos"]
ver_long$col_pos=nuc_pos[as.character(ver_long$col_nuc_pos),"pos"]
ver_long[(ver_long$col_pos>=(ver_long$row_pos-4))&(ver_long$col_pos<=(ver_long$row_pos+4)),]$value=NA #Mask +/- 4nt around diagonal

#Scale within each nucleotide
ver_long[ver_long$col_nuc=="A",]$value=scale(ver_long[ver_long$col_nuc=="A",]$value)
ver_long[ver_long$col_nuc=="T",]$value=scale(ver_long[ver_long$col_nuc=="T",]$value)
ver_long[ver_long$col_nuc=="C",]$value=scale(ver_long[ver_long$col_nuc=="C",]$value)
ver_long[ver_long$col_nuc=="G",]$value=scale(ver_long[ver_long$col_nuc=="G",]$value)

#Imputate NA values
ver_long_casted=reshape2::dcast(ver_long[,c("row_nuc_pos","col_nuc_pos","value")],row_nuc_pos~col_nuc_pos)
rownames(ver_long_casted)=ver_long_casted[,1]
ver_long_casted=as.matrix(ver_long_casted[,-1])
ver_long_casted=ver_long_casted[,rev(colnames(ver_long_casted))]
ver_long_casted_imputed=t(impute.knn(t(ver_long_casted))$data)

ver_long_imputed=reshape2::melt(ver_long_casted_imputed) #Long from for plotting
colnames(ver_long_imputed)=c("row_nuc_pos","col_nuc_pos","value")
ver_long_imputed$row_nuc_pos=factor(as.character(ver_long_imputed$row_nuc_pos),levels=nuc_pos_used)
ver_long_imputed$col_nuc_pos=factor(as.character(ver_long_imputed$col_nuc_pos),levels=rev(nuc_pos_used))
ver_long_imputed$row_nuc=nuc_pos[as.character(ver_long_imputed$row_nuc_pos),"nuc_dna"]
ver_long_imputed$col_nuc=nuc_pos[as.character(ver_long_imputed$col_nuc_pos),"nuc_dna"]
ver_long_imputed$row_pos=nuc_pos[as.character(ver_long_imputed$row_nuc_pos),"pos"]
ver_long_imputed$col_pos=nuc_pos[as.character(ver_long_imputed$col_nuc_pos),"pos"]


rdat_mat=matrix(NaN,nrow=nrow(nuc_pos_regb),ncol=nrow(nuc_pos_regb))
rownames(rdat_mat)=nuc_pos_regb$nuc_pos
colnames(rdat_mat)=nuc_pos_regb$nuc_pos

rdat_mat=rbind(NaN,rdat_mat)
rownames(rdat_mat)[1]="WT"

rdat_mat[c("WT",nuc_pos_used_regb),nuc_pos_regb$nuc_pos]=
  2^(rbind(apply(ver_long_casted_imputed[wt_rows,],2,mean,na.rm=T)[nuc_pos_regb$nuc_pos],
           ver_long_casted_imputed[nuc_pos_used_regb,nuc_pos_regb$nuc_pos])*2-2.1) #Tranformation based on comparison with data deposited in RMDB

#File formatting
rdat_mat_badqual=rep("",nrow(rdat_mat))
rdat_mat_badqual[which(apply(is.nan(rdat_mat),1,all))]="\twarning:badQuality"

rdat_mat_out=cbind(paste("REACTIVITY:",1:nrow(rdat_mat),sep=""),rdat_mat)

colnames(rdat_mat_out)=c("SEQPOS",
                         paste(subset(nuc_pos,(pos>=regb_start)&(pos<=regb_end))$nuc_rna,
                               subset(nuc_pos,(pos>=regb_start)&(pos<=regb_end))$pos
                               -regb_start+1,sep=""))

rdat_annot_out=data.frame(cbind(paste("ANNOTATION_DATA:",1:(nrow(rdat_mat_out)),sep=""),
                                paste("mutation:",
                                      c("WT",
                                        paste(nuc_pos_regb$nuc_rna,nuc_pos_regb$pos-regb_start+1,
                                              nuc_pos_regb$nuc_rna_comp,sep="")),
                                      rdat_mat_badqual,sep="")))

regb_seq=paste(nuc_pos_regb$nuc_rna,collapse="")
rdat_header_out=rbind("RDAT_VERSION\t0.34",paste("NAME\t","Csde1_",regb_start,"_",regb_end,sep=""),paste("SEQUENCE\t",regb_seq,sep=""))

#These 3 need to be cat'ed to .rdat file
write.table(rdat_header_out,file=paste("tables/Csde1_",regb_start,"_",regb_end,".seq",sep=""),
            sep="\t",row.names=F,col.names=F,quote=F) #RNA sequence
write.table(rdat_annot_out,file=paste("tables/Csde1_",regb_start,"_",regb_end,".an",sep=""),
            sep="\t",row.names=F,col.names=F,quote=F) #RNA sequence
write.table(format(rdat_mat_out,digits=4),
            file=paste("tables/Csde1_",regb_start,"_",regb_end,".reac",sep=""),
            sep="\t",row.names=F,col.names=T,quote=F,na="nan")

#REEFFIT params
#--structfile bonus_combined.dot
#--decompose motif
#--cstart 16 --cend 75 --expcl 3

