#### craete the manhanttan plot with annotation 
library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(gtable)
library(RColorBrewer)
library(ggrepel)
library(gtable)
#### import the GWAS
AoD_white <- fread("/Volumes/GoogleDrive/My Drive/Aorta Diameter GWAS/GWAS_AOD/CLEANED.results_glm.AoD.all.white.gz")

### import te set of genomic loci 
genLoci <- fread("/Volumes/GoogleDrive/My Drive/Aorta Diameter GWAS/Nature_gen_revision_round2/table1_genomic_loci_with_revision.forplot.csv")

### extract only relevant info
plot.p<- AoD_white %>% select(cptid:POS,PVAL) %>% filter(PVAL<0.3)

### merge my plot with risk loci
genloci.info<- genLoci %>% select(chr,pos,genes)
plot.p<- merge(plot.p,genloci.info,by.x=c("CHR","POS"),by.y=c("chr",'pos'),all=T)

plot.p<- plot.p %>% arrange(CHR,POS) %>% mutate(row.num=1:7810039,logp=-log10(PVAL)) 

## create the subgrouping 
plot.p$CHR<- as.factor(plot.p$CHR)
chr<-plot.p %>% group_by(CHR) %>% summarise(mean= median(row.num))
chr<- chr[order(chr$mean,decreasing = F),]


##choose the color for the plot 
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,74), col=col_vector)
col2=col_vector[c(1:22)]
names(col2)<- unique(chr$CHR)



### get color and select those that represewnt my gene
manhplot<- ggplot(data=plot.p,aes(x=row.num, y=logp)) +
  geom_point(aes(col=CHR,show.legend=F)) + 
  labs(title='',x='Position',y='-log10(pvalue)') +
  geom_hline(aes(yintercept=7.3), colour="#990000", linetype="dashed") +
  scale_y_continuous(breaks=c(0,2,5,10,20,40,60,65),
                     label=c(0,2,5,10,20,40,60,""),expand = c(0.01,0.01)) +
  scale_x_continuous( breaks = as.vector(chr$mean), 
                      label = as.vector(chr$CHR),expand = c(0.01,0.01)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(fill = NA,colour = "black"),
        panel.background = element_blank(),
        axis.text.x=element_text(size=14,colour="Black"),
        axis.text.y=element_text(size=14,colour="Black")) + 
  scale_color_manual(values=col2) + theme(legend.position = "none")

#pdf(file="/Volumes/GoogleDrive/My Drive/Aorta Diameter GWAS/Manuscript_formatting_NG/final_format_figure/pre_format_figure/manhattan_plot_OAD_wite.pdf",)
manhplot + geom_text_repel(aes(label = genes),
                           size=4,min.segment.length = 0, 
                           seed = 42, box.padding = 0.5,
                           fontface = "italic")
