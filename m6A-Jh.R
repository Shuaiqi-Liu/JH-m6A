library(data.table)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(AnnotationDbi)
library(AnnotationHub)
library(clusterProfiler)
library(ggsci)
library(ggvenn)
library(ggbreak)

silk_gene <- fread("silkgene.csv")
silk_gene$id <- silk_gene$id %>% as.character()


#################DEG#######
silkdeg <- fread("jh_20e_treated_data_fpkm.csv") %>%
  .[Gene_ID %in% silk_gene$id,c(1,17:19,38:40)] %>%
  .[,Gene_ID := as.character(Gene_ID)] %>%
  silk_gene[., on = .(id = Gene_ID)] %>% 
  setorder(group)

pheatmap(log2(deg[1,5:10]+1),cluster_rows = FALSE,
         cluster_cols = FALSE,labels_row = deg$name,
         show_colnames = FALSE,cellheight = 10,cellwidth = 13,
         gaps_col = 3,gaps_row = c(15,21))
pheatmap(log2(deg[1,5:10]+1),cluster_rows = FALSE,
         cluster_cols = FALSE,labels_row = deg$name,
         show_colnames = FALSE,cellheight = 10,cellwidth = 13,
         gaps_col = 3)

hub <- AnnotationHub()
orgdb <- hub[["AH102039"]]

depeak <- fread("JH1VSC1_diffPeak_result.csv")
deg <- fread("JH1_inputVSC1_input_Gene_differential_expression.csv") %>%
  .[,gene_id := as.character(gene_id)]
deg
dou_dif <- deg %>% .[gene_id %in% depeak$geneID & pval < 0.05,]
dou_dif$gene_id %>% uniqueN()

write.csv(dou_dif, file= 'dou_diff.csv',row.names = FALSE)
enrichKEGG(dou_dif$gene_id,organism = "bmor",
           pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 1) %>%
  dotplot(orderBy = "y")+
  theme()
kegg_dot <- enrichKEGG(dou_dif$gene_id,organism = "bmor",
           pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 1) %>%
  as.data.frame() %>% setorder(.,-Count) %>%.[1:8,] %>%
  setorder(.,-Description) %>% .[] 
kegg_dot$Description <- factor(kegg_dot$Description,
                               levels = kegg_dot$Description[order(kegg_dot$Description,decreasing = TRUE)])
kegg_dot
ggplot(kegg_dot,aes(y = Description, x = GeneRatio))+
  geom_point(aes(size = Count,color = pvalue))+
  theme(axis.text.y  = element_blank(),
           axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  theme_bw()
  
###########JH signaling pathway analysis #########

fread("C1_IPVSC1_input_peak.csv") %>% .[geneID %in% fread("jh通路基因.txt")[,id] & log2FoldChange > 2 & pvalue < 0.05,c(1:15)] %>% 
  write.csv(.,file = "JH_signaling_pathway_with_m6A.csv",row.names = FALSE)

fread("m3_rnai.csv")%>% .[Gene_ID %in% fread("jh通路基因.txt")[,id] &`Qvalue_(siM3_/_NC)` < 0.05,] %>% 
  write.csv(.,file = "JH_siM3.csv")

kegg_san <- fread("kegg_dou_dif.csv") %>% .[,Description:= paste0(Description, "(",ID,")")] %>%
  .[,2:17] %>% melt(.,id.vars = "Description",measure.vars = .SD,value.name = "gene_id") %>% 
  .[,c(1,3)] %>% na.omit() %>% .[,gene_id := as.character(gene_id)]

kegg_san <- deg[,c(5,19)] %>% .[kegg_san,on = .(gene_id = gene_id)]
kegg_san[,value := 1]
kegg_san <- setorder(kegg_san,Description) %>% .[]
kegg_san %>% uniqueN(.,by = "gene_id")

library(ggalluvial)
ggplot(kegg_san,aes(axis1=gene_id,
                  axis2=Description,fill = regulation))+
  geom_flow(aes(fill=regulation),
            color = "steelblue",alpha = .7)+
  geom_stratum(width = .3)+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size = 1)+
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank())
  
########four-quadrant plot###

jh_silk <- fread("JH_20E_treat.txt") %>% .[,c(1,14,23)]

colnames(jh_silk) <- colnames(jh_silk) %>% gsub(" ","",.)
jh_silk[, GeneID := as.character(GeneID)]

vol_dt <- jh_silk %>% depeak[,c(8,9,13)][.,on = .(geneID = GeneID)] %>% na.omit()

vol_dt <- vol_dt %>% .[`Qvalue(JH_psg/control_psg)` < 0.05,group := "sign"] %>% 
  .[`Qvalue(JH_psg/control_psg)` >= 0.05,group := "nosign"] %>% unique(., by = "geneID")
ggplot(vol_dt,aes(x = `DiffModLog2FC`, y = `log2(JH_psg/control_psg)`))+
  geom_point(aes(color= group),size = .7,alpha = .8)+
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(family = "Times"))+
  ylim(-5,5)+
  geom_hline(yintercept = c(-1,1),linetype = "dotted",size = .7)+
  geom_vline(xintercept = c(-1,1),linetype = "dotted",size = .7)+
  scale_color_manual(values = c("grey","steelblue"))
 depeak %>% names()

ven_degup <- jh_silk[`Qvalue(JH_psg/control_psg)` < 0.05 & `log2(JH_psg/control_psg)` > 1,GeneID]

ven_degdown <- jh_silk[`Qvalue(JH_psg/control_psg)` < 0.05 & `log2(JH_psg/control_psg)` < 1,GeneID]

ven_meup <- depeak[DiffModLog2FC > 1, geneID] %>% unique(.,by = geneID)
ven_medown <- depeak[DiffModLog2FC < 1 , geneID] %>% unique(.,by = geneID)

ven <- list(ven_degdown,ven_degup,ven_medown,ven_meup)
ggvenn(ven,c(1,2,3,4),text_size = 2.5)+
  scale_fill_igv()


####################Dif peak gene filtration##################

jh_silk2 <- jh_silk[`Qvalue(JH_psg/control_psg)` < 0.05,]
silk_gene[id %in% depeak$geneID,] %>% .[id %in% jh_silk2$GeneID,]


#########number of m6A peaks on m6A-containing genes#############

C1 <- fread("C1_IPVSC1_input_peak.csv") %>%
  .[pvalue <= 0.01 & log2FoldChange > 2,] %>%
  .[,geneID] %>% table() %>% as.data.frame() %>% setDT() %>% 
  .[,Freq] %>% table() %>% as.data.frame()
  
C2 <- fread("JH1_IPVSJH1_input_peak.csv") %>%
  .[pvalue <= 0.01 & log2FoldChange > 2,] %>% .[,geneID] %>% table() %>% as.data.frame() %>% setDT() %>% 
  .[,Freq] %>% table() %>% as.data.frame()

dt <- rbind(C1,C2) %>% setDT() %>% setnames(.,".","num") %>% .[!num %in% c(21,15),] %>% 
  .[1:12,group := "Con"] %>% .[13:22,group := "JH"] %>% .[]
dt %>% write.csv(.,file = "m6a_num.csv",row.names = FALSE)
fread("m6a_num.csv") %>% ggbarplot(.,x = "num",y = "Freq",
                 fill = "group",palette = "Paired",
                 position = position_dodge(),
                 label = TRUE,lab.pos = "out",alpha = .8)+
  theme_bw()+
  scale_fill_manual(values = c("grey","steelblue"))+
  xlab("")+ylab("")+
  theme(legend.position = "none")
fread("m6a_num.csv") %>% 
ggbarplot(.,x = "group", y = "Freq",
                                   position = position_fill(),
                                   fill = "num",rotate = TRUE,
                                   width = .7)+
  scale_fill_aaas()+
  theme_bw()+
  xlab("")+ylab("")

#############CIRCOS###########
library(circlize)
n <- 1000
df <- data.frame(
  sectors = sample(letters[1:8], n, replace = TRUE),
  x = rnorm(n), y = runif(n)
)

circos.par("track.height" = 0.1)
circos.initialize(df$sectors, x = df$x)

circos.track(df$sectors, y = df$y,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(7), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })

bgcol <- rep(c("#fb8072", "#80b1d3"), 4)
circos.trackHist(df$sectors, df$x, bin.size = 0.2, bg.col = bgcol, col = NA)

kyo <- fread("bombyxmori_karyotype.txt") %>% .[,Chr := gsub("omosome=","",Chr)]

circos.initializeWithIdeogram(kyo,ideogram.height = 5 )

kyo$End %>% class()

gene_density <- fread("Bom_density.csv") %>% unique(.,by = "GeneID")
circos.initializeWithIdeogram(kyo)
circos.genomicDensity(gene_density,col = c("steelblue"),track.height = .1,count_by = "number")
circos.clear()

  
m6A_density <- fread("C1_IPVSC1_input_peak.csv") %>% .[chr %like% "NC_",c(1,2,3,6,7)] %>% 
    .[`log2FoldChange` > 2&`pvalue` < 0.01,]
fun1 <- function(x){
  x %>% as.character() %>% gsub("NC_0513","",.) %>%
    as.integer() %>% -57 %>% as.integer() %>% paste0("chr",.)
}


m6A_density$chr <- fun1(m6A_density$chr)
m6A_density %>% setnames(.,c("Chr","Start","End","log2FoldChange","pvalue"))
circos.genomicDensity(m6A_density,col = c("#404948"),track.height = .1,count_by = "number")


DMG <- fread("JH1VSC1_diffPeak_result.csv") %>% .[chr %like% "NC_",c(1,2,3)]
DMG$chr <- fun1(DMG$chr)
setnames(DMG,c("Chr","Start","End"))
circos.genomicDensity(DMG,col = c("#D6AF6C"),track.height = .1,count_by = "number")

fread("JH_20E_treat.csv")

circos.genomicHeatmap(bed, col = col_fun, side = "outside",
                      line_col = as.numeric(factor(bed[[1]])))
bom_anno

bom_anno <- fread("Bom_annotation.csv")
DEG <- fread("jh_deg.csv") %>% na.omit() %>% .[`Qvalue (JH_psg / control_psg)` < 0.05,] %>% 
  bom_anno[.,on = .(GeneID = `Gene ID`)] %>% unique(.,by = "GeneID")
setnames(DEG,"log2 (JH_psg / control_psg)","log2FC")
DEG <- DEG %>% na.omit() %>% .[Chr %like% "chr",c(1:3,5)]
col_fun = colorRamp2(c(-3,0,3), c("blue", "white","red"))
circos.genomicHeatmap(DEG, col = col_fun, side = "inside",
                      connection_height = NULL,heatmap_height = .1 )
circos.clear()

df[,value:= seq(-3,3,by = 1)]
df <- seq(-3,3,by = 1) %>% as.data.frame() %>% setDT()
df <- df %>% .[1:3,v2:= "n"] %>% .[4:6,v2:="m"]

ggplot(df,aes(x = V1))+
  geom_point(aes(x = `.`,y = v2,fill = `.`))+
  scale_fill_gradient2(low="blue", high="red", mid="white")

#########Mutaion_qpcr result#######
library(ggsignif)
compar <- list(c("OE-Wt","JH+OE-Wt"))
fread("rip_qPCR_wt.csv") %>% 
  ggbarplot(.,x = "treat",y = "value",facet.by = "gene",add = "mean_sd",fill = "m6a")+
  geom_jitter(position = "dodge")+
  theme_bw()+
  stat_compare_means(comparisons = compar,method = "t.test")+
  xlab("")+ylab("Relative expression level")+
  theme(legend.position = "none")


silk_gene
seq_id <- c("100379325","101736771","692598","692688","692644","101742826")

write.table(.,file = "gene.txt",quote = FALSE,sep ="\n",row.names = FALSE,col.names = FALSE)

AnnotationDbi::select(orgdb,
                      keys = seq_id,
                      columns = c("ENTREZID","REFSEQ"),
                      keytype = "ENTREZID") %>% 
                      setDT() %>% .[REFSEQ %like% "P",2] %>%
                      write.table(.,file = "gene.txt",quote = FALSE,sep ="\n",row.names = FALSE,col.names = FALSE)

AnnotationDbi::select(orgdb,
      keys = seq_id,
      columns = c("ENTREZID","REFSEQ"),
      keytype = "ENTREZID") %>%setDT() %>%  .[REFSEQ %like% "M",] %>% unique(.,by = "ENTREZID") %>% 
  write.csv(.,file = "gene_list.csv",row.names=FALSE,col.names= FALSE)

##########qPCR#######

fread("qpcr_jh.csv") %>% t() %>% as.data.frame() %>% 
  setnames(.,c("seroin","Loc","BmSPI4","BmSPI5")) %>% .[-1,] %>%row.names()

qpcr <- fread("qpcr_jh_2.csv") %>% t() %>% as.data.frame() %>% setDT() %>% 
  setnames(.,c("seroin1","LOC101736771","BmSPI4","BmSPI5","Ldb","sercin2","YTHDF3","METTL3")) %>% .[-1,] %>%
  .[,group := c("Con","Con","Con","2h","2h","2h","4h","4h","4h","6h","6h","6h","8h","8h","8h","12h","12h","12h","24h","24h","24h")] %>% 
  .[] %>% melt(.,value.name = "value",id.vars = "group",measure.vars = .SD,variable.name = "gene" )

qpcr$value <- as.numeric(qpcr$value)

ggbarplot(qpcr,x = "group",y = "value",facet.by = "gene",add = "mean_sd",fill = "steelblue",alpha = .8)+
  geom_jitter(position = "dodge",size = 1,alpha = .8)+
  theme_bw()+
  stat_compare_means(comparisons = compar,method = "t.test",size = 2)+
  xlab("")+ylab("Relative expression level")
compar <- list(c("Con","2h"),c("Con","4h"),c("Con","6h"),c("Con","8h"),c("Con","12h"),c("Con","24h"))


#################thickness and weight of silkworm cocoon##############

library(ggplot2)
library(ggpubr)
library(data.table)
compar = list(c("male_Con","male_JH"),c("female_Con","female_JH"))
cocon <- fread("cocon_weight.csv") 
cocon %<>% .[,treat_2 := paste0(cocon$gender,"_",cocon$treat)] %>% .[]

ggboxplot(cocon, x = "treat_2",y = "weight",fill = "treat",order = c("male_Con","male_JH","female_Con","female_JH"))+
  geom_point(position = "dodge",size = 1, alpha = .8,color = "black")+
  expand_limits(y = 0)+
  theme_bw()+
  stat_compare_means(comparisons = compar,method = "t.test")+
  xlab("")+ylab("Cocon weight")+
  theme(legend.position = "none")+
  theme_bw()


compar_2 = list(c("Con_male","JH_male"),c("Con_female","JH_female"))
cocoon_thick <- fread("cocoon_thickness.csv")
ggboxplot(cocoon_thick, x = "group",y = "thickness",fill = "treat",order = c("Con_male","JH_male","Con_female","JH_female"))+
  geom_point(position = "dodge",size = 1, alpha = .8,color = "black")+
  expand_limits(y = 0)+
  theme_bw()+
  stat_compare_means(comparisons = compar_2,method = "t.test")+
  xlab("")+ylab("Cocoon thickness")+
  theme(legend.position = "none")
theme_bw()


##############cojoint analysis of m6A-seq and RNAi rna-seq data###########
fread("m3_rnai.csv") %>% .[Gene_ID %in% silk_gene$id,] %>% 
  .[`Qvalue_(siM3_/_NC)` < 0.05,] %>% silk_gene[,c(1,3)][., on = .(id =Gene_ID)] %>% .[,]
  melt(.,id.vars = "name", measure.vars = .SD, variable.name = "group",value.name= "FPKM") %>%
  .[group %like% "NC",treat := "Con"] %>% 
  .[group %like% "si",treat := "siM3"] %>% .[] %>% 
  ggbarplot(.,x = "treat", y = "FPKM",add = "mean_sd",facet.by = "name",scales = "free_y",fill = "steelblue")+
  geom_point(color = "black",alpha = .7,position = "dodge")+
  theme_bw()

library(data.table)
fread("") %>% 
  .[`Qvalue_(siM3_/_NC)` < 0.05 ] %>% .[Gene_ID %in% m6A_containing_genes$geneID,]

m6A_containing_genes <- fread('C1_IPVSC1_input_peak.csv') %>% .[padj < 0.01 & log2FoldChange > 1] %>% unique(.,by = "geneID")

silk_gene %>% .[id %in% m6A_containing_genes$geneID,]

m6a_containing_rnai <- fread("m3_rnai.csv") %>% .[Gene_ID %in% m6A_containing_genes$geneID & `Qvalue_(siM3_/_NC)` < 0.05,]

enrichGO(gene = m6a_containing_rnai$Gene_ID,
         OrgDb = orgdb,
         keyType = "ENTREZID",
         ont = "ALL",
         pvalueCutoff = 1,
         qvalueCutoff = 1) %>% 
  dotplot(showCategory = 30)


###########thickness and weight of cocoon#############
cocoon_thick[treat == "JH",thickness] %>% mean() / cocoon_thick[treat == "Con",thickness] %>% mean()

cocon[treat == "JH",weight] %>% mean() / cocon[treat == "Con",weight] %>% mean()



