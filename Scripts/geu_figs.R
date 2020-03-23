library(data.table)
#library(parallel)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
#library(reshape2)
#library(xtable)
#library(rms)
#library(ggbio)
#library(Homo.sapiens)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library("org.Hs.eg.db")
##library(EnsDb.Hsapiens.v75)
library(biovizBase)
library(GenomicAlignments)
library(GenomicFeatures)
library(RColorBrewer)
##library(VennDiagram)
library(ggrepel)
library(hrbrthemes)
library(lemon)
library(ggforce)


##detach("package:org.Hs.eg.db", unload=TRUE)

#source("/home/ev250/Bayesian_inf/trecase/Functions/Btrecase.noGT.R")
source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")
source('/home/ev250/Cincinatti/Functions/various.R')
source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R")
source("/home/ev250/Bayesian_inf/trecase/Functions/various.R")


###############################################################################################################
########## Open files ########################################################################################
#############################################################################################################


## Look at genes run with GT with ref bias correction

## get wildcards to identify files 

rbias <- snakemake@params[['rbias']]
## rbias=c("norefbias","rbias")

## btrecase_dir  <- paste0("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/", c("GT","fisher001/RNA"))
btrecase_dir <- snakemake@params[['btrecase']]
## get start end fo genes
## gene.coord <- fread("/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt")

gene.coord <- fread(snakemake@inputs[['geneStEnd']]) 

## select genes in chrom 22

gt22 <- gene.coord[chrom==22,]

## add tag distance to gene (closest to start or end)
d <- 10^5
maxRhat <- 1.1
max.info <- 1.1
## For PEP col add 1 to the  4000 posterior draws so I can take log
btrecase <- lapply(btrecase_dir, function(i) {
    tmp <- lapply(rbias, function(j) {
        tmp <- comb.files(path=i, pattern=paste0("^",j, "\\.ENSG[0-9]+.*stan.summary.txt"))
        ## remove bad runs as recommended in stan
        tmp <- tmp[Rhat < maxRhat,]
        tmp <- gene.d(tmp, gt22[, .(gene_id, start,end,chrom)])
        if(any(names(tmp)=="info")) tmp <- tmp[info <= max.info,]
        if(any(names(tmp)=="PEP")) tmp[PEP==0, PEP:=1/4001][, minuslog10PEP:=-log(PEP,10)]
        ## add 95% CI based null column
        tmp <- add.null(tmp)
                   })
    names(tmp) <- rbias
    return(rbindlist(tmp, idcol="rbias", fill=T))})

names(btrecase) <- basename(btrecase_dir)

## save
saveRDS(btrecase, "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/btrecase.GT.nGT.all.rds")

## GEt number of genes run for GT/RNA
lapply(btrecase, function(i) i[, unique(Gene_id), rbias][,.N,rbias])


## get number of assoc run by gene.dist
lapply(btrecase, function(i) i[rbias=="rbias" , .N, .(rbias)])
lapply(btrecase, function(i) i[rbias=="rbias"  & abs(gene.dist) <= d, .N, .(rbias)])

## get total significant with rbias correction
lapply(btrecase, function(i) i[rbias=="rbias" & null.99 == "no", unique(Gene_id), .(rbias)])
lapply(btrecase, function(i) i[model != "trec" & rbias=="rbias" & null.99 == "no", unique(Gene_id), .( model,rbias)])

## get sig associations with rbias by gene.dist
lapply(btrecase, function(i) i[rbias=="rbias" & null.99 == "no", .N, .(rbias)])
lapply(btrecase, function(i) i[rbias=="rbias" & null.99 == "no" & abs(gene.dist) <= d, .N, .(rbias)])

## total with d
lapply(btrecase, function(i) i[abs(gene.dist) <= d, unique(Gene_id), rbias][,.N,rbias])

## Get tags run with GT to combine with RNA (same tags run with or without ref panel bias correction, use only one)

gt.tags <- comb.files(path=btrecase_dir[1], pattern=paste0("^",rbias[2], "\\.ENSG[0-9]+.*tags.lookup.txt"))
gt.tags.run <- merge(gt.tags, btrecase[[1]][rbias=="rbias",.(Gene_id,tag)], by= c("Gene_id","tag"))

##############################
###  out dir ################
############################

## out.dir <- '/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures'
out.dir <- snakemake@params[['out']]

############################################################################################################################
######## Compare gt vs rna within 100kB window #############################################################################
###########################################################################################################################


## merge gt and rna considering different tags, effect size in RNA already corrected by op.dir, only use trec-ase to compare same model

gt.rna <- comp.ngt(path.s=btrecase_dir[2], 
                   pattern.s = paste0("^",rbias[2], "\\.ENSG[0-9]+.*stan.summary.txt"),
                   pattern.t= paste0("^",rbias[2], "\\.ENSG[0-9]+.*tags.lookup.txt"), #"eqtl.tags.lookup.txt",
                   lm.sum=btrecase[["GT"]][rbias == "rbias" & model=="trec-ase", ] ,
                   gt.sum= btrecase[["GT"]][rbias == "rbias" & model == "trec-ase", ],
                   s=c(".rem", ".gt", ".rna"),
                   tags.m2=gt.tags.run)


## effect size in rna is based on gt tag maf (op.dir column)

## remove ".rem" col, artifact to run conp.ngt function

gt.rna[ , grep("rem$", names(gt.rna), value=T):=NULL]

## remove associations high info from gt.rna

gt.rna <- gt.rna[info.rna<=max.info,]
gt.rna <- gt.rna[Rhat.rna<maxRhat & Rhat.gt<maxRhat,]

## add null.95.rna col to gt.rna

gt.rna <- add.null(gt.rna, suffix=".rna")

## make -log10PEP col, if PEP==0, add 1 count to posterior draw of 4000 samples

gt.rna[PEP.rna==0, PEP.rna:=1/4001][, minuslog10PEP.rna:=-log(PEP.rna,10)]

## Look for inconsistencies due to changes in eaf being too close to 0.5, effect size in RNA already corrected by op.dir 

gt.rna[(null.95.gt=='no' & null.95.rna == 'no' & sign(log2_aFC_mean.gt) != sign(log2_aFC_mean.rna) ), .(tag.gt, tag.EAF.gt, tag.rna, tag.EAF.rna, Gene_id)]


## correct if necessary, mean and CIs

gt.rna[(null.95.gt=='no' & null.95.rna == 'no' & sign(log2_aFC_mean.gt) != sign(log2_aFC_mean.rna)), log2_aFC_mean.rna:= -log2_aFC_mean.rna]
gt.rna[(null.95.gt=='no' & null.95.rna == 'no' & sign(log2_aFC_mean.gt) != sign(log2_aFC_mean.rna)),
       grep("log2_aFC_[0-9].*rna", names(gt.rna), value=T):= lapply(grep("log2_aFC_[0-9].*rna", names(gt.rna), value=T), function(i) -get(i))]


## save gt.rna
saveRDS(gt.rna, "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/btrecase.GT.ngt.wide.rds")


## Select association within 100kB, null 99  and plot 

sig <- 99



cols <- c(None="#999999", `obs-GT`="yellow3", `hidden-GT`="#0072B2", Both= "#D55E00")



rna.gt.Sig <-lapply(sig, function(i) {
    
    gt.rna<- add.signif(gt.rna, x1=paste0("null.", i, ".gt"), x2=paste0("null.",i, ".rna"), col=c("obs-GT","hidden-GT"))

    gt.rna.tab <- tab2bplot(dt=gt.rna[ abs(gene.dist.gt) <= d  ,], colors=cols)

    rna.gt <- btrecase.plot(dt=gt.rna[abs(gene.dist.gt) <= d & Signif != "None",],
                        x1=c('log2_aFC_mean.gt', 'log2_aFC_0.5%.gt','log2_aFC_99.5%.gt', paste0('null.',i,'.gt')),
                        x2= paste0(c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', paste0('null.', i)), ".rna"),
                        xl='eQTL-effect (observed GT)', yl='eQTL-effect (hidden GT)',
                        col=c("obs-GT","hidden-GT"),axis.title=12, axis.text=10,
                        legend.title=12, legend.text=10,
                        legend.symbol=4, point.size=3 ,
                        title="Observed vs hidden genotypes", title.size=12) +
    annotation_custom(tableGrob(gt.rna.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(gt.rna.tab$color, rep("black", 4)))))), xmin=-1.1, xmax=-0.3, ymin=0.35, ymax=0.45)

    print(cor.test(gt.rna[abs(gene.dist.gt) <= d & Signif != "None" ,log2_aFC_mean.gt],
              gt.rna[abs(gene.dist.gt) <= d & Signif != "None"  ,log2_aFC_mean.rna]))
    return(rna.gt)
    })


##################################################################################################
## Compare gt trec-ase and rna in terms of associations tested and Signif by gene distance
###################################################################################################

## by info:
common.sig.a <- sapply(seq(0.3,1, 0.1), function(i) gt.rna[abs(gene.dist.gt) <= d & Signif == "Both" & op.dir=="no" & info.rna>=i,.N])

total.sig.a.ngt <- sapply(seq(0.3,1, 0.1), function(i) gt.rna[abs(gene.dist.gt) <= d & null.99.rna=="no" & op.dir=="no" & info.rna>=i,.N])

gt.rna.assoc.info <- data.table(info=seq(0.3,1, 0.1), common=common.sig.a, total=total.sig.a.ngt )
gt.rna.assoc.info[, p:=common/total]

## get trec-ase entries for  gt100 to format long to look at SNP distance to gene, all tested and signif

## choose null

level <- 99
null <- paste0("null.",level)


## all associations
gt.rna.long <- rbindlist(lapply(btrecase, function(i) i[rbias == "rbias" ,]), idcol="source", fill=T)


gt.rna.long <- gt.rna.long[abs(gene.dist) <= d,]

## Add Signif associations only for GT and RNA to make facet plot, indicate that in source col, use 99%CI 

gt.rna.sig <- rbindlist(lapply(unique(gt.rna.long$source), function(i) {
    temp <- gt.rna.long[source == i & get(null) =="no",]
    temp[, source:=paste0(source, ".Sig")]
    return(temp)
    }))


gt.rna.long <- rbind(gt.rna.long,gt.rna.sig)

## recode to have facet by Sig/All overlapping density for GT and RNA

gt.rna.long[, type:="All"][grep("Sig", source), type:="Signif"]
gt.rna.long[grep("GT", source), source:="obs-GT"][grep("RNA", source), source:="hidden-GT"]

## prepare tables to add to plot

table2d <- lapply(unique(gt.rna.long$type), function(i) {
    dt <- gt.rna.long[type==i,.N,source] #[order(1/N),]
    names(dt) <- c("", "SNPs")
    return(dt)
    
})


## change table2d[[2]][2,2] to make it compatible with gt.rna (tags issue)

#table2d[[2]][2,2] <- nrow(gt.rna[ Signif =="Both" | Signif == "hidden-GT" ,])



gl <- lapply(table2d, function(i) tableGrob(i, rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c("yellow3", "#0072B2", rep("black", 2)))))))


dt <- data.table(type=unique(gt.rna.long$type), grob = gl )

cols.d <- c("#0072B2","#F0E442")

dist.gt.ngt <- ggplot(gt.rna.long[abs(gene.dist) <=d,], aes(abs(gene.dist)/1000, color=source, fill=source))  +      
    theme_bw() +
    xlab("Distance to gene (kB)") +
    ylab("Density") +
    theme(legend.title = element_blank()) +     
    guides(colour=guide_legend(override.aes = list(size = 3))) +
    geom_density(alpha = .1) +
    scale_color_manual(values=cols.d) +
    scale_fill_manual(values=cols.d) +
    facet_wrap(~type) +
    theme(strip.text.x = element_text(size=12),
          strip.background=element_rect(fill="white")) + 
    geom_custom(data=dt, aes(grob=grob), x = 0.6, y = 0.9)

dist.gt.ngt.hist <- ggplot(gt.rna.long[abs(gene.dist) <=d,], aes(abs(gene.dist)/1000, color=source, fill=source))  +      
    theme_bw() +
    xlab("Distance to gene (kB)") +
    ylab("Frequency") +
    theme(legend.title = element_blank()) +     
    guides(colour=guide_legend(override.aes = list(size = 3))) +
    geom_histogram(alpha = .1, position="dodge") +
    scale_color_manual(values=cols.d) +
    scale_fill_manual(values=cols.d) +
    facet_wrap(~type, scales="free") +
    theme(strip.text.x = element_text(size=12),
          strip.background=element_rect(fill="white")) + 
    geom_custom(data=dt, aes(grob=grob), x = 0.6, y = 0.9)



################################################################################
## Compare Obs-GT and Hidden-Gt by genes
################################################################################


u.sig <- lapply(unique(gt.rna.long$source), function(i) unique(gt.rna.long[type== "Signif"& null.99=="no" & source ==i & abs(gene.dist) <= d, Gene_id]))

## Get the number of common genes tested in both datasets

common.sig <- lapply(u.sig, function(i) i[i %in% unique(gt.rna$Gene_id)])
names(common.sig) <- unique(gt.rna.long$source)

gt.only <- Reduce(setdiff, common.sig)
rna.only <- Reduce(setdiff, rev(common.sig))
common <- Reduce(intersect, common.sig)

## Distribution of significant SNPs by gene
lapply(unique(gt.rna.long$source), function(i) gt.rna.long[type== "Signif" & null.99=="no" & source ==i & abs(gene.dist) <= d, .N,Gene_id])


df.venn <- data.table(x=c(-1,-1), y=rep(0,2), labels=c("Obs-GT", "Hidden-GT"))

v <- ggplot(df.venn) +
    geom_circle(aes(x0=x,y0=y,r=c(2, 1), fill=labels, color=labels),  alpha=0.25, size=.5) +
    geom_rect(mapping=aes(xmin=-3.5,xmax=2, ymin=-2.7, ymax=2.7, fill="Total genes"), alpha=0.01, color="black") +
    coord_fixed() +
    theme_void() +
    scale_fill_manual(values = c(cols.d, "white") ) +
    scale_colour_manual(values = c(cols.d, "black"), guide = FALSE) +         
    labs(fill = NULL) +
    annotate("text", x= c(-1, -1, -1), y=c(1.3,0, 2.4), label=c(unlist(lapply(list(gt.only, common),length)),length(unique(gt.rna$Gene_id))), size=5, colour=c("orange2", "#0072B2", "black")) +
    
    ggtitle("Overlap of genes with significant eQTLs") +
    theme(plot.title = element_text(hjust = 0.5))


###################################################################################
## Look at individual genes with gt and rna only ##
###################################################################################


## work with 100,000 kb

## to use previous functions need to add ".ngt" data, I will just duplicate rna for simplicity

rna.ngt <- merge(btrecase$RNA[rbias=="rbias",], btrecase$RNA[rbias=="rbias",], by=c("Gene_id", "tag", "tag.EAF", "gene.dist", "rbias"), suffixes=c(".ngt", ".rna"))

rna.ngt100 <- rna.ngt[abs(gene.dist)<=d,]

gt.rna.100 <- gt.rna[abs(gene.dist.gt) <= d,]

## add gene name

gt.rna.100 <- add.name(gt.rna.100)


## add ngt cols to gt.rna.100

all.100 <- merge(gt.rna.100, rna.ngt[,c("Gene_id", "tag", grep(".ngt$", names(rna.ngt), value=T)), with=F], by.x=c("Gene_id", "tag.rna"), by.y=c("Gene_id", "tag"))

## add col tag.ngt

all.100[ , tag.ngt:=tag.rna]


## select MAPK1 as example and add gene track

## gene="NAGA"
gene="MAPK1"
gene_id =unique(gt.rna.100[Gene_name == gene, Gene_id])

data(genesymbol, package = "biovizBase")

cis.w <- unlist(gt22[gene_id == unique(gt.rna.100[Gene_name == gene, Gene_id]), .(start, end)]) + c(-d -10, d+10)

wh <- range(genesymbol[seqnames(genesymbol) == "chr22" & start(genesymbol)  > cis.w[1] & end(genesymbol) < cis.w[2]],
            ignore.strand = T)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

gene_symbol<- data.table(AnnotationDbi::select(org.Hs.eg.db, keys=genes(txdb)$gene_id, 

                                    columns="SYMBOL", keytype="ENTREZID"))

temp <- crunch(txdb, which = wh)#genesymbol[c("NAGA","NDUFA6", "SEPT3" )]) #, type = "reduce")

## remove missing gene_id

temp <- temp[temp$gene_id != ""]

entrez <- data.table(ENTREZID=temp$gene_id)
entrez <- merge(entrez, gene_symbol, by = "ENTREZID", all.x = T, sort=F)
              
temp$symbol = entrez$SYMBOL

tmp <- split(temp, temp$symbol)

## get start for first exon to place legend
all.st <- as.data.table(temp)
all.st <- all.st[type=="exon",][, .SD[1], symbol]
all.st[, col:="springgreen4"][symbol==gene, col:="black"]
setkey(all.st, start)
all.st[, y:= c(1.4, 2.4,1.4, 1.45)]
all.st[, hjust:= c(0,0,0,0)]
##all.st[strand=="+", y:= c(0.6,1.5)][strand=="-", y:=c(2.5,1.8)]


## window for rectangle
st.end <- unlist(gt22[gene_id == unique(gt.rna.100[Gene_name == gene, Gene_id]), .(start, end)])

w <- 2500

## need to add gene names manually as otherwise they overlap or arent the right size

p.txdb <- autoplot(tmp, aes(type=type), label=F, label.color = "black", color = "springgreen4",fill = "springgreen4")+ xlim(cis.w) + scale_x_sequnit("Mb") + geom_rect(mapping=aes(xmin= st.end[1]-w, xmax=st.end[2] + w, ymin=0.7, ymax=1.3), fill=NA, colour="gray") +
    theme(##legend.position = "none",
        ##panel.grid = element_blank(),
        ##axis.title = element_blank(),
        axis.text = element_text(size=10) ##,
        ##axis.ticks.x = element_blank())
        ) +
    
        annotate("text", label = all.st[["symbol"]], x=all.st[["start"]], y=all.st[["y"]], colour=all.st[["col"]], hjust = all.st[['hjust']], size=2)

## Add top cis-SNP to plots

## snp22 <- "/mrc-bsu/scratch/ev250/reference_genome/built37/variation_homo_sapiens-chr22.gvf.gz"

snp22 <- snakemake@input[['var22']]


## get PEP to add to mapk1

mapk.pep <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/PEP/GT/ENSG00000100030.stan.summary.txt")
mapk.pep <- merge(btrecase$GT[rbias == "rbias" & abs(gene.dist) <= d,], mapk.pep[,.(Gene_id,tag,PEP)], by=c("Gene_id", "tag"))
## adding 1 to PEP if 0, then take -log10
mapk.pep[PEP==0, PEP:=1/4001][, minuslog10PEP:=-log(PEP,10)]



mapk1 <- gene.plot2b(x=mapk.pep,   #btrecase$GT[rbias == "rbias" & abs(gene.dist) <= d,],
                     ci="no",
                     ci.x=level/100,
                     null.x=null,
                     y=rna.ngt100,
                     ci.y=level/100,
                     null.y=null,
                     ci.z=level/100,
                     yvar="minuslog10PEP",
                     yaxis=bquote(-log[10](PEP)),
                     xcoord=0.05,
                     ycoord=0.9,
                     z=all.100,
                     gene= gene_id,
                     fisher=NULL,
                     info.s=NULL,
                     gene.track=p.txdb@ggplot,
                     rsid=snp22,
                     hline=-log(0.005, 10))

## p <- lapply(unique(btrecase$RNA[["Gene_id"]]) , function(i)  gene.plot2b(x=btrecase$GT[rbias == "rbias" & abs(gene.dist) <= d,],
##                    ci.x=level/100,
##                    null.x=null,
##                    y=rna.ngt100,
##                    ci.y=level/100,
##                    null.y=null,
##                    ci.z=level/100,
##                    z=all.100,
##                    gene=i,
##                    gene.track=NULL,
##                    rsid=NULL) )


######################################################################################
## Make Figure

top <- plot_grid(dist.gt.ngt, ncol=2, rna.gt.Sig[[1]],  labels='auto')
top2 <- plot_grid(dist.gt.ngt.hist, ncol=2, rna.gt.Sig[[1]],  labels='auto')

bot <- plot_grid(v, mapk1, ncol=2, labels=c("c","d"))

p <- plot_grid(top, NULL, bot, ncol=1, rel_heights=c(0.8,0.1,1))
p2 <- plot_grid(top2, NULL, bot, ncol=1, rel_heights=c(0.8,0.1,1))

## get width and height using par("din") after selecting a windonw that looks good
ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/Fig4.png', p, width = 11.28, height = 10.1)

ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/Fig4alt.png', p2, width = 11.28, height = 10.1)

ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/obs_hidden.png',rna.gt.Sig[[1]], width=5.5, height=3.7)

ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/Venn_obs_hidden.png',v, width=5.5, height=3.7)

####################################################################################################################
## External validity/calibration  figure 
####################################################################################################################

## Use Gtex-ebv as gold standard and compare Deseq, Rasqual and Btrecae-GT

######################################################################################################################################
##' # Use Gtex-ebv as gold standard and compare to Trec, 
##' # DEseq2, Rasqual and Trecase. Choose FDR in Dseq that gives a similar number of
##' # associations as Gtex-ebv
########################################################################################################################################

##' Open and format datasets: Gtex-ebv, DEseq and Rasqual

########################################
## Look at Gtex ebv for chromosome 22 ##
########################################

## ebv <- fread(cmd="zcat /mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.allpairs.txt.gz | awk -F'\t' 'NR == 1 || $2~/^22_/' ", header=T, sep="\t")
ebv <- fread(cmd=paste("zcat", snakemake@input[['gtex']], "| awk -F'\t' 'NR==1 || $2~/^22_/' "))

## ebv.sig <- fread(cmd="zcat /mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.v7.signif_variant_gene_pairs.txt.gz")

ebv.sig <- fread(cmd=paste("zcat", snakemake@input[['sigGtex']]))

ebv[, Gene_id:=gsub("\\..*","",gene_id)]
ebv <- ebv[Gene_id %in% unique(btrecase$GT$Gene_id),]
ebv[, SNP:=gsub("^22_|_b37$", "", variant_id)][, SNP:=gsub("_", ":",SNP)]

ebv.sig <- ebv.sig[variant_id %in% ebv$variant_id,][, null:="no"]
ebv <- merge(ebv,ebv.sig[,.(gene_id, variant_id, null)], by=c("gene_id", "variant_id"), all.x=T)
ebv[is.na(null), null:="yes"]

## add FDR by R for comparison

ebv[,p.adj:= p.adjust(pval_nominal,method = "BH")]
ebv[,null.fdr5:="yes"][p.adj<=0.05, null.fdr5:="no"]


######################################
## FDR range for frequentist tests ###
######################################

r.fdr <- sort(c(10^(seq(-5,-1,1)), 0.05, 0.5, 0.8, 1))

###########
## DEseq ##
##########
## Open and format DEseq2 output (run by Wei-Yu)

## dseq <- rbindlist(lapply(list.files("/mrc-bsu/scratch/wyl37/ElenaData/RunNBmodel", pattern="RunNBmodelbatch[0-9]+_chr22.nbmodelRes.csv", full.names=T), fread))

dseq <- rbindlist(lapply(snakemake@input[['dseq']], fread))

dseq[, SNP:=NULL]

## remove pval NA

dseq <- dseq[!is.na(pvalue),]

## add log2_aFC to dseq

dseq[, log2_aFC:=log2FoldChange*2]

## add BH correction to dseq

setkey(dseq, pvalue)
dseq[,p.adj:= p.adjust(pvalue,method = "BH")]
setkey(dseq, p.adj)

## Add null column for dseq based on 5% FDR
dseq[,null.fdr5:="yes"][p.adj<=0.05, null.fdr5:="no"]

## Add gene distance

dseq <- gene.d(dseq, gt22[, .(gene_id, start,end,chrom)], snp="tag")



#############
## Rasqual ##
#############

## rasq1 <- list.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/outCounts/cis5_10_5/", "ENSG[0-9]+.*txt", full.names=T)
## rasq.header <- list.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/output", "rasqual.header.txt", full.names=T)
## running rasqual small cis-window, only 1 gene excluded 
## rasq2 <- list.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/outCounts/cis1_10_5", "ENSG[0-9]+.*txt", full.names=T)

rasq1 <- list.files(snakemake@params[['rasqual']], "ENSG[0-9]+.*txt", full.names=T)
rasq2 <- list.files(snakemake@params[['rasqual01']], "ENSG[0-9]+.*txt", full.names=T)


rasq.header <- paste0(snakemake@params[['rasqual']], "/rasqual.header.txt")

rasq <- function(x){ ## format rasqual input
    rasqual <- rbindlist(lapply(x, format_rasqual, top.hits="no", header=rasq.header))

    ## rasqual didnt run some genes "Estimated computational time is too long...", remove from output
    nrow(rasqual[rs_id == "SKIPPED",])
    rasqual <- rasqual[rs_id != "SKIPPED",]
    rasqual <-  merge(dseq[,.(Gene_id, tag)], rasqual, by.x=c("Gene_id", "tag"), by.y=c("gene_id", "rs_id"))
    rasqual[, p_adjust:= p.adjust(p, method="BH")]
    rasqual[, log2_aFC:=log2(Fold_change)]
    for (i in r.fdr){
        rasqual[,eval(paste0("null.fdr",i*100)):= "yes"][p_adjust<=i, paste0("null.fdr",i*100):="no"]
    }
    return(rasqual)
}

rasqual <- rasq(rasq1)    

rasqual2 <- rasq(rasq2)

#########
## LM
########
### lm.f <-  list.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/lm/log_counts",  "ENSG[0-9]+.*summary.txt", full.names=T)
lm.f <- list.files(snakemake@params[['lm']], "ENSG[0-9]+.*summary.txt", full.names=T)
lm <- rbindlist(lapply(lm.f, fread))
lm[, p_adjust:= p.adjust(log2.pvalue, method="BH")]

for (i in r.fdr){
    lm[,eval(paste0("null.fdr",i*100)):= "yes"][p_adjust<=i, paste0("null.fdr",i*100):="no"]
    
}

#############################################################
## Select same associations in rasqual, btrecase-GT, dseq and lm

## 5*10^5 cis-window
assoc <- Reduce(fintersect, list(rasqual[,.(Gene_id,tag)], dseq[,.(Gene_id, tag)], btrecase$GT[rbias=="rbias" , .(Gene_id,tag)], lm[,.(Gene_id, tag)]))

## 10^5 cis-window
assoc2 <- Reduce(fintersect, list(rasqual2[,.(Gene_id,tag)], dseq[,.(Gene_id, tag)], btrecase$GT[rbias=="rbias" , .(Gene_id,tag)], lm[,.(Gene_id, tag)]))

common.all <- lapply(list(Rasqual=rasqual, Btrecase=btrecase$GT[rbias=="rbias",], Deseq=dseq, lm=lm), function(i) merge(i, assoc, by=names(assoc)) )

common.all2 <- lapply(list(Rasqual=rasqual2, Btrecase=btrecase$GT[rbias=="rbias",], Deseq=dseq, lm=lm), function(i) merge(i, assoc2, by=names(assoc2)) )

## Convert rasqual, dseq2 and lm to long format by FDR and btrecase using null for common associations only to ease comparison with Gtex-ebv

rasq.l <- reshape(common.all$Rasqual[,c("Gene_id", "tag","log2_aFC", grep("null.fdr",names(common.all$Rasqual), value=T)), with=F],
                  direction="long",
                  varying=list( grep("null.fdr",names(common.all$Rasqual), value=T)),
                  v.names="null.Fdr",
                  times= as.numeric(gsub("null.fdr", "", grep("null.fdr",names(common.all$Rasqual), value=T)))/100,
                  timevar="Fdr")

rasq.l2 <- reshape(common.all2$Rasqual[,c("Gene_id", "tag","log2_aFC", grep("null.fdr",names(common.all2$Rasqual), value=T)), with=F],
                  direction="long",
                  varying=list( grep("null.fdr",names(common.all2$Rasqual), value=T)),
                  v.names="null.Fdr",
                  times= as.numeric(gsub("null.fdr", "", grep("null.fdr",names(common.all2$Rasqual), value=T)))/100,
                  timevar="Fdr")

dseq.l <- rbindlist(lapply(r.fdr, function(i) {
    tmp <- common.all$Deseq[, null.fdr:="yes"][p.adj<=i, null.fdr:="no"]
    tmp[, Fdr:= i]
    tmp <- tmp[, .(Gene_id, tag, log2FoldChange,log2_aFC, p.adj, null.fdr, Fdr, gene.dist)]
    return(tmp)
}))

dseq.l2 <- rbindlist(lapply(r.fdr, function(i) {
    tmp <- common.all2$Deseq[, null.fdr:="yes"][p.adj<=i, null.fdr:="no"]
    tmp[, Fdr:= i]
    tmp <- tmp[, .(Gene_id, tag, log2FoldChange,log2_aFC, p.adj, null.fdr, Fdr, gene.dist)]
    return(tmp)
}))


lm.l <- reshape(common.all$lm[,c("Gene_id", "tag","log2.aFC_mean", grep("null.fdr",names(common.all$lm), value=T)), with=F],
                  direction="long",
                  varying=list( grep("null.fdr",names(common.all$lm), value=T)),
                  v.names="null.Fdr",
                  times= as.numeric(gsub("null.fdr", "", grep("null.fdr",names(common.all$lm), value=T)))/100,
                timevar="Fdr")

lm.l2 <- reshape(common.all2$lm[,c("Gene_id", "tag","log2.aFC_mean", grep("null.fdr",names(common.all2$lm), value=T)), with=F],
                  direction="long",
                  varying=list( grep("null.fdr",names(common.all2$lm), value=T)),
                  v.names="null.Fdr",
                  times= as.numeric(gsub("null.fdr", "", grep("null.fdr",names(common.all2$lm), value=T)))/100,
                timevar="Fdr")


btrecase.gt.l <- rbindlist(lapply(c(95, 99), function(i) {
    dt <- rej.recode(a=0,b=common.all$Btrecase, c=i/100)
    dt[,PEP:=1-post.out]
    null <- paste0( "null.", i)
    tmp <- dt[, c("Gene_id", "tag", "log2_aFC_mean", "PEP", "null.rej", null), with=F]
    setnames(tmp, eval(null), "null")
    tmp[, PIP:= i/100]
    
}))

btrecase.gt.l2 <- rbindlist(lapply(c(95, 99), function(i) {
    dt <- rej.recode(a=0,b=common.all2$Btrecase, c=i/100)
    dt[,PEP:=1-post.out]
    null <- paste0( "null.", i)
    tmp <- dt[, c("Gene_id", "tag", "log2_aFC_mean", "PEP", "null.rej", null), with=F]
    setnames(tmp, eval(null), "null")
    tmp[, PIP:= i/100]
    
}))

## extend btrecase.gt.l to PIP=0.85 and 0.9, null.rej corresponds to normal approx., copy to null to merge with btrecase.gt.l

post.dt <- rbindlist(lapply(c(0.85, 0.9), function(i) {
    dt <- rej.recode(a=0,b=common.all$Btrecase, c=i)
    dt[,PEP:=1-post.out]
    dt[,PIP:=i]
    dt[, null:=null.rej]
    dt <- dt[,.(Gene_id, tag, log2_aFC_mean,PEP,null.rej,null, PIP)]
}))

post.dt2 <- rbindlist(lapply(c(0.85, 0.9), function(i) {
    dt <- rej.recode(a=0,b=common.all2$Btrecase, c=i)
    dt[,PEP:=1-post.out]
    dt[,PIP:=i]
    dt[, null:=null.rej]
    dt <- dt[,.(Gene_id, tag, log2_aFC_mean,PEP,null.rej,null, PIP)]
}))

btrecase.gt.l <- rbind(btrecase.gt.l, post.dt)   
btrecase.gt.l2 <- rbind(btrecase.gt.l2, post.dt2)   

### Merge with gtex-ebv

rasq.ebv <- merge(rasq.l, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
rasq.ebv <- add.signif(rasq.ebv, x1="null.Fdr", x2="null", col=c("Rasqual", "Gtex-ebv"))

dseq.ebv <- merge(dseq.l, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
dseq.ebv <-  add.signif(dseq.ebv, x1="null.fdr", x2="null", col=c("DEseq","Gtex-ebv") )

btrecase.gt.ebv <- merge(btrecase.gt.l, ebv,  by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"), suffixes=c(".btrecaseGT", ".ebv"))
btrecase.gt.ebv <- add.signif(btrecase.gt.ebv, x1="null.btrecaseGT", x2="null.ebv", col=c("Btrecase", "Gtex-ebv"))

rasq.ebv2 <- merge(rasq.l2, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
rasq.ebv2 <- add.signif(rasq.ebv2, x1="null.Fdr", x2="null", col=c("Rasqual", "Gtex-ebv"))

dseq.ebv2 <- merge(dseq.l2, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
dseq.ebv2 <-  add.signif(dseq.ebv2, x1="null.fdr", x2="null", col=c("DEseq","Gtex-ebv") )

btrecase.gt.ebv2 <- merge(btrecase.gt.l2, ebv,  by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"), suffixes=c(".btrecaseGT", ".ebv"))
btrecase.gt.ebv2 <- add.signif(btrecase.gt.ebv2, x1="null.btrecaseGT", x2="null.ebv", col=c("Btrecase", "Gtex-ebv"))

lm.ebv <- merge(lm.l, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
lm.ebv <- add.signif(lm.ebv, x1="null.Fdr", x2="null", col=c("lm", "Gtex-ebv"))

lm.ebv2 <- merge(lm.l2, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
lm.ebv2 <- add.signif(lm.ebv2, x1="null.Fdr", x2="null", col=c("lm", "Gtex-ebv"))

## calculate FDR for btrecase using normal approximation

btrecase.fdr <- lapply(btrecase, function(i) rbindlist(lapply(c(0.99, 0.95, .9,.85), function(j){
    dt <- rej.recode(0,i[rbias=="rbias",],c=j)
    dt[,PEP:=1-post.out]
    rej <- nrow(dt[null.rej=="no",])
    false.pos <- dt[null.rej=="no",sum(PEP)]
    return(data.table(post.level=j, total.rej=rej, total.false.pos=false.pos, FDR=false.pos/rej))
    })))

btrecase.fdr


## number of rejections using normal approx is higher than using posterior except for RNA 99CI

lapply(btrecase, function(i) rbindlist(lapply(c("null.99","null.95"), function(j){
    i[,.N, get(j)][,null:=j]
    })))

## calculate FDR using only associations shared with RASQUAL, Deseq2 and Gtex, very similar to FDR with all associations

fdr.shared <- rbindlist(lapply(unique(btrecase.gt.l$PIP), function(i) {
    dt <- btrecase.gt.l[PIP==i & null.rej=="no",]
    tmp <- data.table(post.level=i, total.rej=nrow(dt), total.fp=dt[,sum(PEP)], FDR=dt[,sum(PEP)]/nrow(dt))
    return(tmp)
}))

fdr.shared <- fdr.shared[order(FDR),]

###################################################
#### Make tables/plot  by number of associations
###################################################

all.ebv <- list(Rasqual=rasq.ebv, DEseq=dseq.ebv,Btrecase=btrecase.gt.ebv, lm=lm.ebv)
all.ebv2 <- list(Rasqual=rasq.ebv2, DEseq=dseq.ebv2,Btrecase=btrecase.gt.ebv2, lm=lm.ebv2)


tabs.dt <- function(l) { ## takes all.ebv or all.ebv2 and format data table for plotting
    
    tabs <- mapply(function(a,b,c) {
        dt <- tab2bplot(dt=a,var=c("Signif", b))
        tables <- lapply(sort(as.numeric(unique(dt[[b]]))), function(i) dt[get(b) ==i,])
        rbindlist(lapply(tables, function(i) format.tab(i, c)))

    },
    a=l,
    b=list("Fdr", "Fdr",  "PIP", "Fdr"),
    MoreArgs=list(c="Gtex-ebv"), SIMPLIFY=F)

    tabs.dt <- rbindlist(tabs, idcol="Method", fill=T)
    tabs.dt[, Significant.level:=as.numeric(Fdr)][!is.na(PIP),Significant.level:=as.numeric(PIP)]
    tabs.dt[, Discoveries:=Rasqual][!is.na(DEseq), Discoveries:=DEseq][!is.na(Btrecase), Discoveries:=Btrecase][!is.na(lm), Discoveries:=lm]
}

tabs1 <- tabs.dt(all.ebv)

tabs2 <- tabs.dt(all.ebv2)



##scaleFUN <- function(x) sprintf("%.2f", x)  ## for making consistent scales across plots

lab <- c("BaseQTL", "DESeq2","lm", "RASQUAL")


######## Plot by FDR ######

fdr.tabs <- lapply(list(tabs1, tabs2), function(i) {
    dt <- merge(i, btrecase.fdr$GT[,.(post.level, FDR)], by.x="Significant.level", by.y="post.level", all.x=T)
    dt[is.na(Fdr), Fdr:=FDR]
    return(approx.fdr(dt))
})


fdr.pa <- lapply(fdr.tabs, function(i) fdr.plot(i, lab, type="associations", col="Gtex-ebv"))

###############################
#### Make tables/plot by genes
###############################


sig.g <- function(l){ ## format for eGenes, input all.ebv or all.ebv2
    sig.genes <- mapply(function(a,b,c,d, f, j, k){   
        rbindlist(lapply(sort(as.numeric(unique(a[[b]]))), function(i) {
            r <- unique(a[get(b) == i & get(j) =="no", Gene_id])
            e <- unique(a[get(b) == i & get(k) =="no", Gene_id])
            dt <- data.table(Signif=c(d, f, c), sig.level=i, SNPs=c(length(setdiff(r,e)), length(intersect(r,e)), length(setdiff(e,r))))
            return(format.tab(dt,c))
        }))},
        a=l,
        b=list("Fdr", "Fdr",  "PIP", "Fdr"),
        j=list("null.Fdr", "null.fdr", "null.btrecaseGT", "null.Fdr"),
        k=list("null", "null", "null.ebv", "null"),
        d=names(all.ebv),
        MoreArgs=list(c="Gtex-ebv", f="Both"), SIMPLIFY=F)

    tabs.g <- rbindlist(sig.genes, fill=T, idcol="Method")

    tabs.g[, Discoveries:=Rasqual][!is.na(DEseq), Discoveries:=DEseq][!is.na(Btrecase), Discoveries:=Btrecase][!is.na(lm), Discoveries:=lm]

    setorder(tabs.g, Method,  sig.level)
}

tabs.g <- sig.g(all.ebv)
tabs.g2 <- sig.g(all.ebv2)


######## Plot by FDR ######

fdr.tabs.g <- lapply(list(tabs.g, tabs.g2), function(i) {
    dt <- merge(i, btrecase.fdr$GT[,.(post.level, FDR)], by.x="sig.level", by.y="post.level", all.x=T)
    dt[Method=="Btrecase", sig.level:=FDR][, Fdr:=sig.level]
    return(approx.fdr(dt))
})


fdr.pg <- lapply(fdr.tabs.g, function(i) fdr.plot(i,lab, "eGenes","Gtex-ebv" ))



## get legend to share between the 2 plots
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

myleg <- g_legend(fdr.pa[[1]])


## p.ebv.fdr <- mapply(function(a,b) {
##     plot_grid(a + theme(legend.position="none"),
##               NULL,
##               b + theme(legend.position="none"),
##               NULL,
##               myleg2,
##               ncol=5,
##               rel_widths=c(.8, 0.05, .8,0.05, 0.25),
##               labels=c('a', '', 'b', '', '')
##     ##label_x=c(0.1,0, 0.12)
##     )
## },
## a=fdr.pa,
## b=fdr.pg,
## SIMPLIFY=F
## )



##########################################################################################
## Comparing associations to GEUVADIS bigger sample size
########################################################################################


## GEUVADIS analysis from Chris
## select  chrom 22 and format compatible with gt
## sig.qtl <- fread('/home/ev250/rds/rds-cew54-wallace-share/Data/GEUVADIS/Geuvadis_sig_eqtl')
sig.qtl <- fread(snakemake@input[['geu_chris']])
sig.qtl <- sig.qtl[CHROM==22,][,Gene_id:=gsub("\\..*","", Gene)]
sig.qtl[, id := paste(POS,REF,ALT, sep=":")][, "null":="no"]
sig.qtl[, aFC:=2*Beta]


## Assess whether significant associations are in GEUVADIS bigger dataset with same direction of effects
geu.btrecase <- lapply(unique(gt.rna.long$source), function(i) {
    merge(gt.rna.long[source==i & type =="Signif" & rbias=="rbias" & abs(gene.dist) <= d,],
          sig.qtl,
          by.x=c("Gene_id", "tag"),
          by.y=c("Gene_id", "id"))
    })
names(geu.btrecase) <- unique(gt.rna.long$source)

common.assoc <- lapply(geu.btrecase, function(i) i[sign(log2_aFC_mean) == sign(Beta),.N,])
total.assoc <- lapply(unique(gt.rna.long$source), function(i) gt.rna.long[source==i & type =="Signif" & rbias=="rbias",.N,])

## 90% of assoc in GT are in GEU and 60% of no-GT are in GEU

mapply(function(a,b) round(a*100/b),
       a=common.assoc,
       b=total.assoc)

## look in no genotypes at different info cut-offs

common.sig.assoc.info <- sapply(seq(0.3,.9, 0.1), function(i) geu.btrecase[[2]][sign(log2_aFC_mean) == sign(Beta) & info>=i ,.N,])
total.sig.ngt <- sapply(seq(0.3,.9, 0.1),
                        function(i) gt.rna.long[source=='hidden-GT'  & type =="Signif" & rbias=="rbias" & abs(gene.dist) <= d &info>=i,.N])

info.geu <- ggplot(data=data.table(common=common.sig.assoc.info, total=total.sig.ngt, PPV=common.sig.assoc.info/total.sig.ngt,  info=seq(0.3,.9, 0.1)),
       aes(info, PPV)) + geom_point(shape=1) + geom_line(linetype="dashed", color="gray") + ylim(0,NA) +
    xlab("Imputation information threshold") +
    ylab("PPV") +
    ggtitle("Associations") +
    geom_label_repel(mapping= aes(label=total),show.legend=FALSE,size=2)

## look by eGenes


common.genes <- lapply(geu.btrecase, function(i) length(unique(i$Gene_id)))
total.genes <- lapply(unique(gt.rna.long$source), function(i) length(unique(gt.rna.long[source==i & type =="Signif" & rbias=="rbias",Gene_id])))

mapply(function(a,b) round(a*100/b),
       a=common.genes,
       b=total.genes)





total.sig.g.ngt <- lapply(seq(0.3,.9, 0.1),
                        function(i) gt.rna.long[source=='hidden-GT'  & type =="Signif" & rbias=="rbias" & abs(gene.dist) <= d &info>=i ,.N, Gene_id]$Gene_id)

##common.g.info <- sapply(total.sig.g.ngt, function(i) sum(i %in% unique(sig.qtl$Gene_id)))


common.g.info <- sapply(total.sig.g.ngt, function(i) sum(i %in% unique(geu.btrecase[[2]]$Gene_id)))

info.geu.g <- ggplot(data=data.table(common=common.g.info,
                                     total=sapply(total.sig.g.ngt, length),
                                     PPV=common.g.info/sapply(total.sig.g.ngt, length),
                                     info=seq(0.3,.9, 0.1)),
                     aes(info, PPV)) + geom_point(shape=1) + geom_line(linetype="dashed", color="gray") + ylim(0,NA) +
    xlab("Imputation information threshold") +
    ylab("PPV") + 
    ggtitle("eGenes")+
    geom_label_repel(mapping= aes(label=total),show.legend=FALSE,size=2)


geu.ngt.p <- plot_grid(info.geu, info.geu.g, ncol=1,labels="auto")

### SAve SuppFig6: hidden GT with GEUVADIS
ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/SuppFig6.png', geu.ngt.p, width=4.5, height=5)

##########################################################################
## Compare all methods to GEU bigger dataset, selecting same associations
##########################################################################


## start from rasq.l, dseq.l and btrecase.gt.l (same gene-SNP associations)

common.all.l <- list(Rasqual=rasq.l, DEseq=dseq.l, Btrecase=btrecase.gt.l, lm=lm.l)

common.all.l2 <- list(Rasqual=rasq.l2, DEseq=dseq.l2, Btrecase=btrecase.gt.l2, lm=lm.l2)


## get number of sig associations by FDR
tot.sig.a <- function(l) {##input common.all.l/l2
    rbindlist(mapply(function(a,b,c) {
        dt <- a[get(c) == "no",.N, get(b)]
        setnames(dt, "get", b)
    },
    a=l,
    b=c("Fdr", "Fdr", "PIP", "Fdr"),
    c=c("null.Fdr", "null.fdr", "null", "null.Fdr"),
    SIMPLIFY=F), idcol="Method", fill=T)
}


tot.sig.assoc <- tot.sig.a(common.all.l)
tot.sig.assoc2 <- tot.sig.a(common.all.l2)

    
                        
## merge with sig.qtl (bigger sample size)
all.geu <-  lapply(common.all.l , function(i) merge(i, sig.qtl, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "id")))

all.geu2 <-  lapply(common.all.l2 , function(i) merge(i, sig.qtl, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "id")))

## Count associations in same direction only, only interested in sig associations

sig.count <-function(l) {
    rbindlist(mapply(function(a, b, c, d) {
    dt <- a[sign(get(b)) == sign(get(c)), .N, d ]
    setnames(dt, names(dt)[2], "null")
    dt <- dt[null == "no",]
    return(dt)
},
a=l,
b=c("log2_aFC", "log2_aFC", "log2_aFC_mean", "log2.aFC_mean"),
d=list(c( "Fdr", "null.Fdr"),
       c( "Fdr", "null.fdr"),
       c( "PIP", "null.x"),
       c("Fdr", "null.Fdr")),
MoreArgs=list(c="Beta"),
SIMPLIFY=F), fill=T, idcol="Method")
}

sig.geu <- sig.count(all.geu)
sig.geu2 <- sig.count(all.geu2)


## Merge total with sig.geu and format, make plot

comp.geu.assoc <- help.geu(tot.sig.assoc, sig.geu,  btrecase.fdr$GT, N=all.geu$DEseq[,.N, Fdr][1,N] )

comp.geu.assoc2 <- help.geu(tot.sig.assoc2, sig.geu2,  btrecase.fdr$GT, N=all.geu2$DEseq[,.N, Fdr][1,N] )


## Try plotting similar FDR and joining FDR with lines

comp.fdr <- approx.fdr(comp.geu.assoc)
comp.fdr2 <- approx.fdr(comp.geu.assoc2)


fdr.a <- lapply(list(comp.fdr,comp.fdr2), function(i) fdr.plot(i, lab, "associations", "sig.geu"))


######################################################################################################################
## By genes
tot.sig.g <- function(l){
    rbindlist(mapply(function(a,b,c) {
    dt <- a[get(c) == "no", unique(Gene_id), get(b)][,.N, get]
    setnames(dt, "get", b)
},
a=l,
b=c("Fdr", "Fdr", "PIP", "Fdr"),
c=c("null.Fdr", "null.fdr", "null", "null.Fdr"), SIMPLIFY=F),
idcol="Method", fill=T)
}

tot.sig.genes <- tot.sig.g(common.all.l)

tot.sig.genes2 <- tot.sig.g(common.all.l2)


sig.geu.gene <- function(l){
    rbindlist(mapply(function(a, b, c, d) {
    dt <- a[sign(get(b)) == sign(get(c)), unique(Gene_id), d ][ ,.N,d]
    setnames(dt, names(dt)[2], "null")
    dt <- dt[null == "no",]
    return(dt)
},
a=l,
b=c("log2_aFC", "log2_aFC", "log2_aFC_mean", "log2.aFC_mean"),
d=list(c( "Fdr", "null.Fdr"),
       c( "Fdr", "null.fdr"),
       c( "PIP", "null.x"),
       c("Fdr", "null.Fdr") ),
MoreArgs=list(c="Beta"),
SIMPLIFY=F), fill=T, idcol="Method")
}

sig.geu.g <- sig.geu.gene(all.geu)
sig.geu.g2 <- sig.geu.gene(all.geu2)


comp.geu.g <- help.geu(tot.sig.genes, sig.geu.g, btrecase.fdr$GT, N=all.geu$DEseq[,unique(Gene_id), Fdr][,.N,Fdr][1,N] )

comp.geu.g2 <- help.geu(tot.sig.genes2, sig.geu.g2, btrecase.fdr$GT, N=all.geu2$DEseq[,unique(Gene_id), Fdr][,.N,Fdr][1,N] )


fdr.g <- lapply(list(comp.geu.g , comp.geu.g2), function(i) fdr.plot(approx.fdr(i), lab, "eGenes", "sig.geu")) 


##############################################
## Combine p.a, p.g, plots.a and plots.g

ex.all.fdr <- mapply(function(a,b,c,d){
    plot_grid(a + theme(legend.position="none"),
              NULL,
              b + theme(legend.position="non"),
              myleg,
              c + theme(legend.position="none"),
              NULL,
              d + theme(legend.position="non"),
              NULL, ncol=4,
              rel_widths=c(4.5,0.3, 4.5,1.5),
              labels=c('a', '', 'b', '', 'c', '', 'd', ''),
              label_x=rep(c(0.1,0, 0.12,0),2))
    
          },
    a= fdr.a,
    b=fdr.g,
    c=fdr.pa,
    d=fdr.pg,
    SIMPLIFY=F)
        


##ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/Fig3.png', ex.all.fdr, width = 6.4, height = 5.3)
## ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/Fig3alt.png', ex.all2, width = 12, height = 7)

mapply(function(i,j) ggsave(paste0(out.dir, "/",j, ".png"), i, width = 6.4, height = 5.3),
       i=ex.all.fdr,
       ##j=c("05", "01"),
       j=c("Fig3", "SuppFig2"),
       SIMPLIFY=F)


#######################
## eGENES as in Rasqual##
#######################

egenes <- fread(cmd=" zcat /mrc-bsu/scratch/ev250/EGEUV1/array_express_eqtl//EUR373.gene.cis.FDR5.best.rs137.txt.gz", header=F)

## from arrayX documentation I got header

names(egenes) <- c("SNP_ID", "ID", "Gene_id", "PROBE_ID", "chr_snp", "chr_gene", "SNPpos", "TSSpos", "Distance", "rvalue", "pvalue", "mlog10pvalue")

## select chrom 22

egenes <- egenes[chr_snp==22,]
egenes[, Gene_id:=gsub("\\.[0-9]*", "", Gene_id)]


## get genes by fdr and add col wether the gene is true positive or false positive

egenes.f <- function(l,assoc, egenes){
    ## l long list with methods with different FDR
    ## assoc data table with common associations across methods to define tested genes
    ## egenes data table with genes in gold standard

    ## get true positives: Genes in assoc and in egenes
    tp <- unique(assoc$Gene_id)[unique(assoc$Gene_id) %in% unique(egenes$Gene_id)]
    ## get true negative: Genes in assoc and not in egenes
    tn <- unique(assoc$Gene_id)[!unique(assoc$Gene_id) %in% unique(egenes$Gene_id)]

    rbindlist(mapply(function(a,b,c) {
        ## select method positives
        u <- unique(a[get(c)== "no" , c("Gene_id", b), with=F])
        ## get positive genes in gold standard and add aux col
        dt <- merge(u, data.table(Gene_id = unique(egenes$Gene_id), col=1:length(unique(egenes$Gene_id))), by="Gene_id", all.x=T)
        ## na values for col are false positives
        TP=dt[!is.na(col), .N, get(b)]
        names(TP) <- c(b, "TP")
        FP <- dt[is.na(col), .N, get(b)]
        names(FP) <- c(b, "FP")
        tmp <- merge(TP, FP, by=b, all=T)
        ## if tmp has NA convert to 0
        tmp[is.na(tmp)] <- 0
        ## add number of gold standard positives and negatives to calculate sensitivity and specificity
        tmp[,GEUpos:=length(tp)][,GEUneg:=length(tn)]
        tmp[, Sensitivity:=100*TP/GEUpos][, FPR:=100*FP/GEUneg]

        return(tmp)  
        
    },
    a=l,
    b=c("Fdr", "Fdr", "PIP", "Fdr"),
    c=c("null.Fdr", "null.fdr", "null", "null.Fdr"), SIMPLIFY=F),
    idcol="Method", fill=T)
}

tpfp <- egenes.f(common.all.l, assoc, egenes)
tpfp2 <- egenes.f(common.all.l2, assoc2, egenes)


ggplot(tpfp, aes(FPR, Sensitivity, color=Method))+ geom_line()
roc1 <- ggplot(tpfp2, aes(FPR, Sensitivity, color=Method))+ geom_line() + geom_point()

## eGenes with sig.qtl

tpfp.geu <- egenes.f(common.all.l, assoc, sig.qtl)
tpfp2.geu <- egenes.f(common.all.l2, assoc2, sig.qtl )

ggplot(tpfp.geu, aes(FPR, Sensitivity, color=Method))+ geom_line()
roc2 <- ggplot(tpfp2.geu, aes(FPR, Sensitivity, color=Method))+ geom_line() + geom_point()

plot_grid(roc1, roc2)

###########################
## QC on genotype errors ##
###########################

dna <- snakemake@input[['dna']]
rna <- snakemake@input[['rna']]

dna <- vcf_w("/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/DNA/chr22.ASE.allsamples.vcf.gz")#, f.arg='"%CHROM %POS %ID %REF %ALT[ %GT]\\n"')
rna <- vcf_w("/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/RNA/chr22.ASE.allsamples.vcf.gz")#, f.arg='"%CHROM %POS %ID %REF %ALT[ %GT]\\n"')

## select fsnps for genes run in RNA, either with ASE counts or used for imputation:

fsnp <- fread(snakemake@input[['eSNPs']])
## fsnp <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt")
fsnp[,id:=paste(POS, REF,ALT, sep=":")]

fsnp <- fsnp[gene_id %in% unique(btrecase$RNA[['Gene_id']]),]

## get fsnps used in inference in GT

fsnp.i <- lapply(list.files(btrecase_dir[1], pattern="rbias.ENSG[0-9]+.*.\\.GT.fsnps.with.counts.rds", full.names=T), readRDS)

fsnp.i <- unique(unlist(fsnp.i))

##remove .n .m ending

fsnp.i <- unique(unlist(lapply(fsnp.i, function(i) gsub(".[n:m]$", "", i))))

## select fsnp used in inference for DNA and in genes run with RNA, both in DNA and RNA

fsnp.i <- fsnp[id %in% fsnp.i,]

dna <- merge(dna, fsnp.i[,.(gene_id, id)], by="id")
rna <- merge(rna, fsnp.i[,.(gene_id, id)], by="id")

## For the snps called in dna and rna, see concordance in GT

rna.dna.c2 <- convert2(dna,rna,chr=22)

rna_errors <- new_err(rna.dna.c2, x="median")

## Look at the genotype for fsnps according to fisher test
  
## fisher <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/QCfSNP/fisher.fSNPs.txt")
fisher <- fread(snakemake@input[['het_prop']])


## add Fisher p value to rna
rna.f <- merge(rna, fisher[pvalue>0.01,] , by.x=c("id", "gene_id"), by.y=c("fsnp", "gene_id"))

## Compare rna_errors before Fisher QC with rna_errors.f after Fisher QC in a xtable

######################
## SuppTable 2 ##
######################


err.sum <- tot_err(DT=rna.dna.c2)

st2 <- xtable(err.sum, caption= "Error rates and missing genotypes for genotyping by RNA-seq. For the fSNPs used in inference with observed genotypes, the number of erroneous calls and missing genotypes across all samples is shown.")

print(st2 ,include.rownames = FALSE, booktabs = TRUE)
              

##comb <- comp.x(dt1=rna_errors, dt2=rna_errors.f, cap="Error rates and missing genotypes for genotyping by RNA-seq. For each sample it was calculated the number of correct genotype calls, number of errors in one chromosome or both and number of mssing genotypes. The table shows the median, interquartile range and percentage across samples before and after we controlled for sample frequency of heterozygosity to the reference panel (QC, Online Methods)")


##################
## Genotyping QC
##################

##  For each snp add the number of samples with hom->het error with ASE >=5 (cut-off) ####

### get ASE counts per individual (n, m cols)
c.ase <- tot.ase_counts(x=rna.f)
dna.ase <- tot.ase_counts(x=dna)

## remove snps  according to cut=offs
c.ase <- filt.fsnp(c.ase,ase=5,min.ase.snp=5,n=5, rem=NULL)

### get m, total ase counts

## rna
m.ase <- t(c.ase[, grep("m", colnames(c.ase), value=T)])
rownames(m.ase) <- sub("\\.m", "", rownames(m.ase))
colnames(m.ase) <- paste0(colnames(m.ase), "_ASE")

## dna

m.dna <- t(dna.ase[, grep("m", colnames(dna.ase), value=T)])
rownames(m.dna) <- sub("\\.m", "", rownames(m.dna))
colnames(m.dna) <- paste0(colnames(m.dna), "_ASE")

### get p = n/m
n.ase <- t(c.ase[, grep("n", colnames(c.ase), value=T)])
rownames(n.ase) <- sub("\\.n", "", rownames(n.ase))
colnames(n.ase) <- paste0(colnames(n.ase), "_p")

p.ase <- n.ase/m.ase


## dna
n.dna <- t(dna.ase[, grep("n", colnames(dna.ase), value=T)])
rownames(n.dna) <- sub("\\.n", "", rownames(n.dna))
colnames(n.dna) <- paste0(colnames(n.dna), "_p")

p.dna <- n.dna/m.dna


### replace AS columns by ASE and p colums

## rna
rna.f[ , grep("AS", names(rna.f), value=T):=NULL]

rna.f <- merge(rna.f, data.table(cbind(m.ase, p.ase), keep.rownames=T), by.x="id", by.y="rn")

## ASE, p from RNA
rna.dna.f <- merge(rna.f, dna, by=c("CHROM","POS","REF","ALT", "id", "gene_id", "ID"), suffixes=c("_RNA","_DNA"))

## dna
## select fsnps from p.dna

p.dna <- p.dna[rownames(p.dna) %in% fsnp$id,]
m.dna <- m.dna[rownames(m.dna) %in% fsnp$id,]


## Add m.dna and p.dna to dna and merge dna with rna

## replace AS columns by ASE and p colums
dna[ , grep("AS", names(dna), value=T):=NULL]

## only keep fSNPs 
dna <- merge(dna, data.table(cbind(m.dna, p.dna), keep.rownames=T), by.x="id", by.y="rn")

## get number of samples with errors per error type ("hom2het", "het2hom", "hom2hom"), missing Gt and percentage of correct GT (correc.per) excluding missing.
rna.dna.err  <- convert3(rna.dna.f)

error <- c("hom2het", "het2hom", "hom2hom", "missing")

################################################################
## look at overall errors for all fsnps including missing data

miss2correct <- ggplot(rna.dna.err[!is.na(gene_id),], aes(correct.per, missing*100/length(grep("AS$", names(rna.dna.err), value=T)))) +
    ylab("Missing GT (%)") +
    xlab("Correct GT (%)" ) +
    geom_rug(col=rgb(.7,0,0,alpha=.2))+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=12)) +
    geom_point(shape=1)


fsnps.p <- ggplot(gt.rna.100[,.SD[1],Gene_id] , aes(x=as.numeric(n.fsnps.gt),y=n.fsnps.rna)) + geom_count(aes(color = ..n.., size = ..n..)) +
    guides(color = 'legend')+
    geom_abline(slope=1, intercept=0, linetype="dashed", color = "gray")+
    xlab("fSNPs (DNA-sequencing)") +
    ylab("fSNPs (RNA-seq)" ) + theme(legend.title = element_blank())+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=12)) +
    scale_colour_gradientn(colours = heat.colors(8))

## alternative plot, select fsnps for gt and rna for each gene (83 in total)

dt.f.100 <- gt.rna.100[,.SD[1],Gene_id][,.( gt=as.numeric(n.fsnps.gt),rna=n.fsnps.rna)]

## make categories 


dt.f.100[, gt.cat:=cut(gt, breaks=sort(c(2, seq(1,max(gt)+5, 5))), right=F, labels=c("1", "2-5", "6-10", "11-15", "16-20", ">20"))]
dt.f.100[, rna.cat:=cut(rna, breaks=sort(c(2, seq(1,max(rna)+5, 5))), right=F, labels=c("1", "2-5", "6-10", "11-15", "16-20", ">20") )]
dt.f.100[, rna.cat.inv:=factor(rna.cat, levels=rev(levels(dt.f.100$rna.cat)))]
dt.f <- dt.f.100[,.N, .(gt.cat, rna.cat, rna.cat.inv)]

f.plot <- ggplot(dt.f, aes(x=gt.cat, y=N, alpha=rna.cat)) + geom_bar(stat='identity',  color = "#209f1b", fill = "#35665C") +
    geom_text(aes(y =8, label = N), color = "red", alpha=1, size =4) +
    facet_wrap(~rna.cat.inv, ncol=1, strip.position = "left") +
    labs(title = "Number of genes",
         y = "fSNPs RNA-seq",
         x = "fSNPs DNA-seq") +
    theme_minimal() +
    theme(strip.background = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 12),
          strip.text.y = element_text(angle = 180),
          strip.text = element_text(size = 12),
          legend.position = "none")



## percentage of excluded fSNPs in RNA due to Fisher correction

## look at dt.f.100
## nf <- colSums(dt.f.100[,.(gt,rna)])

## 100 - nf['rna']*100/nf['gt']

###############################################################################
## look at GT errors by depth (from RNA.genotyping.QC)
## justify why DP=10
##############################################################################

## rna_dp <- name('/mrc-bsu/scratch/ev250/EGEUV1/quant/call_vars/b37/HG96_2215_chr22_q20.txt')
rna_dp <- name(snakemake@input[['rna_dp']])


## Add fisher p-value to fsnps:

rna_dp[, id:=paste(POS, REF, ALT, sep=":")]
rna_dp <- merge(rna_dp,  fisher, by.x="id", by.y="fsnp")

rna.dna.dp <- convert2(dna, rna_dp, chr=22)

##rna.dna.dp2 <-  convert3(rna.dna.dp)

hets.p <- err2.dp(DT=rna.dna.dp,GT="het", x=50)
##setnames(hets.p, "p", "Het")
hets.p[,GT:="Heterozygous"]
hom.p <- err2.dp(rna.dna.dp,GT="homo", x=50)
##setnames(hom.p, "p", "Hom")
hom.p[, GT:="Homozygous"]

## get cumulative proportion of fSNPs by depth
cum.fsnps <- fsnps.depth(DT=rna.dna.dp)

## format to make it compatible with hets.p and hom.p

cum.fsnps <- cum.fsnps[,.(fsnps.cum.prop, Depth)]
names(cum.fsnps) <- c("Calls", "DP")
cum.fsnps[, Calls:= 1- Calls]
cum.fsnps[,lab:="Called fSNPs"]
cum.fsnps <- cum.fsnps[DP <=51,]

comb <- rbind(hets.p, hom.p)

## comb <- Reduce(function(x,y) merge(x, y, by="DP"), list(hets.p, hom.p,cum.fsnps ))

DP.p <- ggplot() +
    geom_line(data=comb, mapping=aes(x=DP, y=p, color=GT)) +
    geom_line(data=cum.fsnps, mapping=aes(x=DP, y=Calls, color=lab)) +
    scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Calls (proportion)"))+
    scale_colour_manual(values = c("black", "blue", "red" )) +
    ##ggtitle("Genotype concordance") +
    labs(x="Depth",
         y="Correct calls (proportion)",
         colour="")  +
    geom_vline(xintercept=10, linetype="dashed", color="gray") +
    guides(colour = guide_legend(reverse=TRUE)) +
    theme(legend.position = c(0.4, 0.8))

### genotyping errors by fisher p-value for het proportions selecting fSNPs with DP>=10

err.fisher <- err2.var(col="pvalue", x=10, DT=rna.dna.dp)

error.fisher.prop <- ggplot(data=err.fisher) +
    geom_point(mapping=aes(y=1-prop.conc, x=-log(pval, 10))) +
    labs(y="Proportion of errors",
         x=expression(paste(-log[10], "(p-value)"))) + geom_vline(xintercept=-log(0.012, 10), linetype="dashed", color="gray")
   

    
qc.g <- plot_grid(DP.p, NULL, error.fisher.prop, nrow=1, labels=c('a','', 'b'), rel_widths=c(1,.2,1))
###################################################################################
##  suppFig3
##################################################################################
ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/SuppFig3.png', qc.g, width = 8.2, height = 3.5)


## add  c("hom2het", "het2hom", "hom2hom", "missing", "id", "correct.per") from rna.dna to rna.dna.p, so I can have missing data, etc with DP

## rna.dna.dp2 <- merge(rna.dna.dp, rna.dna[, c(error, "correct.per", "id"), with=F], by.x="id_RNA", by.y="id")


######################################################################################################################################
## Look at false positives and false negatives comparing rna eQTL effects to equivalent run but with correct genotypes of fSNPs (QC2)

qc1_dir <- snakemake@params[['QC1']]
qc2_dir <- snakemake@params[['QC2']]
## qc1_dir <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/QC1_F/RNA'
## qc2_dir <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/QC2_F/RNA'

qc1 <- comb.files(path=qc1_dir, pattern="ENSG[0-9]+.*.noGT.stan.summary.txt")
qc2 <- comb.files(path=qc2_dir, pattern="ENSG[0-9]+.*.noGT.stan.summary.txt")

## remove bad reads, bad info calls
qc1 <- qc1[Rhat<maxRhat,]
qc1 <- qc1[info<=max.info,]
qc1 <- add.null(qc1)

qc2 <- qc2[Rhat<maxRhat,]
qc2 <- qc2[info<=max.info,]
qc2 <- add.null(qc2)

##
saveRDS(qc1, "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/qc1.rds")
saveRDS(qc2, "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/qc2.rds")


s <- c("fSNP.DNA","fSNP.RNA")
qc2.rna <- merge(btrecase$RNA[rbias=="rbias",], qc2, by=c("Gene_id", "tag"), suffixes=rev(s))
qc2.rna <- add.signif(qc2.rna, "null.99fSNP.DNA" , "null.99fSNP.RNA" , s)


cols <- c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', 'null.99')
xlab <- "eQTL-effect (fSNP DNA)"
ylab <- "eQTL-effect (fSNP RNA)"

qc2.rna.tab <- tab2bplot(dt=qc2.rna, colors=setNames(c("#999999", "yellow3", "#0072B2", "#D55E00"),
                                                     c("None",s, "Both")))


## Look at genes with false positive associations and compare with genes that were correctly identified

fpg <- unique(qc2.rna[Signif=="fSNP.RNA", Gene_id])

correct <- unique(qc2.rna[Signif=="Both", Gene_id])[!unique(qc2.rna[Signif=="Both", Gene_id]) %in% fpg]

fn <- unique(qc2.rna[Signif=="fSNP.DNA", Gene_id])[!unique(qc2.rna[Signif=="fSNP.DNA", Gene_id]) %in% c(correct, fpg)]

## Get missing fSNPs due to low depth in all samples:
miss.f <- lapply(list(fpg, correct, fn), function(j){
    tmp <- unique(qc2.rna[Gene_id %in% j, .(Gene_id, n.fsnpsfSNP.DNA, n.fsnpsfSNP.RNA)])
    tmp[, miss.fsnp:=n.fsnpsfSNP.DNA-n.fsnpsfSNP.RNA]
    })


## for each potential fgp/fn look at errors in fsnps

err.fp <-rbindlist(mapply(function(i,j,l) {
    ## for each list (i) take each gene (k)
    rbindlist(lapply(i, function(k){
        fpid <- fsnp.i[gene_id ==k, id]
        dt <- tot_err(DT=rna.dna.c2[id_RNA %in% fpid,])
        dt[, Gene_id:=k]
        dt[, type:=l]
        ## add missing fSNP (x86 samples)
        dt[Errors=="Missing",N:=N + j[Gene_id==k, miss.fsnp*86]]
        dt[, `%`:=round(N*100/sum(N),2)]
    })
    )
},
i=list(fpg, correct, fn),
j=miss.f,
l=c("Potential false positive", "Correct", "Potential false negative"),
SIMPLIFY=F))



    
## fpg[1] shows somehow biassed estimates for some associations, indicate this gene in col fpb

## qc2.rna[,fpb:=ifelse(Gene_id==fpg[1], "yes", "no")]

e.fsnp <- btrecase.plot(dt=qc2.rna[ Signif !="None",],
                      x1=paste0(cols, s[1]),
                      x2= paste0(cols, s[2]),
                      xl=xlab, yl=ylab, col=s, 
                      title="Effect of fSNP genotype error\n on eQTL effects"
                       ) +
    annotation_custom(tableGrob(qc2.rna.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(qc2.rna.tab$color, rep("black", 4)))))), xmin=-1.2, xmax=-.6, ymin=0.3, ymax=0.5) +
    ## add arrow to indicate fp and fn
    annotate("segment", x=c(0.09, .12, -.765), xend=c(0.09, -.04, -.765), y=c(.3, -.18, -.15), yend=c(.20, -.15, -.065),
             colour=c(rep("#0072B2",2),"orange2" ), arrow=arrow(length = unit(0.1, "cm"))) 
    
print(cor.test(qc2.rna[ Signif != "None" ,log2_aFC_meanfSNP.RNA], qc2.rna[ Signif != "None" ,log2_aFC_meanfSNP.DNA]))

## eGenes, no false positives as fpg are in correct

df.venn <- data.table(x=c(-1,-1), y=rep(0,2), labels=rev(s))
common <- unique(qc2.rna[Signif=="Both", Gene_id])
pos.dna <- unique(qc2.rna[null.99fSNP.DNA=="no", Gene_id])

vqc2 <- ggplot(df.venn) +
    geom_circle(aes(x0=x,y0=y,r=c(2, 1), fill=labels, color=labels),  alpha=0.25, size=.5) +
    geom_rect(mapping=aes(xmin=-3.5,xmax=2, ymin=-2.7, ymax=2.7, fill="Total genes"), alpha=0.01, color="black") +
    coord_fixed() +
    theme_void() +
    scale_fill_manual(values = c(cols.d, "white") ) +
    scale_colour_manual(values = c(cols.d, "black"), guide = FALSE) +         
    labs(fill = NULL) +
    annotate("text", x= c(-1, -1, -1), y=c(1.3,0, 2.4), label=c(unlist(lapply(list(pos.dna, common),length)),length(unique(qc2.rna$Gene_id))), size=5, colour=c("orange2", "#0072B2", "black")) +
    
    ggtitle("Overlap of genes with significant eQTLs") +
    theme(plot.title = element_text(hjust = 0.5))



######################################################################################################################################
## Look at false positives and false negatives comparing imputed vs GT for cis-SNP (same genotypes of fSNPs, QC1)


## merge gt.rna with qc1, select gt cols and tag.rna from gt.rna

qc1.gt <- merge(gt.rna[ ,c("Gene_id", "tag.rna", "op.dir", grep("gt$", names(gt.rna), value=T)), with=F], qc1, by.x=c("Gene_id", "tag.rna"), by.y=c("Gene_id", "tag"))


sig.lab <- c("Observed-GT", "Hidden-GT")
qc1.gt <- add.signif(qc1.gt, "null.99.gt" , "null.99" , sig.lab)


## correct sign of qc1 when op.dir == "yes"

qc1.gt[op.dir=="yes", log2_aFC_mean:= -log2_aFC_mean]



## Look for inconsistencies due to changes in eaf being too close to 0.5

change <- qc1.gt[(Signif!="None"  & sign(log2_aFC_mean.gt) != sign(log2_aFC_mean) ), .(tag.gt, tag.EAF.gt, tag.rna, tag.EAF, Gene_id)]

## qc1.gt[Gene_id == "ENSG00000198951" & tag.rna == "42405922:T:C" & tag.gt== "42370991:T:C", log2_aFC_mean:= -log2_aFC_mean]
## qc1.gt[Gene_id == "ENSG00000198951" & tag.rna == "42405922:T:C" & tag.gt== "42370991:T:C",
##        grep("[0-9]%$", names(qc1.gt), value=T):= lapply(grep("[0-9]%$", names(qc1.gt), value=T), function(i) -get(i))]

for (r in 1:nrow(change)){
    qc1.gt[Gene_id == change[r,Gene_id] & tag.rna == change[r,tag.rna] & tag.gt == change[r,tag.gt], log2_aFC_mean:= -log2_aFC_mean]
    qc1.gt[Gene_id == change[r,Gene_id] & tag.rna == change[r,tag.rna] & tag.gt == change[r,tag.gt],
           grep("[0-9]%$", names(qc1.gt), value=T):= lapply(grep("[0-9]%$", names(qc1.gt), value=T), function(i) -get(i))]
}


cols <- c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', 'null.99')
xlab <- "eQTL-effect Obs-GT"
ylab <- "eQTL-effect Imp-GT"

qc1.gt.tab <- tab2bplot(dt=qc1.gt, colors=setNames(c("#999999", "yellow3", "#0072B2", "#D55E00"),
                                                    c("None", sig.lab, "Both")))

e.cis <- btrecase.plot(dt=qc1.gt[ Signif !="None" ,],
                      x1=paste0(cols, ".gt"),
                      x2= paste0(cols, ""),
                      xl=xlab, yl=ylab, col=sig.lab, 
                      title="Effect of cis-SNP imputation\n on eQTL effects") +
    annotation_custom(tableGrob(qc1.gt.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(qc1.gt.tab$color, rep("black", 4)))))), xmin=-1, xmax=-.2, ymin=0.1, ymax=0.5)


## plot e.cis.info for different info cut-offs
inf <- seq(0.3, .9, 0.2)

f.info <- function(inf){
    qc1.gt.info.tab <- tab2bplot(dt=qc1.gt[info >=inf ,], colors=setNames(c("#999999", "yellow3", "#0072B2", "#D55E00"),
                                                                         c("None", sig.lab, "Both")))
    ##print(qc1.gt[info >=inf & Signif !="None",.N, .(Gene_id,Signif)])
    
    
    e.cis.info <- btrecase.plot(dt=qc1.gt[ Signif !="None" & info >=inf ,],
                                x1=paste0(cols, ".gt"),
                                x2= paste0(cols, ""),
                                xl=xlab, yl=ylab, col=sig.lab, 
                                title=paste("Effect of imputation on eQTL effects\n info >=", inf)) +
        ylim(NA,0.6) +
        
        annotation_custom(tableGrob(qc1.gt.info.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(qc1.gt.info.tab$color, rep("black", 4)))))), xmin=-1, xmax=-.2, ymin=0.1, ymax=0.6)
}

info.p <- lapply(inf, function(i) f.info(i))


info.all <- plot_grid(plotlist=info.p, ncol=2)
ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/SuppFig5info.png', info.all, width = 11, height = 9.8)

## e.f.c <- plot_grid(e.fsnp, e.cis, e.cis.info, ncol=3, labels='auto')
e.f.c <- plot_grid(e.fsnp, e.cis, ncol=2, labels='auto')

ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/SuppFig5.png', e.f.c, width = 15, height = 5)


##################################################################################################
############  Look at observed AI in hets from GT versus RNA, need ASE from DNA as well  ###################

## Select fSNPs with at least m = 20, binom.test(0, 20)$p.value = 1.907349e-06
#### rna

m <- 20

rna_GT <- grep("GT_RNA", names(rna.dna.c2), value=TRUE)
prefix <- gsub("GT_RNA", "",rna_GT)   ## sample

## ## allelic imbalance for all hets in rna, correct p by GT

hets.rna <-  unlist(lapply(prefix, function(i) c(unlist(rna.dna.f[( get(paste0(i,"GT_RNA")) == "0|1") & get(paste0(i,"ASE")) >= m  ,paste0(i,"p"), with=F] ),
                                                 1- unlist(rna.dna.f[get(paste0(i,"GT_RNA")) == "1|0" & get(paste0(i,"ASE")) >= m  ,paste0(i,"p"), with=F] ) )))


## RNA: look at allelic imbalance in  hets: total and by error type

## hom2het.ai <- unlist(lapply(prefix, function(i) c(unlist(rna.dna[( get(paste0(i,"GT_RNA")) == "0|1") & get(paste0(i,"ASE")) >= m  & (get(paste0(i,"GT_DNA")) == "0|0" |  get(paste0(i,"GT_DNA")) == "1|1" )  ,paste0(i,"p"), with=F] ),
##                                                  1- unlist(rna.dna[get(paste0(i,"GT_RNA")) == "1|0" & get(paste0(i,"ASE")) >= m  & (get(paste0(i,"GT_DNA")) == "0|0" |  get(paste0(i,"GT_DNA")) == "1|1" )  ,paste0(i,"p"), with=F] ) )))


## het.ok.ai <- unlist(lapply(prefix, function(i)  c(unlist(rna.dna[( get(paste0(i,"GT_RNA")) == "0|1") & get(paste0(i,"ASE")) >= m  & (get(paste0(i,"GT_DNA")) == "1|0" | get(paste0(i,"GT_DNA")) == "0|1")    ,paste0(i,"p"), with=F] ),
##                                                   1-unlist(rna.dna[ get(paste0(i,"GT_RNA")) == "1|0"  & get(paste0(i,"ASE")) >= m  & (get(paste0(i,"GT_DNA")) == "1|0" | get(paste0(i,"GT_DNA")) == "0|1")    ,paste0(i,"p"), with=F] ))))



## rna.ai <- data.table(type=rep(c("hom2het", "het"), time= c(length(hom2het.ai), length(het.ok.ai))),
##                         AI=c(hom2het.ai, het.ok.ai))



hets.dna <- unlist(lapply(prefix, function(i) c(unlist(dna[( get(paste0(i,"GT")) == "0|1") & get(paste0(i,"ASE")) >= m  ,paste0(i,"p"), with=F] ),
                                                1- unlist(dna[get(paste0(i,"GT")) == "1|0" & get(paste0(i,"ASE")) >= m  ,paste0(i,"p"), with=F] ) )))



rna.dna2 <- merge(rna, dna, by=c("CHROM","POS","REF","ALT", "id"), suffixes=c("_RNA","_DNA"))

## rna, also look at AI for hets called hom (need to look at AI from dna, as coded 0 in rna (hom))

het2hom.ai <- unlist(lapply(prefix, function(i) c(unlist(rna.dna2[( get(paste0(i,"GT_RNA")) == "0|0" | get(paste0(i,"GT_RNA")) == "1|1") & get(paste0(i,"ASE")) >= m  & ( get(paste0(i,"GT_DNA")) == "0|1" )  ,paste0(i,"p"), with=F] ),
                                                  1-unlist(rna.dna2[( get(paste0(i,"GT_RNA")) == "0|0" | get(paste0(i,"GT_RNA")) == "1|1") & get(paste0(i,"ASE")) >= m  & (get(paste0(i,"GT_DNA")) == "1|0" )  ,paste0(i,"p"), with=F] )
                                                  )))

rna.ai <- rbind(rna.ai, data.table(type="het2hom", AI=het2hom.ai))

## plots, AI in rna, AI rna vs dna

s <- 0.8 ## size line in plot
## AI.rna <- ggplot(rna.ai, aes(qlogis(AI), colour=type)) +
##     theme_bw() +
##     xlab("Allelic Imbalance (logit(Alt. reads/ Total reads))") +
##     ylab("Density") + theme(legend.title = element_blank()) +     
##     guides(colour=guide_legend(override.aes = list(size = 1))) +
##     scale_colour_manual(values=c("red", "#0072B2", "black"))+
##     geom_vline(xintercept=qlogis(0.5), linetype="dashed", color="gray") +
##     geom_density(size=s) 

p.dna.rna <- data.table(type=rep(c("DNA", "RNA"), time= c(length(hets.dna), length(hets.rna))),
                        AI=c(hets.dna, hets.rna))

AI.dna.rna <- ggplot(p.dna.rna, aes(qlogis(AI), colour=type, linetype=type)) +
    theme_bw() +
    xlab("AI (logit(Alt. reads/ Total reads))") +
    ylab("Density") + theme(legend.title = element_blank()) +
    xlim(-2,2) +
    guides(colour=guide_legend(override.aes = list(size = 1))) +
    geom_density(size=s) +
    geom_vline(xintercept=qlogis(0.5), linetype="dashed", color="gray") +
    scale_colour_manual(values=c("#FFCC33", "#0072B2"))

## AI.dr <- plot_grid(miss2correct,fsnps.p,
##                    NULL, NULL,                           
##                    AI.dna.rna, AI.rna,
##                    nrow=3, labels=c("a", "b", "", "", "c", "d"),  rel_heights=c(1,0.2,1))

## alternative

gterr.missing <- plot_grid(miss2correct,f.plot, AI.dna.rna, ncol=1, labels='auto')
                   ## NULL, NULL,                           
                   ## AI.dna.rna, AI.rna,
                   ## nrow=3, labels=c("a", "b", "", "", "c", "d"),  rel_heights=c(1,0.2,1))


ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/SuppFig4.png', gterr.missing, width = 5.6, height = 6.7)


################################################################################################################################################
#######################################  Ref bias #################################################################################

## Distribution of allelic imbalance estimates

AI <- fread(snakemake@input[['AI']])
## AI <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/post_remap/pre_post_AI.txt")

## Calculate 99% confidence interval testing whether AI_post is different from 0.5

AI[Total_post >0 , c("lowCI", "hCI") := lapply(1:2, function (i) binom.test(x=NALT_post, n=Total_post, p = 0.5, alternative = "two.sided", conf.level = 0.99)$conf.int[i]),by=.(NALT_post,Total_post) ]

## Add 99% CI relative to 0.5

AI[ , c("null_lCI", "null_hCI") := lapply(c("lowCI", "hCI"), function(i)  0.5 + get(i) -AI_post) ]

AI[, Discard:= ifelse(AI_post > null_hCI, "yes", "no")][Total_pre<100, Discard := "yes"]

AI[,.N,Discard]

## Indicate fSNPs in DNA and RNA that can be used in estimation (hets with ASE)

dna.f <- colSums(dna.ase)[colSums(dna.ase) >0]
dna.f <- dna.f[grep(".m$", names(dna.f))]
dna.f <- setNames(dna.f, gsub("\\.m", "", names(dna.f)))

AI[, id:=paste(POS,REF,ALT, sep=":")]

AI <- merge(AI, data.table(id=names(dna.f), DNA.fsnp=names(dna.f)), by="id", all.x=T)
AI[!is.na(DNA.fsnp) , DNA.fsnp:="yes",][is.na(DNA.fsnp), DNA.fsnp:="no",]

## Add column indicating whether fSNP was used in inference for DNA (Discard =no and DNA.fsnp=yes)

AI[ , Inference.DNA:="No"][Discard=="no" & DNA.fsnp == "yes", Inference.DNA:="Yes"]

## same for RNA, use rna.f

AI <- merge(AI, data.table(id=rna.f$id, RNA.fsnp=rna.f$id), by="id", all.x=T)
AI[!is.na(RNA.fsnp) , RNA.fsnp:="yes",][is.na(RNA.fsnp), RNA.fsnp:="no",]

## Add column indicating whether fSNP was used in inference for DNA (Discard =no and DNA.fsnp=yes)

AI[ , Inference.RNA:="No"][Discard=="no" & RNA.fsnp == "yes", Inference.RNA:="Yes"]


labs <- setNames(c("Excluded", "Inference"), c("No", "Yes"))

## Draw lines for 99% CIs from the null value
## to plot AI_post in logit scale I convert 1 to 0.99 and 0 to 0.01
AI[, AI_post.plot:=AI_post][AI_post>0.99, AI_post.plot:=0.99][AI_post<0.01, AI_post.plot:=0.01]

EAI.dna <- ggplot(AI[Total_post>0 & Inference.DNA == "Yes" ,] ) +
    geom_line(aes(x=log(Total_pre), y=qlogis(null_hCI)),  colour="gray", linetype="dashed", na.rm=TRUE) +
    geom_point(aes(x=log(Total_pre), y=qlogis(AI_post.plot)), shape=1, color="blue") +
    #geom_point(aes(x=log(Total_pre), y=AI_post, colour=Inference.DNA, shape=Inference.DNA)) +
    ##scale_shape_manual(values=c(1,16)) +
    geom_vline(xintercept=log(100), linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    ylab("Reference bias estimate") +
    xlab("Depth (log)") +
    scale_y_continuous(breaks=c(-3,0,3,6), limits=c(-4,7)) +
    scale_x_continuous(breaks=seq(6,12,2), limits=c(1, 14)) +
    geom_rug(aes(x=log(Total_pre), y=qlogis(AI_post.plot)), color="blue", alpha=0.1) + #rgb(.2,.8,.3,alpha=.2)) +
    #xlab("Reads before re-alignment (log)") +
    ## scale_colour_manual(values=c("grey72", "#CC79A7")) +
    ## labs(colour = "Inference", shape ="Inference") +
    ## facet_grid(~ Inference.DNA, labeller=labeller(Inference.DNA=labs))  +
    theme(strip.background =element_rect(fill="white", color="black", linetype="solid"))


EAI.rna <- ggplot(AI[Total_post>0 & Inference.RNA == "Yes",] ) +
    #geom_line(aes(x=log(Total_pre), y=null_lCI), colour="gray", linetype="dashed") +
    geom_line(aes(x=log(Total_pre), y=qlogis(null_hCI)),  colour="gray", linetype="dashed") +
    geom_point(aes(x=log(Total_pre), y=qlogis(AI_post.plot)), shape=1, color="blue") +
    ##geom_point(aes(x=log(Total_pre), y=AI_post, colour=Inference.RNA, shape=Inference.RNA)) +
    scale_shape_manual(values=c(1,16)) +
    geom_vline(xintercept=log(100), linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    ylab("Reference bias\nestimate") +
    xlab("Depth (log)") +
    scale_y_continuous(breaks=c(-3,0,3,6), limits=c(-4,7)) +
    scale_x_continuous(breaks=seq(6,12,2), limits=c(1, 14)) +
    geom_rug(aes(x=log(Total_pre), y=qlogis(AI_post.plot)), color="blue", alpha=0.1) +
    ## facet_grid(~ Inference.RNA, labeller=labeller(Inference.RNA=labs))  +
    theme(strip.background =element_rect(fill="white")) 
    ## scale_colour_manual(values=c("grey72", "#CC79A7")) +
    ## labs(colour = "Inference", shape ="Inference") #+
    #ggtitle("Allelic imbalance estimates")



## ## plot histograms for hets.dna, hets.rna separately to add to ref panel figure
## hets.p <- mapply(function(i, j) {
   
##     ggplot(p.dna.rna[type==i,], aes(AI)) +
##         theme_bw() +
##         xlab(j) +
##         ylab("Density") + #theme(legend.title = element_blank()) +     
##                                         #guides(colour=guide_legend(override.aes = list(size = 1))) +
##         geom_density(color="red") +
##         geom_vline(xintercept=0.5, linetype="dashed", color="gray") +
##         coord_cartesian(ylim=c(0, 6))
##     },
##     i= unique(p.dna.rna$type),
##     j=c("", "Observed AI"),
##     SIMPLIFY=F)



########################################################################################
## Compare eQTL estimates with and without reference bias correction, GT only ASE assoc

r.bias <- lapply(btrecase, function(i) {
    u <- unique(i$rbias)
    print(u)
    names(u) <- c("Without", "With")
    print(u)
    x <- paste("null.99", u, sep=".")
    i <- i[model != "trec",]
    dt <- merge(i[rbias == u[1],], i[rbias ==u[2], ], by=c("Gene_id", "tag", "tag.EAF", "gene.dist"), suffixes= paste0(".",u))
    dt <- add.signif(dt, x1=x[1], x2=x[2], col=names(u))
    return(dt)
    })


cols <- c("log2_aFC_mean", "log2_aFC_0.5%", "log2_aFC_99.5%", "null.99")
colors.p <- c("#999999", "green4", "orchid4", "salmon")

## Plots 

rbias.p <- mapply(function(i,j,k) {
    dt <- i[abs(gene.dist) <=d,]
    tab <- tab2bplot(dt, colors=setNames(colors.p, levels(dt$Signif)))
    p <- btrecase.plot(dt[Signif !="None",],
                       x1=paste(cols, rbias[1], sep="."),
                       x2= paste(cols, rbias[2], sep="."),
                       xl='eQTL-effect (no correction)', yl='eQTL-effect\n(correction)',
                       col=c("Without", "With"),axis.title=12, axis.text=10,
                       legend.title=12, legend.text=10,
                       legend.symbol=4, point.size=3 ,
                       title=paste("Reference panel bias correction\n", j), title.size=12,
                       colors=colors.p) +
        annotation_custom(tableGrob(tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=k[1], xmax=k[2], ymin=k[3], ymax=k[4])
    return(p)
},
i=r.bias,
j=c("Observed-GT", "Hidden-GT"),
k=list(c(-1.0, -.4, .5, 1.2), c(-.4,-.2, .25, .5)),
SIMPLIFY=F)


## lm correction vs no correction (y ~ x no intercept)

 mapply(function(i,j,k) {
    dt <- i[abs(gene.dist) <=d & Signif != "None",]
    cr <- lm(dt$log2_aFC_mean.rbias ~ dt$log2_aFC_mean.norefbias-1)
    
    return(summary(cr)$coefficients)    
},
i=r.bias,
SIMPLIFY=F)

####################################################################


## select rows that m.dna arent all 0

keep <- which(rowSums(m.dna)>0)

m.dna.snp <- rowSums(m.dna[keep,])

n.dna.snp <- unlist(lapply(names(keep), function(j) sum(unlist(lapply(prefix, function(i) c(unlist(dna[id == j & ( get(paste0(i,"GT")) == "0|1")  ,paste0(i,"ASE"), with=F] ),
                                                                                            m.dna[j,paste0(i,"ASE") ] - unlist(dna[id == j & get(paste0(i,"GT")) == "1|0" ,paste0(i,"ASE"), with=F] ) ))))))


AI.hets <- merge(AI, data.table(id=names(m.dna.snp), raw.AI=n.dna.snp/m.dna.snp, depth=m.dna.snp), by="id")


## plot fsnps used for inference
AI.hets.p <- ggplot(AI.hets[ Inference.DNA=='Yes',], aes( log(depth),qlogis(raw.AI) )) +
    geom_point(shape=1, color="#5ab4ac") +
    geom_vline(xintercept=log(100), linetype="dashed", color="gray") +
    geom_hline(yintercept=0, linetype="dashed") +
    ylab("AI (raw estimate)") +
    xlab("Depth (log)") +
    scale_y_continuous(breaks=c(-3,0,3,6), limits=c(-4,7)) +
    scale_x_continuous(breaks=seq(6,12,2), limits=c(1, 14)) +
    #scale_x_continuous(breaks=seq(6,12,2))+
    #scale_y_continuous(breaks=c(-3,0,3,6), limits=c(-3,4))+
    geom_rug(col=rgb(.3,.5,.5,alpha=.2))

AI.hets.p.rna <- ggplot(AI.hets[ Inference.RNA=='Yes',], aes( log(depth),qlogis(raw.AI) )) +
    geom_point(shape=1, color="#5ab4ac") +
    geom_vline(xintercept=log(100), linetype="dashed", color="gray") +
    geom_hline(yintercept=0, linetype="dashed") +
    ylab("AI(raw estimate)") +
    xlab("Depth (log)") +
    scale_y_continuous(breaks=c(-3,0,3,6), limits=c(-4,7)) +
    scale_x_continuous(breaks=seq(6,12,2), limits=c(1, 14)) +
    geom_rug(col=rgb(.3,.5,.5,alpha=.2))



##########################################################################################################
##  Fig2 ##
######################
AI.est.gt <- plot_grid(AI.hets.p, EAI.dna , rbias.p[[1]],  nrow=1,  rel_widths=c(0.8, 0.8, 1))

ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/Fig2bottom.png', AI.est.gt, width = 12, height = 3)                    

######################################################################################################
## Supp Fig2 ##
######################
AI.est.ngt <- plot_grid( AI.hets.p.rna, EAI.rna,  rbias.p[[2]], labels="auto", ncol=1 )#,rel_widths=c(0.6, 1, 0.8) )


ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/SuppFig2.png', AI.est.ngt, width = 6, height = 6.8)

###################################################################################################################
## Prior control, compare rbias estimates with GT to flat prior N(0,10)
## Supp Fig 1 ##
##################################################################################################################

#####################################################
## plot distribution for priors

infp=snakemake@params[['inf_prior']]
noninf=snakemake@params[['non_infprior']]

## infp <- list(mean=c(0,0), sd=c( 0.0309, 0.3479), mix=c(0.97359164, 0.02640836))
## noninf <- list(mean=0, sd=10, mix=1)

## helper function
f <- function(prior){
    k=length(prior)/3 ## number of gaussians
    s <- seq(1,length(prior),k)
    l <- lapply(1:3, function(i) as.numeric(prior[s[i]: (s[i]+k-1)]))
    names(l) <- c("mean", "sd", "mix")
    return(l)
}

infp <- f(infp)
noninf <- f(noninf)


## mix normal for stat_function for ggplot

## mix.n <- function(x,prior){
##     ## number of components
##     n <- unique(sapply(prior, length))
##     tmp <- Reduce("+", lapply(seq_along(n), function(i) prior$mix[i]*dnorm(x, prior$mean[i], prior$sd[i])))
##     return(tmp)
## }


## mix.f <-function(x, i){prior$mix[i]*dnorm(x, prior$mean[i], prior$sd[i])}
## sdnorm  <- function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}

## infp <- mix.n(infp)
## hist(mix.n(noninf))
    
## ggplot(mix.n(infp), aes(values)) + geom_density()
## plot(mix.n(infp))

## xlim=30
## x=seq(-xlim,xlim,.01)
## colors=c("#999999", "#F0E442", "#0072B2", "#D55E00")

## ## Use same colors as btrecase.plot

## inf.p <- ggplot(data.table(x=x), aes(x)) +
##     stat_function(fun=mix.f, args=list(i=1), color=colors[3]) +
##     xlab("") +
##     ylab("Density") +
##     ggtitle("Informative prior") +
##     theme(plot.title = element_text(face="plain"))

## uninf.p <- ggplot(data.table(x=x), aes(x)) +
##     #stat_function(fun=mix.n, args=list(infp), color="blue") +
##     stat_function(fun=mix.n, args=list(noninf), color=colors[2]) +
##     xlab("") +
##     ylab("Density") +
##     ggtitle("Uninformative prior")+
##     theme(plot.title = element_text(face="plain"))


## From Chris
d0 <- function(x) 
  dnorm(x, mean=0, sd=0.03) #* 0.97
d1 <- function(x) 
    dnorm(x, mean=0, sd=0.35) #* 0.03

labs <- data.frame(x=log(c(1.05,1.3)),
                   fill=factor(c(0,1)),
                   lab=c("97% 'null'\nN(0,0.03)","3% 'non-null'\nN(0,0.35)"))

labs$y <- c(d0(labs$x[1]),d1(labs$x[2]))
df0 <- data.frame(x=c(seq(-0.15,0.15,0.001)),fill=factor(0,levels=c(0,1)))
df0$y <- d0(df0$x)
df1 <- data.frame(x=c(seq(-1.2,1.2,0.01)),fill=factor(1,levels=c(0,1)))
df1$y <- d1(df1$x)
xe <- c(1/3, 1/2, 1/1.5, 1, 1.5, 2, 3)

inf.p <- ggplot(mapping=aes(x=x,y=y)) +
  geom_polygon(aes(fill=fill),data=df0,col="grey20") +
  geom_polygon(aes(fill=fill),data=df1,col="grey20") +
  geom_label_repel(mapping=aes(label=lab,fill=fill),data=labs,nudge_y=0.2,nudge_x=0.1,size=6/.pt) +
  scale_x_continuous("Odds ratio", breaks=log(xe),labels=signif(xe,2)) +
  scale_fill_manual(values=c("0"="#99999933","1"="#00990033")) +
  ylab("Density") +
  #theme_cowplot(font_size=6) +
  theme(legend.position="none",
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        plot.title=element_text(face="plain")) +
    ggtitle("Mixture prior on eQTL effects")



##################################################################

## fpriordir <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/weakprior/GT"
fpriordir <- snakemake@params[['flatPrior']]


fprior <- comb.files(path=fpriordir, pattern="ENSG[0-9]+.*stan.summary.txt")
fprior <- fprior[Rhat < 1.1,]

rbas.prior <- merge(btrecase$GT[rbias=="rbias" & abs(gene.dist)<=d, ], fprior, by=c("Gene_id", "tag"), suffixes=c(".Informative", ".Uninformative"))

cols <- c(None="#999999", Uninformative="yellow3", Informative="#0072B2", Both= "#D55E00")

rbas.prior <- add.signif(rbas.prior, x2=paste0("null.99",  ".Informative"), x1=paste0("null.99", ".Uninformative"), col=rev(c("Informative","Uninformative")))

prior.tab <- tab2bplot(dt=rbas.prior, colors=cols)

prior.p <- btrecase.plot(dt=rbas.prior[Signif != "None"  ,],
                        x2=paste(c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', 'null.99'), "Informative", sep="."),
                        x1= paste0(c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', 'null.99'), ".Uninformative"),
                        yl='eQTL effect\n(informative prior)', xl='eQTL effect (uninformative prior)',
                        col=rev(c("Informative","Uninformative")),axis.title=12, axis.text=10,
                        legend.title=12, legend.text=10,
                        legend.symbol=4, point.size=3 ,
                        title="Effect of prior on eQTL estimates", title.size=12) +
    annotation_custom(tableGrob(prior.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(prior.tab$color, rep("black", 4)))))), xmin=-35, xmax=-29, ymin=0.5, ymax=2.5)

cor.test(x=rbas.prior[Signif=="Both" | Signif == "Informative", log2_aFC_mean.Informative],
       y=rbas.prior[Signif=="Both" | Signif == "Informative", log2_aFC_mean.Uninformative])

## plotting
prior.comb <- lapply(list(inf.p, prior.p), function(p) ggplot_gtable(ggplot_build(p)))
## make same witdths

prior.comb[[1]]$widths  <- prior.comb[[2]]$widths
tmp <- Reduce(function(a,b) a$widths  <- b$widths, prior.comb)

##top.sup1 <- plot_grid(inf.p,labels='auto')
sup1 <- plot_grid(plotlist=prior.comb, ncol=1) ##inf.p, prior.p, ncol=1, labels="auto")

ggsave('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/SuppFig1.png', sup1, width = 6, height = 6.3)


###########################################################################################################
### Prepare data for SuppTable 1 for simulations of reference panel, extracted/adapted from btrecase.R ###
########################################################################################################

## sim.dir <- "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_50"
sim.dir <- snakemake@params[['sim_dir']]
## sim.pat <- c('diff_0.05_hfq_0.1_pop', 'diff_0.05_hfq_0.1_beta')
## sim.pat2 <- c('diff_0_hfq_0.1_pop', 'diff_0_hfq_0.1_beta')
sim.pat <- snakemake@params[['sim_pattern']]
sim.pat2 <- snakemake@params[['sim_pattern2']]


dif05 <- sim.for(sim.dir, sim.pat)
dif0 <-  sim.for(sim.dir, sim.pat2)

## Print table

table.sim(dif05,cap="Effect of external panel on eQTL estimates. A population (Pop) of 50,000 haplotypes of a cis-eQTL and 3 fSNPs were simulated, with a $\\beta_{aFC}=0.4$ (Online Methods). From this population a random sample of 1000 haplotypes was extracted as was used as reference panel (RP). Samples of 100 haplotypes were also extracted from the population of haplotypes if the sum of the square difference of haplotype frequencies between the population and the sample relative to the haplotype frequency on the population was equal or higher than 0.05, for those haplotypes with frequency above 0.1 in the population. This procedure was repeated 100 times. For each of the 100 samples eQTL effects were estimated either with the full model (modelling both between-individual (BI) and ASE signals) or ASE signals only. Each model was run with either known sample haplotypes (True haplotypes) or treating phasing as latent as estimating sample haplotypes using haplotypes from the population, the reference panel or the sample itself. The table shows the mean  $\\widehat{\\beta_{aFC}}$, the proportion of times $\\beta_{aFC}$ is in the 95\\% credible interval and the proportion of times than the null value is within the 95\\% credible interval.")



### Sup Table with samples used in analysis
## get date for dowloaded samples
## get name of samples

downD <- fread(snakemake@input[["sample_down_date"]])
sampN <- fread(snakemake@input[["sample_info"]], header=F)

##downD <- fread("/mrc-bsu/scratch/ev250/EGEUV1/sample_info/sample.downdate")
## sampN <- fread("/mrc-bsu/scratch/ev250/EGEUV1/sample_info/GBR.samples", header=F)

names(downD) <- paste0("V", 1:ncol(downD))
sampN[, Bname:=basename(V2)]

sampD <- merge(sampN, downD[, .(V6,V7,V8,V9)], by.x="Bname", by.y="V9")

## select samples used for analysis

sampU <- gsub("_GT", "", grep("_GT", names(rna), value=T))

sampD <- sampD[V1 %in% sampU,]

sampD[, Downloaded_date:=paste(V6,V7,V8, sep="_")]

sampD <- sampD[,.(V1, V2, Downloaded_date)]

setnames(sampD, c("V1", "V2"), c("Name", "URL"))

setkey(sampD, Name, URL)

write.table(sampD, '/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/SuppTableX.txt', row.names=F)

