## Load R packages and functions:


stan.dir <- snakemake@params[['stan_dir']]
skin <- snakemake@params[['skin']]
rbias <- snakemake@params[['rbias']]
d <- as.numeric(snakemake@params[['d']])
info.max <- as.numeric(snakemake@iparams[['infoMax']])
gene_coord <- snakemake@input[['gene_coord']]
gene_names <- snakemake@params[['genes2follow']]
y=list(snakemake@params[['y_G1']], snakemake@params[['y_G2']])
j=list(snakemake@params[['j_G1']], snakemake@params[['j_G2']])
r=list(snakemake@params[['r_G1']], snakemake@params[['r_G2']])
w=as.numeric(snakemake@params[['w']])
tagsdir=snakemake@params[['tags_dir']]
out=unlist(snakemake@output)
Gtex.dir=snakemake@params[['gtex_dir']]
sig.exp=snakemake@input[['gtexSigExp']]
sig.noexp=snakemake@input[['gtexSigNoexp']]
rsfiles=snakemake@input[['var']]
cor <- readRDS(snakemake@input[['cor']])
f2 <- snakemake@params[['stan1M2T_dir']] ## same fSNPs as when running individual models
f.sum2 <- list.files(f2, "2Tissues.*.txt", full.names=T)
drg=data.table(read.xlsx2(snakemake@input[['drg']],
                          sheetIndex=1 ,
                          colClasses=snakemake@params[['colclass']]))
gwas=snakemake@input[['gwas']]
gwasr2=snakemake@input[['gwasr2']]

if(!exists("tagsdir")){
    tagsdir <- NULL
}

if(!exists("infoMax")){
    info.max <- NULL
}



library(data.table)
library(ggforce)
library(grid)
library(pBrackets)
library(xlsx)


source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")
source("/home/ev250/Bayesian_inf/trecase/Functions/various.R")
source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R", verbose=FALSE)

#' Figures for eQTL no genotypes paper: psoriasis vs normal skin
#'
#' @param stan.dir path to dir with stan output or running model with and without ref bias correction
#' @param rbias prefix to identify BaseQTL runs with reference panel bias correction
#' @param skin character with names used to identify skin type
#' @param gene_coord full name to data table  with gene_id, gene start and end columns
#' @param d distance to use for cis window, defaults to 10^5
#' @param info.max upper limit for info score to consider associations, defaults to NULL
#' @param Gtex.dir path to dir with with significant associations for Gtex skin
#' @param sig.exp full name with significant associations for exposed skin in gtex
#' @param sig.noexp full names with signficant associations for non exposed skin in gtex
#' @param gene_names character vector with the names of genes to focus
#' @param y list with y coordinates for gene labels, each element for each gene
#' @param j list with just labels for each gene
#' @param r list with yaxis coord for rectangle
#' @param w numeric to adjust the size of the window for making rectangle for gene, defaults to 2000
#' @param tagsdir path to dir with tags for simple model, NULL if the same as stan.dir
#' @param cor full name to file with r/r2 summaries for snps
#' @param rsfiles full name for files with variant information to extract rsid as formatted as ensembl ftp.ensembl.org/pub/grch37/current/variation/gvf/homo_sapiens/homo_sapiens-chr22.gvf.gz for appropiate built and chromosome, for each gene to plot
#' @param f.sum2 names of files with stan output with 2T into joint model
#' @param drg file with DRG genes in pso vs normal skin, downloaded from
#' @param gwas file with psoriasis gwas hits used to select genes to run
#' @param gwasr2 file with square correlation of gwas hit with rsnps run in model
#' @param out full file names to save figures
#' @keywords figures eQTL no-genotypes
#' @export
#' @return save ggplots to files
#' 
#' pso_nor1()

pso_nor1 <- function(stan.dir, rbias, skin, gene_coord, d=10^5, info.max=NULL, Gtex.dir, sig.exp, sig.noexp, gene_names, y,j, r,  w=2000, tagsdir=NULL, cor,  rsfiles, f.sum2, drg, gwas, gwasr2, out){

    ## get summaries by skin with ref bias correction

    geneStEnd <- fread(gene_coord)
    Skin <- unlist(lapply(skin, function(i) gsub("_", " ", paste0(toupper(substr(i,1,1)), substr(i, 2, nchar(i))))))
   
    ## remove bad runs: Rhat>1.01
     ## make PEP 0 to 1/4001 and take log for plotting
    
    all <- lapply(skin, function(j) {
            dt <- comb.files(stan.dir,pattern=paste0("^", rbias,'\\.ENSG[0-9]+.', j, '.noGT.stan.summary.txt'))
            dt <- dt[Rhat <1.01,]
            dt <- dt[info<=info.max,]
            dt[, skin:=j]
            dt <- add.name(sum=dt)
            dt <- gene.d(dt, geneStEnd[, .(gene_id, start,end,chrom)])
            dt <- dt[PEP==0, PEP:=1/4001][, minuslog10PEP:=-log(PEP,10)]
            dt <- merge(dt, geneStEnd[, .(gene_id, chrom)], by.x="Gene_id", by.y="gene_id")
            return(dt)
        })
    names(all) <- skin

    all.l <- rbindlist(all)
    ## common assoc
    all.w <- Reduce(function(a,b){
        dt <- merge(a,b, by=c("Gene_id", "chrom", "tag", "tag.EAF", "gene.dist"), suffixes=paste0(".",skin))
        dt <- add.signif(dt, paste("null.99", skin[1], sep="."), paste("null.99", skin[2], sep="."), Skin)
        return(dt)
       },
        all)
   
   
    ## Plot
    pso.norm.tab <- tab2bplot(all.w[abs(gene.dist) <=d,], colors=setNames(c("#999999","yellow3", "#0072B2","#D55E00"), c("None", Skin, "Both")))
    cols <- c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', 'null.99')
    
    pso.norm.p <- btrecase.plot(all.w[abs(gene.dist)<=d & Signif != "None",],
                                x1=paste(cols, skin[1], sep="."),
                                x2= paste(cols, skin[2], sep="."),
                                xl=paste0('eQTL-effect (',Skin[1],')'), yl=paste0('eQTL-effect (',Skin[2],')'),
                                col=Skin,
                                axis.title=12, axis.text=10,
                                legend.title=12, legend.text=10,
                                legend.symbol=4, point.size=3 ,
                                title="Psoriasis vs normal skin", title.size=12) +
        annotation_custom(tableGrob(pso.norm.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(pso.norm.tab$color, rep("black", 4)))))), xmin=-1.8, xmax=-.3, ymin=0.7, ymax=1.8)

    ## correlation of effects
    print(cor.test(all.w[abs(gene.dist) <= d & Signif != "None" ,get(paste0("log2_aFC_mean.", skin[1]))],
              all.w[abs(gene.dist) <= d & Signif != "None"  , get(paste0("log2_aFC_mean.", skin[2])) ]))

    ##  genes tested in each model and in common, data to add to manuscript text

    each <- lapply(all, function(i) unique(i[abs(gene.dist)<=d ,Gene_id]))
    ## number_of_genes_each_model
    lapply(each, length)
    ## number of common genes
    length(Reduce(intersect, each))
    ## number of ind genes per skin
    lapply(each, function(i) length(i) - length(Reduce(intersect, each)))
    
##########################################################################################################
#' ## Compare with Gtex data for skin tissues
#################################################

    
    gtex.skin <- list.files(Gtex.dir, full.names=T)

    gtexList <- mclapply(gtex.skin, function(i) {
        dt=fread(cmd=paste0("zcat ", i, " | grep -E '", paste(unique(all.l$Gene_id), collapse="|"), "' "),header=F)
        names(dt)  <- names(fread(cmd=paste0("zcat ", i, " | head -1 " )))
        return(dt)
    })

    ##names(gtexList) <- c("no_sun", "sun")
    ## Read signif associations in gtex skin 

    sigGtex <- lapply(list(sig.noexp, sig.exp), function(i) {
        dt <- fread(cmd=paste0("zcat ", i, " | grep -E '", paste(unique(all.l$Gene_id), collapse="|"), "' "),header=F)
        names(dt)  <- names(fread(cmd=paste0("zcat ", sig.exp, " | head -1 " )))
        dt[, null:="no"]
    })

    ## merge sig associations with 'all gtex' and select strongest association when combining exp and non exp skin

    all.gtex <- rbindlist(mapply(function(a,b) {
        dt <- merge(a,b, by=names(a) ,all =T)
        dt[is.na(null), null:="yes"]           
        ## format to match sum Gene_id and tag col
        dt[, Gene_id:=gsub("\\..*","", gene_id)]
        ##dt[, Chrom:=gsub( "_.*", "", variant_id)]
        dt[, SNP:=sub("^[0-9]*_", "", variant_id)]
        dt[, SNP:=sub("_b37", "", SNP)]
        dt[, SNP:=gsub("_", ":", SNP)]
        return(dt)
    }, a=gtexList,
    b=sigGtex,
    SIMPLIFY=F))

    setkey(all.gtex, variant_id, gene_id, min_pval_nominal)
    all.gtex <- all.gtex[,.SD[1], by=.(variant_id, gene_id)]

    ## Prepapre table with Gene_name, Gene_id, tag, tag.EAF, log2mean-CI for normal and psoriasis, slope and pval when available in Gtex
   
    skin.wide.all <- Reduce(function(a,b) merge(a[,.(Gene_name, Gene_id, chrom, tag, tag.EAF, log2_aFC_mean, `log2_aFC_0.5%`, `log2_aFC_99.5%`, null.99, PEP, minuslog10PEP)],
                                                b[,.(Gene_name, Gene_id, chrom, tag, tag.EAF, log2_aFC_mean, `log2_aFC_0.5%`, `log2_aFC_99.5%`, null.99,  PEP, minuslog10PEP)],
                                                by=c("chrom", "Gene_name", "Gene_id", "tag", "tag.EAF"),
                                                all=T,
                                                suffixes=paste0(".", names(all))),
                            x=all)
    ## Add gtex info,

    skin.wide.all.gtex <- merge(skin.wide.all, all.gtex[, .(Gene_id, slope, pval_nominal, null, SNP)], by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"), all.x=T)
    reorder <- c("chrom", "Gene_id","Gene_name")
    setcolorder(skin.wide.all.gtex, c(reorder, names(skin.wide.all.gtex)[!names(skin.wide.all.gtex) %in%  reorder]))
    setnames(skin.wide.all.gtex, c("slope", "pval_nominal"), paste0(c("slope", "pval_nominal"), ".GTEx_skin"))

    ## Add whether gene is DRG and r2 with gwas hit
    gwas <- fread(gwas)
    skin.wide.all.gtex <- merge(skin.wide.all.gtex, drg, by.x="Gene_name", by.y="Gene.Symbol", all.x=T)

    gr2 <- fread(gwasr2)

    skin.wide.all.gtex <- merge(skin.wide.all.gtex, gr2[, CHROM:=as.character(CHROM)], by.x=c("Gene_id", "tag", "chrom"), by.y= c("Gene_id", "tag", "CHROM"), all.x=T)
    
############################################################
################# Save SuppTable 3 #########################
############################################################

    write.table(x= skin.wide.all.gtex, file=grep("SuppTable3.csv$", out, value=T), row.names=F, sep=",")
    
    ## Compare associations in normal & psoriasis with reference panel bias correction with Gtex, using d

    skin.gtex <- lapply(all, function(i) {
        dt <- merge(i[abs(gene.dist) <=d, ], all.gtex, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
        dt <- add.signif(dt, "null.99", "null", col=c("BaseQTL", "Gtex-skin"))
        return(dt)
    })
 

    tab <- lapply(skin.gtex, function(i) {
        tab2bplot(i, var="Signif", colors=setNames(c("#999999","yellow3","#D55E00","#0072B2"), c("None", "BaseQTL", "Both", "Gtex-skin")))
    })

#################
## Venn diagram
################

######### make distinction between genes not run in one model and sig/not sig with the other.
    ## ENSG00000196805 has 1 fsnp with missing value in normal skin ###

    ## look at intersect and per skin sig genes

    common <- sort(Reduce(intersect, lapply(all, function(i) unique(i[null.99=="no" & abs(gene.dist)<=d, Gene_name]))))
    nor <- sort(Reduce(setdiff, lapply(all, function(a) unique(a[null.99=='no' & abs(gene.dist)<=d, Gene_name]))))                            
    pso <-  sort(Reduce(setdiff, lapply(rev(all), function(a) unique(a[ null.99=='no' & abs(gene.dist)<=d, Gene_name]))))

    ## genes run in both skins but significant in either
    tot.both <-  sort(Reduce(intersect, lapply(all, function(i) unique(i[ abs(gene.dist)<=d, Gene_id]))))
    nor.both <- sort(Reduce(setdiff, lapply(all, function(a) unique(a[Gene_id %in% tot.both & null.99=='no' & abs(gene.dist)<=d, Gene_name]))))
    pso.both <-  sort(Reduce(setdiff, lapply(rev(all), function(a) unique(a[Gene_id %in% tot.both & null.99=='no' & abs(gene.dist)<=d, Gene_name]))))

    pso.only <- pso[!pso %in% pso.both]

    ## indicate pso.only with * for Venn diagram to write in figure legend
    pso[pso %in% pso.only] <- paste0(pso.only, "*")

    ## After examining the plots some associations are not very convincing, remove

    ## Get genes sig by BaseQTL and with associations in gtex

    all.tested <- mapply(function(i,j) unique(i[Gene_name %in% unique(j[['Gene_name']]) & null.99=="no"& abs(gene.dist)<=d, Gene_name]),
           i=all,
           j=skin.gtex,
           SIMPLIFY=F)

    ## Get gene-SNP associations significant by BaseQTL and Gtex
    sig.both <- lapply(skin.gtex, function(i) unique(i[Signif== "Both", Gene_name]))

    sig.skin <- mapply(function(i,j) unique(i[Signif== "BaseQTL" & !Gene_name %in% j , Gene_name]),
                       i=skin.gtex,
                       j=sig.both,
                       SIMPLIFY=F)
     

    ## Get genes sig in BaseQTL but not tested in Gtex using the same associations

    sig.nogtex <- mapply(function(i,j) unique(i[null.99=="no" & abs(gene.dist)<=d & !Gene_name %in% unique(j[['Gene_name']]) , Gene_name]),
                         i=all,
                         j=skin.gtex,
                         SIMPLIFY=F)

    df.venn <- data.table(x=c(-.5,.2), y=rep(0,2), labels=Skin)
    cols.d <- c( "#0072B2" ,"#F0E442")

    ## for braces
    bracketsGrob <- function(...){
        l <- list(...)
        e <- new.env()
        e$l <- l
        grid:::recordGrob(  {
            do.call(grid.brackets, l)
        }, e)
    }

    ## prepare dts for normal, common or pso only genes, color black if no sig in gtex-skin, blue if no tested and pink if sig
    dts <- lapply(list(nor, common, pso), function(i) {
        tmp <- data.table(genes= i)
        tmp[, color:="black"][genes %in% unique(unlist(sig.both)), color:="maroon3"][genes %in% unique(unlist(sig.nogtex)), color:="blue"]
        setkey(tmp, color)
        return(tmp)
    })

    names(dts) <- c("nor", "common", "pso")
            
    
    v <- ggplot(df.venn) +
        geom_circle(aes(x0=x,y0=y,r=c(1,1.3), fill=labels, color=labels), alpha=0.15, size=1)  +
        geom_rect(mapping=aes(xmin=-3.5,xmax=3.5, ymin=-2.6, ymax=2, fill="Total genes"), alpha=0.01, color="black") +
        coord_fixed() +
        theme_void() +
        scale_fill_manual(values = c(cols.d, "white") ) +
        scale_colour_manual(values = c(cols.d, "black"), guide = FALSE) +         
        labs(fill = NULL) +
        annotate("text", x= c(-1.3, -0.6, 1, 0),
                 y=c(.1,0.1,.1,1.7),
                 label=c(sapply(dts, nrow) ,paste("Total genes:",length(unique(all.l$Gene_id)))),
                 size=6, colour=c("#0072B2","green4","orange2", "black")) +
        annotate("text", x=rep(2, nrow(dts[['pso']])), y=rev(seq(from=-1.05, by=.23, length.out=nrow(dts[['pso']]))),
                 label=dts[['pso']]$genes,
                 size=3,
                 hjust=0,
                 colour=dts[['pso']]$color) +
        annotate("text", x=rep(-2, nrow(dts[['nor']])), y=rev(seq(from=-.38, by=.23, length.out=nrow(dts[['nor']]))),
                 label=dts[['nor']]$genes,
                 size=3,
                 hjust=1,
                 colour=dts[['nor']]$color) +
        annotate("text", x=rep(seq(from=-2.2, length.out=4, by=.9),ceiling(nrow(dts[['common']])/4))[1:nrow(dts[['common']])], ## make suffcient number of rows of 4 elements (ceiling()) but then select the right number of genes 1:nrow() )
                 y=sort(rep(seq(from=-2.3, by=.23, length.out=3),ceiling(nrow(dts[['common']])/3)), decreasing=T)[1:nrow(dts[['common']])],
                 label=dts[['common']]$genes,
                 size=3,
                 hjust=0,  
                 colour=dts[['common']]$color)

    ## add braces and arrow:
    
    b1 <- bracketsGrob(0.27, .67, 0.27, 0.47, h=0.02, lwd=2, col="#0072B2")
    b2 <- bracketsGrob(0.74, .34, 0.74, 0.81, h=0.02, lwd=2, col="orange2")
    b3 <- bracketsGrob(0.22, .25, 0.62, 0.25, h=0.02, lwd=2, col="green4")

    v <- v + annotation_custom(b1) +
        annotation_custom(b2) +
        annotation_custom(b3) +
        annotate("segment", x=-0.62, xend=-.6, y=-0.2, yend=-1.3, colour="green4", arrow=arrow(length = unit(0.03, "npc")))


    ggsave("/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/vennIndTissues.png", v, width=7, height=7)


    ## look at sig genes by gwas proximity

    gwas.prox <- skin.wide.all.gtex[Gene_name %in% c(common,nor, pso) & !is.na(rs_id) & (null.99.Psoriasis_skin == "no" | null.99.normal_skin =="no"), .(proxyR2=max(r2TagGwas)), .(Gene_name, null.99.Psoriasis_skin, null.99.normal_skin, rs_id,Alleles, POS) ]

    ## add rsid to eQTL and look in lare studies other tissues (blood, etc)

    

#########################################################################################################
######################## individual genes #############################################################

    #### IND MODEL ###
    ## Add var Skin for nice labels in plot
    all.l[ ,Skin:= gsub("_", " ", paste0(toupper(substr(skin,1,1)), substr(skin, 2, nchar(skin))))]

    ## select relevant genes
    names(gene_names) <- sapply(gene_names, function(i) unique(all.l[Gene_name == i, Gene_id]))

    ind.genes <- all.l[gene.dist<=d & Gene_id %in%   names(gene_names) ,]
    
   

#### JOINT MODEL
    sum2 <- rbindlist(lapply(grep("summary", f.sum2,value=T), fread))

    ## make sure Rhat >1.1 in 2T 1model, select rows with Rhat for all parameters < 1.1
    keep2 <- apply(sum2[, grep("Rhat", names(sum2)), with=F], 1, max) < 1.1
    sum2 <- sum2[keep2,]

    ## make significance columns based on coefficient estimates and select cis-SNPs by distance to gene
    ##sum2 <- add.sig(sum2, type="Signif", cols=c("Signif.bp", "Signif.bn"), lab=c("bp", "bn"), newCol="Signif.bp.bn")
    sum2 <- add.name(sum2)
    sum2 <- gene.d(sum2, geneStEnd[, .(gene_id, start,end,chrom)])
    sum2 <- merge(sum2, geneStEnd[, .(gene_id, chrom)], by.x="Gene_id", by.y="gene_id")
    sum2 <- sum2[abs(gene.dist)<=d,]

    ## make sum2 long to ease plotting, select ba and bd

    param <- c("ba", "bd")
    names.param <- lapply(param, function(i) grep(i, names(sum2), value=T))
    com <- names(sum2)[!names(sum2) %in% c(unlist(names.param))]

    sum2l <- rbindlist(mapply(function(a,d,e) {
        dt <- sum2[, c(com, a), with=F]
        setnames(dt,  c(com, a), c(com,d))
        ## Sort associations by PEP, if PEP==0, add 1 count 1/4001 (total is the number of posterior samples plus 1)
        dt[PEP==0, PEP:=1/4001]
        ## add col -log10(PEP)
        dt[, minuslog10PEP:=-log(PEP,10)]
        ##dt[, skin:=c]
        dt[,Param:=e]
        ##dt[, null.99:="no"][Signif=="no", null.99:="yes"]
    },
    a=names.param,
    ##b=infocol,
    ##c=rev(skin),
    e=param,
    MoreArgs=list(d=c(gsub("\\.ba", "", names.param[[1]]))), ## new names
    SIMPLIFY=F))

    

    ## ease for formatting, I will select top snp based on ba.d
    ## sum2l[,null.99:="yes"][Signif.ba=="yes", null.99:="no"]
    ## sum2l[, log2_aFC_da:=0][log2_aFC_mean.ba >0 & Signif.ba=="yes", log2_aFC_da:=`log2_aFC_0.5%.ba`][log2_aFC_mean.ba <0 & Signif.ba=="yes", log2_aFC_da:=abs(`log2_aFC_99.5%.ba`)]

    


    ### FORMAT Ind and Joint models for plotting
    cor <- readRDS(cor)
    cor <- cor[names(gene_names)]
    
    
    l2plot <- mapply(for.plot,
                     sum=list(ind.genes, sum2l[Gene_id %in% names(gene_names),]),
                     sig.col=c("null.99", "Signif"),
                     var=c("skin", "Param"),
                 MoreArgs=list(
                     rs=rsfiles,
                     r=cor ,
                     col="minuslog10PEP",
                     r.col="r"),
                 SIMPLIFY=F)
    
    names(l2plot) <- c("ind", "joint")

  ### Plotting
    
    plots1t <- lapply(names(gene_names), function(i) {
        sum=l2plot[['ind']][Gene_id == i,]
        gene_plot_skin(sum99=sum,
                       gene=i,
                       pos="tag",
                       geneStEnd=NULL, 
                       title="gene",
                       yaxis="-log10(PEP)",
                       null='null.99',
                       gwas=NULL,
                       var="Skin",                                                                          
                       ci=NULL,
                       yvar="minuslog10PEP",
                       colvar= "r.cat",
                       colors=setNames(c("blue","yellow","green", "orange", "red"),levels(sum$r.cat)),
                       shapevar="info.cat",
                       shape=setNames(c(0,1,2,5),levels(sum$info.cat)),
                       sizevar="rsid",
                       size=setNames(ifelse(levels(as.factor(sum$rsid)) == "", 1.5, 3), levels(as.factor(sum$rsid))),
                       size.leg=FALSE,
                       d=d,
                       geneCoord=geneStEnd) +geom_hline(yintercept = -log(0.005, 10), linetype="dashed")
        })
        names(plots1t) <- names(gene_names)
        

  

    plots2M <-lapply(names(gene_names), function(i) {
        sum <- l2plot[['joint']][Gene_id == i ,]
        gene_plot_skin(sum99=sum,
                       gene=i,
                       pos="tag",
                       geneStEnd=NULL, 
                       title="gene",
                       yaxis="-log10(PEP)",
                       null='Signif',
                       gwas=NULL,
                       var="Param",                                                                          
                       ci=NULL,
                       yvar="minuslog10PEP",
                       colvar= "r.cat",
                       colors=setNames(c("blue","yellow","green", "orange", "red"),levels(sum$r.cat)),
                       shapevar="info.cat",
                       shape=setNames(c(0,1,2,5),levels(sum$info.cat)),
                       sizevar="rsid",
                       size=setNames(ifelse(levels(as.factor(sum$rsid)) == "", 1.5, 3), levels(as.factor(sum$rsid))),
                       size.leg=FALSE,
                       d=d,
                       geneCoord=geneStEnd)+geom_hline(yintercept = -log(0.005,10), linetype="dashed")
    })
    names(plots2M) <- names(gene_names)

    
    bot <- 1.4
    top <- 2.4
    ylow <- 0.2
  

    ## Combine indpendent with joint model

    ##Need to make a list of nested list with each element a plot from plots and a plot of plots1Mr2 for the same gene, remove legend from plots

    p.2M.r <- mapply(function(a,b) list(a + theme(legend.position="non"),b + theme(legend.position="non")),
                      a=plots1t,
                      b=plots2M,
                      SIMPLIFY=F)

    g.2M.r <- mapply(mix.ggbio,
                  plot= p.2M.r ,
                  gene=gene_names,
                  x=names(gene_names),
                  y.text=list(c(bot, top),
                              rep(bot,2)),
                  y.rect=list( c(bot+0.2, top-0.1),                             
                               c(ylow+.5, bot-0.1)
                              ), 
                  just=list(rep(0,2),
                            rep(1,2)),
                  MoreArgs=list(geneStEnd=geneStEnd,
                                d=d,          
                                w=w,
                                yaxis=expression('-log' [10] * '(PEP)'),
                                relh=c(1,1,0.6)),
                  SIMPLIFY=F)

  

   ## get legend to share between the 2 plots
                                        #extract legend
                                        #https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
    g_legend<-function(a.gplot){
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)}

    myleg <- g_legend(plots1t[[1]])
    ##myleg2 <- g_legend(plots2[[1]])

    top.p <- plot_grid(pso.norm.p, v, nrow=1, labels="auto")
    bot.p <- plot_grid(g.2M.r[[1]] , myleg, g.2M.r[[2]], nrow=1, rel_widths=c(1,.15,1), labels=c("c", "", "d"))
    p <- plot_grid(top.p, bot.p, nrow=2, rel_heights=c(0.6,1))

    ## bot.p2 <- plot_grid(g.2M.eff[[1]] , myleg2, g.2M.eff[[2]], nrow=1, rel_widths=c(1,.15,1), labels=c("c", "", "d"))
    ## p2 <- plot_grid(top.p, bot.p2, nrow=2, rel_heights=c(0.8,1))

    ## save r2 and eff until we decide which one or alternative to plot
    ggsave(grep("Fig5r.png$", out, value=T) , p, width=15.6, height=11.5)
    ## ggsave(grep("Fig5eff.png$", out, value=T) , p2, width=15.6, height=11.5)

    ## Sankey diagram

    bd.s <- unique(sum2l[Param=="bd" & Signif=="yes", Gene_name])
    ba.s <- unique(sum2l[Param=="ba" & Signif=="yes", Gene_name])
    common.j <- ba.s[!ba.s %in% bd.s]
    ## significant by skin, some genes have sig bd for both skins on different associations
    skin.s <- lapply(c("bn", "bp"), function(i) unique(sum2l[get(paste0("Signif.", i))=="yes" & Gene_name %in% bd.s, Gene_name]))
    names(skin.s) <- skin
    ## get gene names for common nor and pso
    j.dts <- list(nor=skin.s[[skin[1]]], common=common.j, pso=skin.s[[skin[2]]])

    ## get number of genes from source (ind model) to target (2M)
    source.target <- lapply(j.dts, function(i) sapply(dts, function(j) sum(j$genes %in% i)))

    ## create a connection data table
    links <- data.table(source=unlist(apply(names(source.target), function(i) rep(i, sum(source.target[[i]]))))

    

}


## gene_names <- c("ERAP1", "PPIF")
## skin <- c( "normal_skin","Psoriasis_skin" )
## rbias <- "refbias"
## d <- 10^5
## gene_coord <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt"
## pso99 <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/results/normal_pso_ci99_summary.txt"
## stan.dir <- "/mrc-bsu/scratch/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna"
## y = list(c(1.4, 2.4), c(1.4,1.6,1.4,2.4, 1.4))
## j <- list(rep(0,2), rep(0,5))
## r <- list(c(1.7, 2.3), c(1.7, 2.3))

## w <- 3000
## Gtex.dir <-"/mrc-bsu/scratch/ev250/psoriasis/eQTL_gtex"
## sig.exp <- '/mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Skin_Sun_Exposed_Lower_leg.v7.signif_variant_gene_pairs.txt.gz'
## sig.noexp <- '/mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Skin_Not_Sun_Exposed_Suprapubic.v7.signif_variant_gene_pairs.txt.gz'
## tagsdir=NULL
## out <- c('/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/Fig5r.png', '/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/Fig5eff.png', '/mrc-bsu/scratch/ev250/bayesian_trecase/Paper/Figures/SuppTable3.csv')
## rsfiles <- c("/mrc-bsu/scratch/ev250/reference_genome/built37/variation_homo_sapiens-chr5.gvf.gz", "/mrc-bsu/scratch/ev250/reference_genome/built37/variation_homo_sapiens-chr10.gvf.gz")
## cor <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/objects/fisher001EURr.rds"
## f.sum2 <- list.files("/mrc-bsu/scratch/ev250/psoriasis/refbias/Btrecase/2Tissues/jointInfo_Fish001_hetInd", "2Tissues.*.txt", full.names=T)
## info.max = 1.1
## drg=data.table(read.xlsx2('/mrc-bsu/scratch/ev250/psoriasis/pso_gwas/DRG_JID2014.xls' ,
##                           sheetIndex=1 ,
## colClasses=c("character", rep("numeric", 7))))
## gwas <- '/mrc-bsu/scratch/ev250/psoriasis/pso_gwas/SupTable2.csv'

## gwasr2 <- "/mrc-bsu/scratch/ev250/psoriasis/refbias/Btrecase/results/r2_spike_gwas.txt"


pso_nor1(stan.dir, rbias, skin, gene_coord, d, info.max, Gtex.dir, sig.exp, sig.noexp, gene_names, y, j, r, w, tagsdir, cor, rsfiles, f.sum2, drg, gwas, gwasr2,  out)


  ## plots2 <- lapply(names(gene_names), function(i) gene_plot_skin(sum99=indr[gene.dist<=d,],
    ##                                                               gene=i,
    ##                                                               pos="tag",
    ##                                                               geneStEnd=NULL, 
    ##                                                               title="gene",
    ##                                                               null='null.99',
    ##                                                               gwas=NULL,
    ##                                                               var="Skin",                                                                          
    ##                                                               ci=c(0.5, 99.5),
    ##                                                               shapevar="rsid",
    ##                                                               shape=setNames(ifelse(levels(as.factor(indr$rsid)) == "", 1, 5), levels(as.factor(indr$rsid))),
    ##                                                               shape.leg=FALSE,
    ##                                                               sizevar="rsid",
    ##                                                               size=setNames(ifelse(levels(as.factor(indr$rsid)) == "", 1.5, 3), levels(as.factor(indr$rsid))),
    ##                                                               size.leg=FALSE,
    ##                                                               d=d,
    ##                                                               geneCoord=geneStEnd))
    ## names(plots2) <- names(gene_names)


   ## g.plots <- mapply(mix.ggbio,
    ##               plot=plots,
    ##               gene=gene_names,
    ##               x=names(plots),
    ##               y.text=list(c(bot, top),
                      
    ##                           c(rep(bot,4), top , rep(bot,2) )
    ##                           ## c(bot, bot, top, bot),
    ##                           ## c(rep(bot, 4), bot+0.1, top, bot),
    ##                           ## c(bot, top,bot,top, rep(bot,3)),
    ##                           ## c(bot, top),
                              
    ##                           ## c(bot,top),
    ##                           ## rep(bot,4),
    ##                           ## c(bot, top, bot, bot +0.2, rep(bot,5)),
    ##                           ## c(rep(bot, 5), bot+0.1, rep(bot,2)),
    ##                           ## rep(bot,8),
    ##                           ## c(bot, top, top+0.3,bot, top, stop, rep(bot,2), bot+0.2, bot, top,bot,top,bot),
    ##                           ## rep(bot,3),
    ##                           ## rep(bot, 8),
    ##                           ## c(rep(bot,3), bot+0.2, bot,top, stop,rep(bot,2)),
    ##                           ## c(bot, top, bot, top+0.1, rep(bot,2), top+0.1, bot,top, stop, bot, top, bot),
    ##                           ## c(bot+0.1, top, bot, top, bot, bot+0.1,bot),
    ##                           ## c(rep(bot,4), bot+0.1, rep(bot,3))
                              
    ##                           ),
    ##               y.rect=list( c(bot+0.2, top-0.1),
    ##                   #c(ylow, bot-0.1),
    ##                           ##c(bot+0.2, top-0.1)
    ##                            c(ylow, bot-0.1)
    ##                           ## c(ylow, bot-0.1),
    ##                           ## c(bot+0.2, top-0.1),
    ##                           ## c(bot+0.2, top-0.1),
                             
    ##                           ## c(bot+0.2, top-0.1),
    ##                           ## c(ylow, bot-0.1),
    ##                           ## c(ylow, bot-0.1),
    ##                           ## c(ylow, bot-0.1),
    ##                           ## c(ylow, bot-0.1),
    ##                           ## c(ylow, bot-0.1),
    ##                           ## c(ylow, bot-0.1),
    ##                           ## c(ylow, bot-0.1),
    ##                           ## c(ylow, bot-0.1),
    ##                           ## c(bot+0.2, top-0.1),
    ##                           ## c(ylow, bot-0.1),
    ##                           ## c(ylow, bot-0.1)
    ##                           ), 
    ##               just=list(rep(0,2),
    ##                  # rep(0, 10),
    ##                         c(rep(1,2),rep(0,5))
    ##                         ## rep(0,4),
    ##                         ## c(rep(1,2), rep(0,5)),
    ##                         ## rep(0, 7),
    ##                         ## rep(0,2),
                            
    ##                         ## rep(0,2),
    ##                         ## rep(0,4),
    ##                         ## rep(0,9),
    ##                         ## rep(0,8),
    ##                         ## rep(0,8),
    ##                         ## c(rep(0,11),1, rep(0,2)),
    ##                         ## rep(0,3),
    ##                         ## rep(0,8),
    ##                         ## c(0, 1, 0, 1, rep(0,5)),
    ##                         ## rep(0,13),
    ##                         ## c(rep(0,4), 1, rep(0,2)),
    ##                         ## rep(0,8)
    ##                         ),
    ##               MoreArgs=list(geneStEnd=geneStEnd,
    ##                             d=d,          
    ##                             w=w,
    ##                             relh=c(1,0.6),
    ##                             yaxis="Distance to null (log2)"),
    ##               SIMPLIFY=F)

    ## ## alternative with eQTL effects
    ## g.plots2 <- mapply(mix.ggbio,
    ##               plot=plots2,
    ##               gene=gene_names,
    ##               x=names(plots),
    ##               y.text=list(c(bot, top),
    ##                           c(rep(bot,4), top , rep(bot,2) )),
    ##               y.rect=list( c(bot+0.2, top-0.1),                             
    ##                            c(ylow, bot-0.1)
    ##                           ), 
    ##               just=list(rep(0,2),
    ##                         c(rep(1,2),rep(0,5))
    ##                         ),
    ##               MoreArgs=list(geneStEnd=geneStEnd,
    ##                             d=d,          
    ##                             w=w,
    ##                             relh=c(1,0.6)),
    ##               SIMPLIFY=F)

     ## g.2M.eff <- mapply(mix.ggbio,
    ##               plot= p.2M.eff ,
    ##               gene=gene_names,
    ##               x=names(gene_names),
    ##               y.text=list(c(bot, top),
    ##                           c(rep(bot,4), top , rep(bot,2) )),
    ##               y.rect=list( c(bot+0.2, top-0.1),                             
    ##                            c(ylow, bot-0.1)
    ##                           ), 
    ##               just=list(rep(0,2),
    ##                         c(rep(1,2),rep(0,5))
    ##                         ),
    ##               MoreArgs=list(geneStEnd=geneStEnd,
    ##                             d=d,          
    ##                             w=w,
    ##                             ##yaxis="Distance to null (log2)",
    ##                             relh=c(1,1,0.6)),
    ##               SIMPLIFY=F)
   
