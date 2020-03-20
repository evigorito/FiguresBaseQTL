library(data.table)
library(ggplot2)

#' Creates auxiliary graphs for Figure 1
#' 
#' @param out full name to save graphs
#' 
#' @return saves file with mapping info of fSNPs to genes
#' @export
#' graphs.aux()

graphs.aux <- function(out){

    ## create plot for between individual variation 

    DT <- data.table(Counts=c(20,30,37),Genotype= c("0/0","0/1","1/1"))

    pt <- ggplot(DT, aes(x=Genotype, y=Counts)) + 
        geom_bar(stat="identity") + theme_bw() + 
        ggtitle("Between individual \nexpression") +
        ylab("Total counts (c)") + 
        theme(plot.title = element_text(size=30) , axis.title = element_text(size=24), axis.text.x = element_text(colour= "black", size = 24),axis.text.y = element_text(colour="black", size = 24))


######## within individual variation ######
    olive.col <- colors()[grep("green", colors())]
    orange.col = colors()[grep("orange", colors())]
    i=12


    AI <-data.table(Geno=sort(rep(c("0/0","0/1","1/1"),2),decreasing = F),  
                    Gene.expression=c(6,6,11,6,11,11), 
                    fill=c( rep(orange.col[6],2) , 
                           olive.col[i], 
                           orange.col[6], 
                           rep(olive.col[i],2)))

    bars <- c(6, NA, 6, NA, 11,NA)


    pw <- ggplot(AI, aes(x=Geno, y=Gene.expression, fill=fill)) +
        
        geom_bar(stat="identity", colour="black") +theme_bw() + 
        ylab("Mapped counts (m)")  + 
        xlab("Genotype") +
        theme(plot.title = element_text(size=30), axis.title = element_text(size=24), 
              axis.text.x = element_text(colour="black", size = 24), 
              axis.text.y = element_text(colour="black", size = 24)) + 
        
                                        #geom_errorbar(aes(y = bars, ymin = bars, ymax = bars)) + 
        theme(legend.position='none') + 
        scale_fill_manual(values=c(olive.col[i], orange.col[6]
                                   )) +
        annotate("text", x = c(1,2,3), y = c(8, 13, 17), 
                 label = paste(rep("pi ==",3), 
                               c(0.5, round(11/18, 1), 0.5)) , 
                 size=8, parse = TRUE) +
        ggtitle("Within individual \nexpression")


    ## save plots
    l <- list(pt,pw)
    names(l) <- unlist(lapply(c("Between", "Within"), function(i) grep(i, out, value=T)))

    for (i in seq_along(l)){
        ggsave(names(l)[[i]], l[[i]], device="png")
    }
}

    


graphs.aux(out = snakemake@output[['outf']])

