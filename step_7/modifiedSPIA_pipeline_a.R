suppressMessages(library(SPIA))
suppressMessages(library(optparse))

options(error=traceback)

option.list <- list(
    make_option(c("-b","--babelomics"), type="character", default=NULL, dest="babelomics", 
                metavar="file", help="[REQUIRED] Input file with the babelomics results."),
    make_option(c("-d","--deseq"), type="character", default=NULL, dest="deseq", 
                metavar="file", help="[REQUIRED] Input file with miRNA DESeq results"),
    make_option(c("-m","--mirna"), type="character", default=NULL, dest="mirna", 
                metavar="file", help="[REQUIRED] File with miRNA names and miRBase accessions"),
    make_option(c("-t","--tdp43"), type="character", default=NULL, dest="tdp43", 
                metavar="file", help="[REQUIRED] File with TDP-43 regulated miRNAs"),
    make_option(c("-x","--outprefix"), type="character", default="result", dest="out.prefix",
                metavar="prefix", help="Prefix for output files")
)

parser <- OptionParser(option_list=option.list)
opt <- parse_args(parser)

out.prefix <- opt$out.prefix
out.dir <- opt$out.dir
# load data
babelomics.table <- read.table(opt$babelomics, header=T, sep="\t", stringsAsFactors=F, comment.char="")
deseq_mirna.table <- read.table(opt$deseq,header=T,sep="\t", stringsAsFactors=F)
mirname2acc <- read.table(opt$mirna, header=T, sep="\t", stringsAsFactors=F)
mirname2acc <- mirname2acc[,c(1,2)] # this is to cover situation where I use miRBase database with extra columns
names(mirname2acc) <- c("miRBaseAccession","miRNA_name")
tdp43.mirnas <- read.table(opt$tdp43, sep="\t")
candidate.mirnas <- as.character(tdp43.mirnas[,1])

# combine the tables into a clean dataset
new.deseq.table <- merge(mirname2acc, deseq_mirna.table, by="miRBaseAccession", all.y=T)
new.overall.table <- merge(new.deseq.table, babelomics.table, by.x="miRBaseAccession", by.y="X.term", all=T)

final.columns <- c("miRBaseAccession","miRNA_name","log2FoldChange","odds_ratio_log","term_size","list2_positives","list2_percentage","list2_positive_ids","padj","adj_pvalue")
final.all.table <- new.overall.table[!is.na(new.overall.table$term_size),final.columns]
names(final.all.table) <- c("miRBaseAccession","ID","log2FoldChange","odds_ratio_log","numPredictedTargets","numDEgenes","percentDEgenes","DEgeneList","pNDE","pFatiscan")

final.all.table[is.na(final.all.table$pNDE),"pNDE"] <- 1
# add the column where the p-values are combined
final.all.table$pG <- combfunc(final.all.table$pNDE, final.all.table$pFatiscan, combine="norminv")
final.all.table$pGFdr <- p.adjust(final.all.table$pG, method="BH")
final.all.table$pGFWER <- p.adjust(final.all.table$pG, method="bonferroni")
final.all.table <- final.all.table[order(final.all.table$pG),]

final.pos.gene.table <- subset(final.all.table, odds_ratio_log >= 0)
final.neg.gene.table <- subset(final.all.table, odds_ratio_log < 0)

# write out the table with the pG included
write.out.spia.tables <- function(table.prefix, final.table) {
    final.pos.table <- subset(final.table, log2FoldChange >= 0)
    final.neg.table <- subset(final.table, log2FoldChange < 0)
    write.table(final.table, paste(table.prefix, "allResults.txt", sep="_"), quote=F, sep="\t", row.names=F)
    write.table(final.table[final.table$pGFdr <= 0.05,], paste(table.prefix, "allSigResultsFDR0.05.txt", sep="_"), quote=F, sep="\t", row.names=F)
    write.table(final.table[which(final.table$ID %in% candidate.mirnas),], paste(table.prefix, "allTDP43MirnaResults.txt", sep="_"), quote=F, sep="\t", row.names=F)

    write.table(final.pos.table, paste(table.prefix, "upMirnaAllResults.txt", sep="_"), quote=F, sep="\t", row.names=F)
    write.table(final.pos.table[final.pos.table$pGFdr <= 0.05,], paste(table.prefix, "upMirnaSigResultsFDR0.05.txt", sep="_"), quote=F, sep="\t", row.names=F)
    write.table(final.pos.table[which(final.pos.table$ID %in% candidate.mirnas),], paste(table.prefix, "upTDP43MirnaResults.txt", sep="_"), quote=F, sep="\t", row.names=F)

    write.table(final.neg.table, paste(table.prefix, "downMirnaAllResults.txt", sep="_"), quote=F, sep="\t", row.names=F)
    write.table(final.neg.table[final.neg.table$pGFdr <= 0.05,], paste(table.prefix, "downMirnaSigResultsFDR0.05.txt", sep="_"), quote=F, sep="\t", row.names=F)
    write.table(final.neg.table[which(final.neg.table$ID %in% candidate.mirnas),], paste(table.prefix, "downTDP43MirnaResults.txt", sep="_"), quote=F, sep="\t", row.names=F)
    invisible()
}

upgene.prefix <- paste(out.prefix,"upGenes",sep="_") 
write.out.spia.tables(upgene.prefix, final.pos.gene.table)
downgene.prefix <- paste(out.prefix,"downGenes",sep="_")
write.out.spia.tables(downgene.prefix,final.neg.gene.table)
