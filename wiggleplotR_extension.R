library(wiggleplotr)
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(biomaRt)

ext_name2wiggleplotr.metadata <- function(ext_name, biomart = mart, combine = F, df2 = "", by.x = "gene_name", by.y = "external_gene_name", all = F){
    df <- getBM(
        filters= "external_gene_name",
        attributes= c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "strand", "transcript_appris"),
        values= ext_name,
        mart = biomart)
    colnames(df) <- c("transcript_id", "gene_id", "gene_name", "strand", "transcript_appris")
    return(df)
}

metadata2wiggleplotr.exons <- function(ensembl_transcript_id, biomart = mart, combine = F, df2 = "", by.x = "gene_name", by.y = "external_gene_name", all = F){
    df <- getBM(
        filters= "ensembl_transcript_id",
        attributes= c("ensembl_transcript_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand"),
        values= ensembl_transcript_id,
        mart = biomart)
    
    df.GRanges <- GRanges(seqnames = df$chromosome_name, ranges = IRanges(start = df$exon_chrom_start, end = df$exon_chrom_end), strand = df$strand)
    mcols(df.GRanges)$transcript_id <- df$ensembl_transcript_id

    return(df.GRanges)
}

metadata2wiggleplotr.utr <- function(ensembl_transcript_id, biomart = mart, combine = F, df2 = "", by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all = F){
df <- getBM(
        filters= "ensembl_transcript_id",
        attributes= c("ensembl_transcript_id", "chromosome_name", "5_utr_start", "5_utr_end", "strand"),
        values= ensembl_transcript_id,
        mart = biomart)

df2 <- getBM(
        filters= "ensembl_transcript_id",
        attributes= c("ensembl_transcript_id", "chromosome_name", "3_utr_start", "3_utr_end", "strand"),
        values= ensembl_transcript_id,
        mart = biomart)

colnames(df) <- c("transcript_id", "chromosome_name", "utr_start", "utr_end", "strand")
colnames(df2) <- colnames(df)

df <- rbind(df, df2)
df <- df[!is.na(df$utr_start) & !is.na(df$utr_end), ]

df.GRanges <- GRanges(seqnames = df$chromosome_name, ranges = IRanges(start = df$utr_start, end = df$utr_end), strand = df$strand)
mcols(df.GRanges)$transcript_id <- df$transcript_id

return(df.GRanges)
}

plotWiggle <- function(gene, log2.transform = F, prior.count = NULL, plot.principal.only = F, heights = c(0.75, 0.25), ...) {
    metadata <- ext_name2wiggleplotr.metadata(gene)
    
    suppressWarnings( if(plot.principal.only == T & !is.na(metadata$transcript_appris)) metadata <- metadata[!metadata$transcript_appris == "", ] )
    
    metadata$transcript_appris <- NULL
    
    # create empty GRangesList for used in creation of exonic and utr regions (for loop)
    all.exons <- GRangesList()
    all.ccds <- GRangesList()
    
    for (i in 1:nrow(metadata)) {
        # get utr regions
        utrs <- metadata2wiggleplotr.utr(metadata$transcript_id[i])
        # get exon regions including utr regions
        exons <- metadata2wiggleplotr.exons(metadata$transcript_id[i])
        
        # convert to GRangesList (required for setdiff)
        utrs <- GRangesList(test = utrs); names(utrs) <- metadata$transcript_id[i]
        exons <- GRangesList(test = exons); names(exons) <- metadata$transcript_id[i]
        
        # find differences beteen exonic and utr regions giving rising to coding regions
        ccds <- setdiff(exons,utrs, ignore.strand = F)
        
        # append transcript exon and utr regions to GRangesList
        all.exons <- append(all.exons, exons, length(all.exons))
        all.ccds <- append(all.ccds, ccds, length(all.ccds))
    }
    
    # Plot coverage 
    if(log2.transform == T) {
        p <- plotCoverage(all.exons, all.ccds, transcript_annotations = metadata, return_subplots_list = T, ...)
        if(is.null(prior.count)) prior.count = 1-min(p$coverage_plot$data$coverage) # add prior count to avoid negative values
        p$coverage_plot$data$coverage <- p$coverage_plot$data$coverage + prior.count
        suppressMessages(p$coverage_plot <- p$coverage_plot + scale_y_continuous(trans=log2_trans())) #scale_y_continuous(trans='log2'))
        suppressMessages(cowplot::plot_grid(p$coverage_plot, p$tx_structure, align = "v", rel_heights = heights, ncol = 1))
    } else {
        plotCoverage(all.exons, all.ccds, transcript_annotations = metadata, ...)    
    }
}
