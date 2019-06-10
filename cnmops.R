writeLines(">>> Loading Libraries")
library(cn.mops)
writeLines(">>> Loading Libraries [DONE]")

writeLines(">>> Script will be launched with following parameters")
snakemake

writeLines(">>> Analysis of Autosomal CNVs")
writeLines(">>> Getting Read Counts")
bamDataRanges <- getReadCountsFromBAM(BAMFiles=snakemake@input$bam,
                                      sampleNames=snakemake@params$sampleNames,
                                      refSeqNames=snakemake@params$autosomes,
                                      WL=snakemake@params$WL,
                                      parallel=snakemake@threads)
writeLines(">>> Getting Read Counts [DONE]")

writeLines(">>> Calculating Initial CNVs")
resCNMOPS <- cn.mops(bamDataRanges,
                     parallel=snakemake@threads,
                     classes=snakemake@params$classes,
                     segAlgorithm=snakemake@params$segAlgorithm,
                     I=snakemake@params$foldChange,
                     returnPosterior=snakemake@params$returnPosterior,
                    )
warnings()
writeLines(">>> Calculating Initial CNVs [DONE]")

writeLines(">>> Calculating Integer Copy Numbers")
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
writeLines(">>> Calculating Integer Copy Numbers [DONE]")

writeLines(">>> Saving Results")
write.table(cnvs(resCNMOPS), file=snakemake@output$cnvs_autosomes, sep='\t', quote=FALSE, row.names=FALSE)
write.table(cnvr(resCNMOPS), file=snakemake@output$cnvr_autosomes, sep='\t', quote=FALSE, row.names=FALSE)
writeLines(">>> Saving Results [DONE]")

writeLines(">>> Analysis of Heterosomal CNVs")
writeLines(">>> Getting Read Counts")
bamDataRanges <- getReadCountsFromBAM(BAMFiles=snakemake@input$bam,
                                      sampleNames=snakemake@params$sampleNames,
                                      refSeqNames=snakemake@params$heterosomes,
                                      WL=snakemake@params$WL,
                                      parallel=snakemake@threads)
writeLines(">>> Getting Read Counts [DONE]")

writeLines(">>> Normalizing Read Counts")
bamDataRanges <- normalizeChromosomes(bamDataRanges, ploidy=snakemake@params$ploidy)
writeLines(">>> Normalizing Read Counts [DONE]")

writeLines(">>> Calculating Initial CNVs")
resCNMOPS <- cn.mops(bamDataRanges,
                     parallel=snakemake@threads,
                     classes=snakemake@params$classes,
                     segAlgorithm=snakemake@params$segAlgorithm,
                     I=snakemake@params$foldChange,
                     returnPosterior=snakemake@params$returnPosterior,
                     norm=FALSE
                    )
warnings()
writeLines(">>> Calculating Initial CNVs [DONE]")

writeLines(">>> Calculating Integer Copy Numbers")
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
writeLines(">>> Calculating Integer Copy Numbers [DONE]")

writeLines(">>> Saving Results")
write.table(cnvs(resCNMOPS), file=snakemake@output$cnvs_heterosomes, sep='\t', quote=FALSE, row.names=FALSE)
write.table(cnvr(resCNMOPS), file=snakemake@output$cnvr_heterosomes, sep='\t', quote=FALSE, row.names=FALSE)
writeLines(">>> Saving Results [DONE]")

writeLines(">>> Completed Analysis")
