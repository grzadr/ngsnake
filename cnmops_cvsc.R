writeLines(">>> Loading Libraries")
library(cn.mops)
writeLines(">>> Loading Libraries [DONE]")

writeLines(">>> Script will be launched with following parameters")
snakemake

writeLines(">>> Analysis of Autosomal CNVs")
writeLines(">>> Getting Case Read Counts")
bamDataRangesCase <- getReadCountsFromBAM(BAMFiles=snakemake@input$caseBAM,
                                          sampleNames=snakemake@params$caseSampleNames,
                                          refSeqNames=snakemake@params$autosomes,
                                          WL=snakemake@params$WL,
                                          parallel=snakemake@threads)
warnings()
writeLines(">>> Getting Case Read Counts [DONE]")

writeLines(">>> Getting Control Read Counts")
bamDataRangesControl <- getReadCountsFromBAM(BAMFiles=snakemake@input$controlBAM,
                                             sampleNames=snakemake@params$controlSampleNames,
                                             refSeqNames=snakemake@params$autosomes,
                                             WL=snakemake@params$WL,
                                             parallel=snakemake@threads)
warnings()
writeLines(">>> Getting Control Read Counts [DONE]")

writeLines(">>> Calculating Initial CNVs")
resCNMOPS <- referencecn.mops(cases=bamDataRangesCase,
                              controls=bamDataRangesControl,
                              parallel=snakemake@threads,
                              classes=snakemake@params$classes,
                              segAlgorithm=snakemake@params$segAlgorithm,
                              I=snakemake@params$foldChange,
                              returnPosterior=snakemake@params$returnPosterior)
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
writeLines(">>> Getting Case Read Counts")

bamDataRangesCase <- getReadCountsFromBAM(BAMFiles=snakemake@input$caseBAM,
                                          sampleNames=snakemake@params$caseSampleNames,
                                          refSeqNames=snakemake@params$heterosomes,
                                          WL=snakemake@params$WL,
                                          parallel=snakemake@threads)
warnings()
writeLines(">>> Getting Case Read Counts [DONE]")

writeLines(">>> Getting Control Read Counts")
bamDataRangesControl <- getReadCountsFromBAM(BAMFiles=snakemake@input$controlBAM,
                                             sampleNames=snakemake@params$controlSampleNames,
                                             refSeqNames=snakemake@params$heterosomes,
                                             WL=snakemake@params$WL,
                                             parallel=snakemake@threads)
warnings()
writeLines(">>> Getting Control Read Counts [DONE]")

writeLines(">>> Normalizing Case Read Counts")
bamDataRangesCase <- normalizeChromosomes(bamDataRangesCase,
                                          ploidy=snakemake@params$casePloidy)
warnings()
writeLines(">>> Normalizing Case Read Counts [DONE]")

writeLines(">>> Normalizing Control Read Counts")
bamDataRangesControl <- normalizeChromosomes(bamDataRangesControl,
                                             ploidy=snakemake@params$controlPloidy)
warnings()
writeLines(">>> Normalizing Control Read Counts [DONE]")

writeLines(">>> Calculating Initial CNVs")
resCNMOPS <- referencecn.mops(cases=bamDataRangesCase,
                              controls=bamDataRangesControl,
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
