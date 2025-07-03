assignChromosomeRegionWithCounts <- function(
    peaks.RD, exon, TSS, utr5, utr3, proximal.promoter.cutoff = c(
      upstream = 2000,
      downstream = 100
    ), immediate.downstream.cutoff = c(
      upstream = 0,
      downstream = 1000
    ), nucleotideLevel = FALSE, precedence = NULL,
    TxDb = NULL, export.file = NULL) {
  if (!is.null(TxDb)) {
    if (!inherits(TxDb, c("TxDb", "EnsDb"))) {
      stop("TxDb must be an object of TxDb or similar such as EnsDb, \n                     try\n?TxDb\tto see more info.")
    }
    if (!inherits(peaks.RD, c("GRanges"))) {
      stop("peaks.RD must be a GRanges object.")
    }
    if (!all(c("upstream", "downstream") %in% names(proximal.promoter.cutoff))) {
      stop("proximal.promoter.cutoff must contain elements upstream and downstream")
    }
    if (!all(c("upstream", "downstream") %in% names(immediate.downstream.cutoff))) {
      stop("immediate.downstream.cutoff must contain elements upstream and downstream")
    }
    if (!is.null(precedence)) {
      if (!all(precedence %in% c(
        "Exons", "Introns", "fiveUTRs",
        "threeUTRs", "Promoters", "immediateDownstream"
      ))) {
        stop("precedence must be a combination of \n                         Exons, Introns, fiveUTRs, threeUTRs, \n                         Promoters, immediateDownstream")
      }
    }
    ignore.strand <- all(as.character(strand(peaks.RD)) == "*")
    exons <- exons(TxDb, columns = NULL)
    introns <- unique(unlist(intronsByTranscript(TxDb)))
    fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
    threeUTRs <- unique(unlist(threeUTRsByTranscript(TxDb)))
    transcripts <- unique(transcripts(TxDb, columns = NULL))
    options(warn = -1)
    try({
      promoters <- unique(promoters(TxDb,
                                    upstream = proximal.promoter.cutoff["upstream"],
                                    downstream = proximal.promoter.cutoff["downstream"]
      ))
      immediateDownstream <- unique(downstreams(transcripts,
                                                upstream = immediate.downstream.cutoff["upstream"],
                                                downstream = immediate.downstream.cutoff["downstream"]
      ))
      promoters <- GenomicRanges::trim(promoters)
      immediateDownstream <- GenomicRanges::trim(immediateDownstream)
    })
    microRNAs <- tryCatch(microRNAs(TxDb), error = function(e) {
      return(NULL)
    })
    tRNAs <- tryCatch(tRNAs(TxDb), error = function(e) {
      return(NULL)
    })
    options(warn = 0)
    annotation <- list(
      exons, introns, fiveUTRs, threeUTRs,
      promoters, immediateDownstream
    )
    if (!is.null(microRNAs)) {
      annotation <- c(annotation, microRNAs = microRNAs)
    }
    if (!is.null(tRNAs)) {
      annotation <- c(annotation, tRNAs = tRNAs)
    }
    annotation <- lapply(annotation, function(.anno) {
      mcols(.anno) <- NULL
      .anno
    })
    names(annotation)[1:6] <- c(
      "Exons", "Introns", "fiveUTRs",
      "threeUTRs", "Promoters", "immediateDownstream"
    )
    # peaks.RD <- formatSeqnames(peaks.RD, exons)
    seqlevelsStyle(peaks.RD) <- seqlevelsStyle(exons)
    peaks.RD <- unique(peaks.RD)
    annotation <- GRangesList(annotation)
    newAnno <- c(unlist(annotation))
    if (ignore.strand) {
      newAnno.rd <- newAnno
      strand(newAnno.rd) <- "*"
      newAnno.rd <- reduce(trim(newAnno.rd))
      Intergenic.Region <- gaps(newAnno.rd, end = seqlengths(TxDb))
      Intergenic.Region <- Intergenic.Region[strand(Intergenic.Region) ==
                                               "*"]
    } else {
      newAnno.rd <- reduce(trim(newAnno))
      Intergenic.Region <- gaps(newAnno.rd, end = seqlengths(TxDb))
      Intergenic.Region <- Intergenic.Region[strand(Intergenic.Region) !=
                                               "*"]
    }
    if (!all(seqlevels(peaks.RD) %in% seqlevels(newAnno))) {
      warning("peaks.RD has sequence levels not in TxDb.")
      sharedlevels <- intersect(seqlevels(newAnno), seqlevels(peaks.RD))
      peaks.RD <- keepSeqlevels(peaks.RD, sharedlevels,
                                pruning.mode = "coarse"
      )
    }
    mcols(peaks.RD) <- NULL
    if (!is.null(precedence)) {
      annotation <- annotation[unique(c(precedence, names(annotation)))]
    }
    names(Intergenic.Region) <- NULL
    annotation$Intergenic.Region <- Intergenic.Region
    anno.names <- names(annotation)
    ol.anno <- findOverlaps(peaks.RD, annotation, ignore.strand = ignore.strand)
    if (nucleotideLevel) {
      jaccardIndex <- unlist(lapply(annotation, function(.ele) {
        intersection <- intersect(.ele, peaks.RD, ignore.strand = ignore.strand)
        union <- union(.ele, peaks.RD, ignore.strand = ignore.strand)
        sum(as.numeric(width(intersection))) / sum(as.numeric(width(union)))
      }))
      jaccardIndex <- jaccardIndex[anno.names]
      names(jaccardIndex) <- anno.names
      jaccardIndex[is.na(jaccardIndex)] <- 0
      newAnno <- unlist(annotation)
      newAnno$source <- rep(names(annotation), lengths(annotation))
      newAnno.disjoin <- disjoin(newAnno,
                                 with.revmap = TRUE,
                                 ignore.strand = ignore.strand
      )
      if (!is.null(precedence)) {
        revmap <- cbind(
          from = unlist(newAnno.disjoin$revmap),
          to = rep(seq_along(newAnno.disjoin), lengths(newAnno.disjoin$revmap))
        )
        revmap <- revmap[order(revmap[, "to"], revmap[
          ,
          "from"
        ]), , drop = FALSE]
        revmap <- revmap[!duplicated(revmap[, "to"]), ,
                         drop = FALSE
        ]
        newAnno.disjoin$source <- newAnno[revmap[, "from"]]$source
      } else {
        revmap <- unlist(newAnno.disjoin$revmap)
        newAnno.disjoin <- rep(newAnno.disjoin, lengths(newAnno.disjoin$revmap))
        newAnno.disjoin$source <- newAnno[revmap]$source
      }
      ol.anno <- findOverlaps(peaks.RD, newAnno.disjoin,
                              ignore.strand = ignore.strand
      )
      queryHits <- peaks.RD[queryHits(ol.anno)]
      subjectHits <- newAnno.disjoin[subjectHits(ol.anno)]
      totalLen <- sum(as.numeric(width(peaks.RD)))
      queryHits.list <- split(queryHits, subjectHits$source)
      lens <- unlist(lapply(queryHits.list, function(.ele) sum(as.numeric(width(unique(.ele))))))
      percentage <- 100 * lens / totalLen
    } else {
      ol.anno.splited <- split(queryHits(ol.anno), anno.names[subjectHits(ol.anno)])
      jaccardIndex <- unlist(lapply(anno.names, function(.name) {
        union <- length(annotation[[.name]]) + length(peaks.RD) -
          length(unique(subjectHits(findOverlaps(peaks.RD,
                                                 annotation[[.name]],
                                                 ignore.strand = ignore.strand
          ))))
        intersection <- length(ol.anno.splited[[.name]])
        intersection / union
      }))
      names(jaccardIndex) <- anno.names
      ol.anno <- as.data.frame(ol.anno)
      ol.anno.splited <- split(ol.anno, ol.anno[, 2])
      hasAnnoHits <- do.call(rbind, ol.anno.splited[names(ol.anno.splited) !=
                                                      as.character(length(annotation))])
      hasAnnoHits <- unique(hasAnnoHits[, 1])
      ol.anno <- ol.anno[!(ol.anno[, 2] == length(annotation) &
                             (ol.anno[, 1] %in% hasAnnoHits)), ]
      if (!is.null(precedence)) {
        ol.anno <- ol.anno[!duplicated(ol.anno[, 1]), ]
      }
      subjectHits <- anno.names[ol.anno[, 2]]
      counts <- table(subjectHits)
      if (!is.null(export.file)) {
        write.csv(as.data.frame(counts), file = export.file, row.names = TRUE)
      }
      percentage <- 100 * counts / length(peaks.RD)
    }
    len <- length(anno.names) - length(percentage)
    if (len > 0) {
      tobeadd <- rep(0, len)
      names(tobeadd) <- anno.names[!anno.names %in% names(percentage)]
      percentage <- c(percentage, tobeadd)
    }
    percentage <- percentage[anno.names]
    return(list(counts = counts, percentage = percentage, jaccard = jaccardIndex))
  }
}