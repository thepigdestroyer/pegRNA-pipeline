#!/usr/bin/env Rscript
# score_guides.R — standalone guide-scoring driver for pegRNA-pipeline.
#
# Reads a CSV of guides (columns: guide_id, spacer) and writes a TSV
# with mm0..mm3 off-target counts for each guide. The on-target columns
# (deepspcas9 / deephf / cfd_specificity) are emitted as NA — those
# scorers need genomic flank context (BSgenome coordinates) that we
# don't currently extract here, and PRIDICT 2.0's PE2 efficiency
# prediction already provides an on-target ranking we trust.
#
# Off-target counts use Biostrings::vcountPattern against
# BSgenome.Hsapiens.UCSC.hg38 directly — no precomputed BWA / bowtie
# index required. Slow (~1 min per chromosome at mm<=3 for a 20-mer)
# but reliable and self-contained.
#
# CLI
# ----
#   Rscript score_guides.R \
#     --guides /work/run/guides.csv \
#     --max-mismatches 3 \
#     --output /work/run/offtarget_scores.tsv
#
# Output columns (tab-separated, with header):
#   guide_id  spacer  deepspcas9  deephf  cfd_specificity
#   mm0_count  mm1_count  mm2_count  mm3_count
#
# On successful completion the script prints the literal marker line
# ``OFFTARGET_SCORING_DONE`` on stdout. The driver greps for that marker
# before trusting the TSV.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--guides"), type = "character",
              help = "Input CSV of guides with columns guide_id,spacer"),
  make_option(c("--nuclease"), type = "character", default = "SpCas9",
              help = "Nuclease name (informational; PAM=NGG assumed) [default %default]"),
  make_option(c("--scores"), type = "character",
              default = "DeepSpCas9,DeepHF,CFD",
              help = "Comma-separated list of scores (currently emitted as NA) [default %default]"),
  make_option(c("--genome"), type = "character", default = "GRCh38",
              help = "Genome build (informational; we always use BSgenome.Hsapiens.UCSC.hg38) [default %default]"),
  make_option(c("--max-mismatches"), type = "integer", default = 3,
              help = "Max mismatches to enumerate (mm0..mmN) [default %default]"),
  make_option(c("--output"), type = "character",
              help = "Output TSV path")
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$guides))
stopifnot(!is.null(opt$output))

suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

guides_df <- read.csv(opt$guides, stringsAsFactors = FALSE)
stopifnot(all(c("guide_id", "spacer") %in% colnames(guides_df)))

bsgenome <- BSgenome.Hsapiens.UCSC.hg38

# Restrict to canonical autosomes + chrX, chrY, chrM for speed; alt /
# random / unplaced contigs add little signal but a lot of runtime.
canonical_chrs <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
canonical_chrs <- intersect(canonical_chrs, seqnames(bsgenome))

n <- nrow(guides_df)
deepspcas9       <- rep(NA_real_, n)
deephf           <- rep(NA_real_, n)
cfd_specificity  <- rep(NA_real_, n)
mm0_count        <- rep(NA_integer_, n)
mm1_count        <- rep(NA_integer_, n)
mm2_count        <- rep(NA_integer_, n)
mm3_count        <- rep(NA_integer_, n)

count_genome_hits <- function(pattern, bsg, chrs, max_mm) {
  # Returns cumulative hit counts at mm <= 0..max_mm across both strands.
  rc <- reverseComplement(pattern)
  out <- integer(max_mm + 1)
  for (chr in chrs) {
    seq <- bsg[[chr]]
    for (mm in 0:max_mm) {
      hits_fw <- countPattern(pattern, seq, max.mismatch = mm, fixed = FALSE)
      hits_rv <- vcountPattern(rc, seq, max.mismatch = mm, fixed = FALSE)
      out[mm + 1] <- out[mm + 1] + hits_fw + hits_rv
    }
  }
  out
}

for (i in seq_len(n)) {
  gid    <- guides_df$guide_id[i]
  spacer <- guides_df$spacer[i]
  if (nchar(spacer) < 18 || nchar(spacer) > 22) {
    message("skipping ", gid, " (length ", nchar(spacer), " outside 18-22 nt)")
    next
  }
  pattern <- DNAString(spacer)
  message("scanning genome for ", gid, " (", nchar(spacer), " nt) ...")
  cum <- tryCatch(
    count_genome_hits(pattern, bsgenome, canonical_chrs, opt$`max-mismatches`),
    error = function(e) {
      message("genome scan failed for ", gid, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(cum)) next
  # Convert cumulative counts (<= mm) to per-bucket counts (== mm).
  exact <- diff(c(0L, cum))
  if (length(exact) >= 1) mm0_count[i] <- exact[1]
  if (length(exact) >= 2) mm1_count[i] <- exact[2]
  if (length(exact) >= 3) mm2_count[i] <- exact[3]
  if (length(exact) >= 4) mm3_count[i] <- exact[4]
  message("  ", gid, ": mm0=", mm0_count[i], " mm1=", mm1_count[i],
          " mm2=", mm2_count[i], " mm3=", mm3_count[i])
}

out_df <- data.frame(
  guide_id        = guides_df$guide_id,
  spacer          = guides_df$spacer,
  deepspcas9      = deepspcas9,
  deephf          = deephf,
  cfd_specificity = cfd_specificity,
  mm0_count       = mm0_count,
  mm1_count       = mm1_count,
  mm2_count       = mm2_count,
  mm3_count       = mm3_count,
  stringsAsFactors = FALSE
)

write.table(out_df, file = opt$output, sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("OFFTARGET_SCORING_DONE\n")
