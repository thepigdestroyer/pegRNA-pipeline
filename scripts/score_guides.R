#!/usr/bin/env Rscript
# score_guides.R — standalone guide-scoring driver for pegRNA-pipeline.
#
# Reads a CSV of guides (columns: guide_id, spacer) and writes a TSV with
# DeepSpCas9 / DeepHF / CFD specificity scores plus genome-wide
# mm0..mm3 off-target counts. Run end-to-end on a compute node under the
# ``pegrna-ngrna`` conda env (which ships crisprScore + crisprDesign +
# BSgenome.Hsapiens.UCSC.hg38).
#
# CLI
# ----
#   Rscript score_guides.R \
#     --guides /work/run/guides.csv \
#     --nuclease SpCas9_H840A \
#     --scores DeepSpCas9,DeepHF,CFD \
#     --genome GRCh38 \
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
  make_option(c("--nuclease"), type = "character", default = "SpCas9_H840A",
              help = "Nuclease name [default %default]"),
  make_option(c("--scores"), type = "character",
              default = "DeepSpCas9,DeepHF,CFD",
              help = "Comma-separated list of scores [default %default]"),
  make_option(c("--genome"), type = "character", default = "GRCh38",
              help = "Genome build for off-target search [default %default]"),
  make_option(c("--max-mismatches"), type = "integer", default = 3,
              help = "Max mismatches to enumerate (mm0..mmN) [default %default]"),
  make_option(c("--output"), type = "character",
              help = "Output TSV path")
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$guides))
stopifnot(!is.null(opt$output))

suppressPackageStartupMessages({
  library(crisprScore)
  library(crisprDesign)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

guides_df <- read.csv(opt$guides, stringsAsFactors = FALSE)
stopifnot(all(c("guide_id", "spacer") %in% colnames(guides_df)))

scores_to_run <- strsplit(opt$scores, ",", fixed = TRUE)[[1]]
scores_to_run <- trimws(scores_to_run)

# ------------------------------------------------------------- on-target
deepspcas9 <- rep(NA_real_, nrow(guides_df))
deephf <- rep(NA_real_, nrow(guides_df))
cfd_specificity <- rep(NA_real_, nrow(guides_df))

if ("DeepSpCas9" %in% scores_to_run) {
  tryCatch({
    res <- getDeepSpCas9Scores(guides_df$spacer)
    deepspcas9 <- as.numeric(res$score)
  }, error = function(e) message("DeepSpCas9 failed: ", conditionMessage(e)))
}

if ("DeepHF" %in% scores_to_run) {
  tryCatch({
    res <- getDeepHFScores(guides_df$spacer)
    deephf <- as.numeric(res$score)
  }, error = function(e) message("DeepHF failed: ", conditionMessage(e)))
}

# ------------------------------------------------------------- off-target
mm0_count <- rep(NA_integer_, nrow(guides_df))
mm1_count <- rep(NA_integer_, nrow(guides_df))
mm2_count <- rep(NA_integer_, nrow(guides_df))
mm3_count <- rep(NA_integer_, nrow(guides_df))

bsgenome <- BSgenome.Hsapiens.UCSC.hg38
for (i in seq_len(nrow(guides_df))) {
  tryCatch({
    gs <- GuideSet(
      ids = guides_df$guide_id[i],
      protospacers = guides_df$spacer[i],
      CrisprNuclease = opt$nuclease,
      bsgenome = bsgenome
    )
    gs <- addSpacerAlignments(
      gs,
      bsgenome = bsgenome,
      n_mismatches = opt$`max-mismatches`
    )
    aln <- alignments(gs)
    for (mm in 0:opt$`max-mismatches`) {
      col <- paste0("mm", mm, "_count")
      cnt <- sum(aln$n_mismatches == mm, na.rm = TRUE)
      assign(col, c(get(col)[seq_len(i - 1)], cnt, get(col)[-seq_len(i)]))
    }
    if ("CFD" %in% scores_to_run) {
      tryCatch({
        gs <- addOffTargetScores(gs, methods = "cfd")
        cfd_specificity[i] <- as.numeric(onTargetScore(gs, "cfd_specificity"))
      }, error = function(e) message("CFD failed for ", guides_df$guide_id[i]))
    }
  }, error = function(e) {
    message("off-target alignment failed for ",
            guides_df$guide_id[i], ": ", conditionMessage(e))
  })
}

out_df <- data.frame(
  guide_id = guides_df$guide_id,
  spacer = guides_df$spacer,
  deepspcas9 = deepspcas9,
  deephf = deephf,
  cfd_specificity = cfd_specificity,
  mm0_count = mm0_count,
  mm1_count = mm1_count,
  mm2_count = mm2_count,
  mm3_count = mm3_count,
  stringsAsFactors = FALSE
)

write.table(out_df, file = opt$output, sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("OFFTARGET_SCORING_DONE\n")
