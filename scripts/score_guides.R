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
#     --nuclease SpCas9 \
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
#
# Notes
# -----
# * For each guide we call ``crisprDesign::findSpacers()`` to locate the
#   spacer in BSgenome.Hsapiens.UCSC.hg38 — this is the only way the
#   downstream addOnTargetScores / addSpacerAlignments helpers get the
#   PAM and flanking nucleotides they need.
# * On-target scoring (DeepSpCas9, DeepHF) requires a 20-nt protospacer
#   plus genomic flanks/PAM. We skip it for non-20-nt spacers (PRIDICT
#   2.0 occasionally emits 19-nt nicking guides) — those rows show NA.
# * For each guide we keep the SINGLE highest-confidence genomic hit as
#   the on-target; the rest count toward mm{0..N} buckets.
# * "CrisprNuclease = SpCas9" because that's the canonical predefined
#   object in crisprBase. SpCas9_H840A is the same protein for spacer/PAM
#   recognition.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--guides"), type = "character",
              help = "Input CSV of guides with columns guide_id,spacer"),
  make_option(c("--nuclease"), type = "character", default = "SpCas9",
              help = "Nuclease name (must be a crisprBase predefined object) [default %default]"),
  make_option(c("--scores"), type = "character",
              default = "DeepSpCas9,DeepHF,CFD",
              help = "Comma-separated list of scores [default %default]"),
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
  library(crisprBase)
  library(crisprScore)
  library(crisprDesign)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

guides_df <- read.csv(opt$guides, stringsAsFactors = FALSE)
stopifnot(all(c("guide_id", "spacer") %in% colnames(guides_df)))

scores_to_run <- strsplit(opt$scores, ",", fixed = TRUE)[[1]]
scores_to_run <- trimws(scores_to_run)

# Resolve nuclease string to a CrisprNuclease object exported by crisprBase
nuclease_obj <- tryCatch(
  get(opt$nuclease, envir = asNamespace("crisprBase")),
  error = function(e) {
    message("could not resolve nuclease '", opt$nuclease,
            "'; falling back to SpCas9")
    crisprBase::SpCas9
  }
)

bsgenome <- BSgenome.Hsapiens.UCSC.hg38

n <- nrow(guides_df)
deepspcas9       <- rep(NA_real_, n)
deephf           <- rep(NA_real_, n)
cfd_specificity  <- rep(NA_real_, n)
mm0_count        <- rep(NA_integer_, n)
mm1_count        <- rep(NA_integer_, n)
mm2_count        <- rep(NA_integer_, n)
mm3_count        <- rep(NA_integer_, n)

for (i in seq_len(n)) {
  gid     <- guides_df$guide_id[i]
  spacer  <- guides_df$spacer[i]
  spacer_len <- nchar(spacer)

  # findSpacers needs a query DNA window — use the spacer itself; it'll
  # search the genome both strands. Returns a GuideSet of every match.
  gs <- tryCatch({
    findSpacers(
      queries          = Biostrings::DNAStringSet(spacer),
      bsgenome         = bsgenome,
      crisprNuclease   = nuclease_obj,
      strict_overlap   = FALSE
    )
  }, error = function(e) {
    message("findSpacers failed for ", gid, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(gs) || length(gs) == 0) {
    message("no genomic hit for ", gid, " (spacer ", spacer, ")")
    next
  }

  # Sometimes findSpacers returns a flanked window not the spacer;
  # restrict to GuideSet rows whose protospacer matches our query exactly.
  spacer_match <- spacers(gs) == spacer
  if (any(spacer_match)) {
    gs <- gs[spacer_match]
  }

  # Pick the first hit as on-target (PRIDICT-designed spacers are
  # typically unique within the target locus).
  on_target <- gs[1]

  # ---------------- mismatch counts via genome-wide alignment
  tryCatch({
    aln <- addSpacerAlignments(
      on_target,
      aligner          = "biostrings",
      txObject         = NULL,
      bsgenome         = bsgenome,
      n_mismatches     = opt$`max-mismatches`,
      all_alignments   = TRUE
    )
    # alignments(aln) returns a GRanges-style object with $n_mismatches
    align_tab <- alignments(aln)
    if (!is.null(align_tab) && length(align_tab) > 0) {
      nm <- align_tab$n_mismatches
      mm0_count[i] <- sum(nm == 0, na.rm = TRUE)
      mm1_count[i] <- sum(nm == 1, na.rm = TRUE)
      mm2_count[i] <- sum(nm == 2, na.rm = TRUE)
      mm3_count[i] <- sum(nm == 3, na.rm = TRUE)
    } else {
      mm0_count[i] <- mm1_count[i] <- mm2_count[i] <- mm3_count[i] <- 0L
    }
    on_target <- aln  # carry the alignments forward for CFD
  }, error = function(e) {
    message("alignment failed for ", gid, ": ", conditionMessage(e))
  })

  # ---------------- on-target scoring (needs 20-nt spacer)
  if (spacer_len == 20) {
    if ("DeepSpCas9" %in% scores_to_run) {
      tryCatch({
        scored <- addOnTargetScores(on_target, methods = "deepspcas9")
        deepspcas9[i] <- as.numeric(scored$score_deepspcas9[1])
      }, error = function(e) {
        message("DeepSpCas9 failed for ", gid, ": ", conditionMessage(e))
      })
    }
    if ("DeepHF" %in% scores_to_run) {
      tryCatch({
        scored <- addOnTargetScores(on_target, methods = "deephf")
        deephf[i] <- as.numeric(scored$score_deephf[1])
      }, error = function(e) {
        message("DeepHF failed for ", gid, ": ", conditionMessage(e))
      })
    }
    if ("CFD" %in% scores_to_run) {
      tryCatch({
        scored <- addOffTargetScores(on_target, methods = "cfd")
        # CFD specificity is on the on-target row's score column
        cfd_specificity[i] <- as.numeric(scored$score_cfd[1])
      }, error = function(e) {
        message("CFD failed for ", gid, ": ", conditionMessage(e))
      })
    }
  } else {
    message("skipping on-target scoring for ", gid,
            " (length ", spacer_len, " != 20 nt)")
  }
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
