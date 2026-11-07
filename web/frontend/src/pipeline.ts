import type { ProgressEvent, Strategy } from "./types";

/** What a step of the cascade is doing, for the dot and the word beside it. */
export type StepState = "idle" | "running" | "passed" | "blocked" | "visited" | "unreached";

/**
 * The gate cascade, described for the person reading the screen.
 *
 * `id` matches the gate ids the CLI emits on its `@@P` progress lines and writes
 * into gates.json. `why` is the part the old UI never said out loud: what the
 * gate is actually protecting you from.
 */
export interface PipelineStep {
  id: string;
  title: string;
  engine: string;
  /** False for the plumbing steps that cannot reject a candidate. */
  gate: boolean;
  why: string;
  /** Substrings used to locate the blocking step in a failure message. */
  keywords?: string[];
  /**
   * Prime Assembly: patterns that recognise this step's rejections among the diagnostic lines the CLI records
   * (one line per rejected pair or pegRNA). A step with these is located from the diagnostics, not from prose.
   */
  failure?: RegExp[];
  /** What a rejection at this step is a rejection of, for the summary of why a run failed. */
  unit?: string;
  /**
   * When one diagnostic line reports a number of rejections ("11 failed the off-target screen"), the pattern
   * whose first group is that number. Without it, each matching line is one rejection.
   */
  countFrom?: RegExp;
}

export const PIPELINE: PipelineStep[] = [
  {
    id: "input",
    title: "Design the guide",
    engine: "PRIDICT2",
    gate: false,
    why: "PRIDICT2 picks the pegRNA spacer at your insertion site, and the candidate nicking guides the ngRNA is chosen from.",
  },
  {
    id: "offtarget",
    title: "Off-target check",
    engine: "cas-offinder",
    gate: true,
    keywords: ["off-target", "offtarget"],
    why: "Rejects a second-nick guide that also matches elsewhere in the genome, which would cut somewhere you did not intend. By default no near-match is allowed; the Fine-tuning panel can allow a few, and the result always states how many the chosen guide has.",
  },
  {
    id: "enum",
    title: "Build the candidate pool",
    engine: "codon enumeration",
    gate: false,
    why: "Enumerates nick positions and PBS/RTT variants, ranked so the cheap checks run before the expensive models.",
  },
  {
    id: "recode",
    title: "Recode the arm",
    engine: "synonymous codons",
    gate: false,
    why: "Rewrites the replaced stretch with synonymous codons so the flap no longer matches the genome it has to displace.",
  },
  {
    id: "inv0",
    title: "No deletion (INV-0)",
    engine: "invariant",
    gate: true,
    keywords: ["inv-0", "microhomology", "byte-exact", "locked"],
    why: "Hard guarantee that nothing is deleted from the genome and the recoding stays protein-silent.",
  },
  {
    id: "cryptic",
    title: "Cryptic splice sites",
    engine: "MaxEntScan",
    gate: true,
    keywords: ["cryptic"],
    why: "Catches a new splice site accidentally created by the edit, which would scramble the mRNA.",
  },
  {
    id: "splice",
    title: "Splicing impact",
    engine: "SpliceAI",
    gate: true,
    keywords: ["splice delta"],
    why: "Cross-checks that the edit does not disturb the existing splice sites around it.",
  },
  {
    id: "importance",
    title: "Protein tolerance",
    engine: "AlphaFold-DB",
    gate: true,
    keywords: ["structurally important", "buried"],
    why: "Checks the tag is going somewhere the protein tolerates — exposed and flexible, not buried in the fold.",
  },
  {
    id: "structure",
    title: "pegRNA folding",
    engine: "ViennaRNA",
    gate: true,
    keywords: ["cross-region", "pairing"],
    why: "Rejects a pegRNA that folds back on itself, which stops the template from doing its job.",
  },
  {
    id: "docking",
    title: "Cas9 binding",
    engine: "Boltz-2",
    gate: true,
    keywords: ["interface iptm", "docking"],
    why: "Checks the extension has not stopped the scaffold binding nCas9, by comparing against what a neutral RNA of the same length scores. Longer RNAs always score lower, so the comparison has to be length-matched. Beyond about 230 nt the model can no longer tell a normal RNA from a non-binder, so a longer design is reported against its own length but not accepted or rejected on it.",
  },
  {
    id: "report",
    title: "Report",
    engine: "text + CSV",
    gate: false,
    why: "The winning design, broken into labelled parts you can order.",
  },
];

/**
 * Prime Assembly's cascade, in the order the walk does it. The ids are the gate ids of `design-pa`. `importance` is
 * judged once, up front. Then the guides are paired, screened by the sequence rules and by the off-target scan, and
 * each pair's donor and flaps are built and checked, all before any model is asked about a pegRNA: those steps
 * (`pairs` to `strand`) stream no progress of their own. Last, each pair in turn, best first, has its pegRNAs
 * folded, is screened for splicing, and has them docked (A, then B): these stream per pair.
 */
export const PA_PIPELINE: PipelineStep[] = [
  {
    id: "importance",
    title: "Protein tolerance",
    engine: "AlphaFold-DB",
    gate: true,
    failure: [/structurally important/i, /importance scorer error/i],
    unit: "site",
    why: "Checks the insertion site once, before any guide is screened: the cargo has to go somewhere the protein tolerates, exposed and flexible, not buried in the fold. It does not depend on which guides are used.",
  },
  {
    id: "pairs",
    title: "Pair the guides",
    engine: "geometry",
    gate: true,
    failure: [/^0 PAM-in pairs/i],
    unit: "window",
    why: "Finds every pair of guides on opposite strands that face each other across the insertion point: one cutting at or before it on the top strand, one at or after it on the bottom strand, 30 to 100 bp apart, with their PAMs clear of each other. The pairs are ranked by how few known risks each carries.",
  },
  {
    id: "rules",
    title: "Guide sequence rules",
    engine: "built-in rules",
    gate: true,
    failure: [/failed the sequence rules/i],
    countFrom: /(\d+) failed the sequence rules/,
    unit: "guide",
    why: "Drops a guide whose spacer holds a run of four T (RNA polymerase III reads that as the end of the transcript, so the guide is cut short), or whose GC content is outside 30 to 75%, or whose primer-binding site would run off the end of your sequence window.",
  },
  {
    id: "offtarget",
    title: "Off-target check",
    engine: "cas-offinder",
    gate: true,
    failure: [/failed the off-target screen/i],
    countFrom: /(\d+) failed the off-target screen/,
    unit: "guide",
    why: "Screens every remaining guide against the whole genome in one batch and rejects any that also matches elsewhere. Both guides of a pair must pass. By default no near-match is allowed; the Fine-tuning panel can allow a few, and the result states how many each chosen guide has. A guide with an exact duplicate elsewhere is never allowed.",
  },
  {
    id: "recode",
    title: "Build the donor arms",
    engine: "synonymous codons",
    gate: false,
    failure: [/arm partition failed/i],
    unit: "pair",
    why: "The donor puts back the DNA between the two nicks, so nothing is deleted. Where that stretch codes for protein it is rewritten with synonymous codons, so each flap does not simply match the genome it has to displace; where it does not, it is copied exactly. A guide's PAM inside the stretch is broken where a synonymous change allows it, and only weakened where it does not.",
  },
  {
    id: "inv0",
    title: "No deletion (INV-0)",
    engine: "invariant",
    gate: true,
    failure: [/INV-0/i, /the product is \d+ nt, expected/i],
    unit: "pair",
    why: "Hard guarantee that nothing is deleted from the genome (the donor restores everything it replaces), that the product is exactly your window plus your insert, and that the recoding stays protein-silent.",
  },
  {
    id: "flaps",
    title: "Write the flaps",
    engine: "built-in rules",
    gate: true,
    failure: [/no usable extension/i, /two ends share a run/i],
    unit: "pair",
    why: "Reads each pegRNA's template off the end of the donor, trying flap lengths from 30 to 50 nt, and rejects a template that holds a run of four T (RNA polymerase III would stop there). It also rejects a pair whose donor ends share a long identical stretch, because the two flaps could then pair with each other instead of with the donor.",
  },
  {
    id: "strand",
    title: "Strand check",
    engine: "executable check",
    gate: true,
    failure: [/PA strand mismatch/i],
    unit: "pair",
    why: "Rebuilds the product from the two pegRNAs and the donor, strand by strand, the way the cell would, and checks it is exactly the intended allele. It catches a reverse-complement mistake before it costs an experiment.",
  },
  {
    id: "structure_a",
    title: "pegRNA A folding",
    engine: "ViennaRNA",
    gate: true,
    failure: [/pegRNA A \([^)]*\): (cross-region pairing|structure scorer error)/i],
    unit: "pegRNA",
    why: "Rejects a pegRNA that folds back on itself, which stops the template from doing its job. One pairing is expected and is not counted: because the donor copies the genome at the nick, the end of the template pairs with the end of the spacer (up to 3 bp). The result states how much there was.",
  },
  {
    id: "structure_b",
    title: "pegRNA B folding",
    engine: "ViennaRNA",
    gate: true,
    failure: [/pegRNA B \([^)]*\): (cross-region pairing|structure scorer error)/i],
    unit: "pegRNA",
    why: "The same check for the second guide's pegRNA.",
  },
  {
    id: "cryptic",
    title: "Cryptic splice sites",
    engine: "MaxEntScan",
    gate: true,
    failure: [/cryptic[ -]splice/i],
    unit: "pair",
    why: "Catches a new splice site accidentally created by the edit, which would scramble the mRNA.",
  },
  {
    id: "splice",
    title: "Splicing impact",
    engine: "SpliceAI",
    gate: true,
    failure: [/splice \(SpliceAI|splice cross-check/i],
    unit: "pair",
    why: "Cross-checks that the edit does not disturb the existing splice sites around it. It binds only when it can use the real genomic sequence around your window (give the transcript); without it the number is shown but advisory.",
  },
  {
    id: "docking_a",
    title: "pegRNA A binds Cas9",
    engine: "Boltz-2",
    gate: true,
    failure: [/pegRNA A \([^)]*\): docking/i],
    unit: "pegRNA",
    why: "Checks the extension has not stopped the scaffold binding nCas9, by comparing against what a neutral RNA of the same length scores. A Prime Assembly pegRNA is about 150 nt, well inside the range where the model can tell a normal RNA from a non-binder, so this is a real accept-or-reject, unlike for a long tj-PE design.",
  },
  {
    id: "docking_b",
    title: "pegRNA B binds Cas9",
    engine: "Boltz-2",
    gate: true,
    failure: [/pegRNA B \([^)]*\): docking/i],
    unit: "pegRNA",
    why: "The same check for the second guide's pegRNA. It runs only for a pair whose first pegRNA passed.",
  },
  {
    id: "report",
    title: "Report",
    engine: "text + FASTA + CSV",
    gate: false,
    why: "Both pegRNAs and the donor, broken into labelled parts you can order, with a FASTA of the three sequences, a table of every nucleotide that changed, and notes on how to make the donor.",
  },
];

/** The cascade to draw for a strategy. */
export function pipelineFor(strategy: Strategy): PipelineStep[] {
  return strategy === "pa" ? PA_PIPELINE : PIPELINE;
}

/** The steps a Prime Assembly run repeats for every pair it tries. The others are decided once, while planning. */
export const PA_PER_PAIR = new Set(["structure_a", "structure_b", "cryptic", "splice", "docking_a", "docking_b"]);

/**
 * Planning sends no progress of its own, and nearly all of its time is one step: the batch scan of every candidate
 * guide against the genome, which takes minutes where the steps around it take milliseconds. So once the insertion
 * site is judged, this is the step shown as running.
 */
export const PA_SLOW_PLANNING_STEP = "offtarget";

const paIndex = (gate: string) => PA_PIPELINE.findIndex((s) => s.id === gate);

/**
 * The progress shown after one more event arrives. tj-PE keeps the newest event of each gate. Prime Assembly tries
 * one guide pair after another, and when a new pair starts what the earlier pairs showed no longer describes the
 * run, so their events go; the insertion-site check, which is judged once for the whole run, stays.
 */
export function mergeProgress(
  prev: Record<string, ProgressEvent>,
  ev: ProgressEvent,
  strategy: Strategy,
): Record<string, ProgressEvent> {
  if (strategy !== "pa") return { ...prev, [ev.gate]: ev };
  const next: Record<string, ProgressEvent> = {};
  for (const [gate, old] of Object.entries(prev))
    if (gate === "importance" || old.cand >= ev.cand) next[gate] = old;
  next[ev.gate] = ev;
  return next;
}

/** The newest progress event of a Prime Assembly run: the highest pair, and within it the furthest step. */
export function paLatestEvent(progress: Record<string, ProgressEvent>): ProgressEvent | null {
  return Object.values(progress).reduce<ProgressEvent | null>(
    (best, ev) =>
      best === null || ev.cand > best.cand || (ev.cand === best.cand && paIndex(ev.gate) >= paIndex(best.gate))
        ? ev
        : best,
    null,
  );
}

/**
 * What each Prime Assembly step is doing while a run streams progress, by step id.
 *
 * A step that has reported is passed or blocked. The step shown as running is the NEXT one, because a docking call
 * takes minutes and the step before it is done; after the insertion-site check that is the guide scan, not the
 * instant step that follows it in the list. Planning steps ahead of the running step (or of the last one that
 * reported) are passed without having reported: the walk only gets that far when they are.
 */
export function paLiveStates(progress: Record<string, ProgressEvent>): Record<string, StepState> {
  const latest = paLatestEvent(progress);
  const latestIdx = latest ? paIndex(latest.gate) : -1;
  const running =
    latest === null
      ? 0
      : !latest.passed
        ? -1
        : latest.gate === "importance"
          ? paIndex(PA_SLOW_PLANNING_STEP)
          : latestIdx + 1;
  const frontier = running >= 0 ? running : latestIdx;
  const out: Record<string, StepState> = {};
  PA_PIPELINE.forEach((step, index) => {
    const ev = progress[step.id];
    out[step.id] = ev
      ? ev.passed
        ? "passed"
        : "blocked"
      : index === running
        ? "running"
        : latest && index < frontier && !PA_PER_PAIR.has(step.id)
          ? "passed"
          : "idle";
  });
  return out;
}

export const COMPARISON_SYMBOL: Record<string, string> = {
  ">=": "\u2265",
  "<=": "\u2264",
  "==": "=",
  "<": "<",
  ">": ">",
};

/**
 * Index of the step that blocked a failed run, or -1.
 *
 * The CLI reports a failure as prose, so this matches the step keywords against
 * the combined stderr + report text, exactly as the original page did.
 */
export function findBlockingStep(
  stderr: string,
  report: string,
  steps: PipelineStep[] = PIPELINE,
  diagnostics: string[] = [],
): number {
  // Prime Assembly records one line per rejected pair or pegRNA. Many pairs are tried, so the step to blame is the
  // furthest one any candidate reached: that is where the best candidate stopped.
  if (diagnostics.length > 0 && steps.some((s) => s.failure)) {
    for (let i = steps.length - 1; i >= 0; i--) {
      const patterns = steps[i].failure;
      if (patterns && diagnostics.some((line) => patterns.some((re) => re.test(line)))) return i;
    }
    return -1;
  }
  const text = `${stderr} ${report}`.toLowerCase();
  return steps.findIndex((s) => s.keywords?.some((k) => text.includes(k)));
}

export interface FailureCount {
  /** The step's id, or "other" for rejections that no step claims. */
  id: string;
  title: string;
  unit: string;
  count: number;
}

/** A line the walk writes about one pair ("pair 101/164: ..."), as opposed to a summary of the whole run. */
const ABOUT_ONE_PAIR = /^pair \d+\/\d+:/;

/**
 * How many rejections each Prime Assembly step accounts for, from the diagnostic lines, in pipeline order. A line
 * that reports a number of guides (a step's `countFrom`) counts that many; any other matching line is one rejection.
 * A line about a pair that no step claims is counted as "other" rather than left out, so the counts can fall short
 * of the record only by the summary lines, never by a reason the page failed to recognise.
 */
export function failureCounts(steps: PipelineStep[], diagnostics: string[]): FailureCount[] {
  const out: FailureCount[] = [];
  for (const step of steps) {
    if (!step.failure) continue;
    let count = 0;
    for (const line of diagnostics) {
      if (!step.failure.some((re) => re.test(line))) continue;
      const n = step.countFrom ? line.match(step.countFrom) : null;
      count += n ? Number(n[1]) : 1;
    }
    if (count > 0) out.push({ id: step.id, title: step.title, unit: step.unit ?? "candidate", count });
  }
  const unclaimed = diagnostics.filter(
    (line) => ABOUT_ONE_PAIR.test(line) && !steps.some((s) => s.failure?.some((re) => re.test(line))),
  ).length;
  if (unclaimed > 0) out.push({ id: "other", title: "Other reasons", unit: "pair", count: unclaimed });
  return out;
}

/** First per-candidate diagnostic out of a fail-loud NoValidPegRNA message. */
export function blockReason(stderr: string): string {
  const m = stderr.match(/cand\d+:\s*([^;]+)/);
  return m ? m[1].trim() : "";
}
