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
    why: "Rejects a guide that also matches elsewhere in the genome, which would cut somewhere you did not intend.",
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
    why: "Checks the extension has not stopped the scaffold binding nCas9, by comparing against what a neutral RNA of the same length scores. Longer RNAs always score lower, so the comparison has to be length-matched.",
  },
  {
    id: "report",
    title: "Report",
    engine: "text + CSV",
    gate: false,
    why: "The winning design, broken into labelled parts you can order.",
  },
];

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
export function findBlockingStep(stderr: string, report: string): number {
  const text = `${stderr} ${report}`.toLowerCase();
  return PIPELINE.findIndex((s) => s.keywords?.some((k) => text.includes(k)));
}

/** First per-candidate diagnostic out of a fail-loud NoValidPegRNA message. */
export function blockReason(stderr: string): string {
  const m = stderr.match(/cand\d+:\s*([^;]+)/);
  return m ? m[1].trim() : "";
}
