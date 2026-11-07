export type Mode = "fake" | "api";

/**
 * How the insert reaches the genome. "tjpe": it rides in the pegRNA's own extension (tags, short inserts).
 * "pa": Prime Assembly, two pegRNAs write short flaps and a separate donor DNA, made by PCR, supplies the rest.
 */
export type Strategy = "tjpe" | "pa";

export interface DesignConfig {
  gene?: string;
  transcript?: string;
  uniprot_acc?: string;
  rtt_aa?: string;
  cargo?: string;
  genome?: string;
  site_s?: number;
  cds_mask?: boolean[];
  pbs_lengths?: number[];
  ngrna_distance_min?: number;
  ngrna_distance_max?: number;
  pegrna_strand?: "Fw" | "Rv";
  arm_side?: "5prime" | "3prime";
  max_run?: number;
  spacer?: string;
  het_snps?: unknown[];
  het_edit?: boolean;
  protein_length?: number;
  structure_threshold?: number;
  /** Off-target sites the second nick (ngRNA) may have; 0 is the strict default. */
  ngrna_max_offtargets?: number;
  // Prime Assembly only
  /** Nucleotides each pegRNA writes, and so how much of each end of the donor must match: 30 to 50. */
  flap_length?: number;
  /** Nick-to-nick distance limits, in bp. */
  min_span?: number;
  max_span?: number;
  /** Off-target sites either guide may have; 0 is the strict default. */
  max_offtargets?: number;
  /** How many guide pairs the walk may try, best first. */
  max_pairs?: number;
  [key: string]: unknown;
}

export interface Health {
  status: string;
  version: string;
  modes: Mode[];
  /** Absent from a server that predates Prime Assembly. */
  strategies?: Strategy[];
  default_mode: Mode;
  repo_ok: boolean;
  api_python_ok: boolean;
  nvidia_key: boolean;
  concurrency: { limit: number };
}

export interface ComponentRow {
  index?: string;
  component?: string;
  role?: string;
  part_of?: string;
  sequence?: string;
  length?: string;
  [key: string]: string | undefined;
}

export interface GateMetric {
  metric?: string;
  value?: number | string | null;
  cutoff?: number | string | null;
  cmp?: string;
  extra?: string;
  note?: string;
}

/** One `@@P` progress line from the CLI, emitted as each gate is evaluated. */
export interface ProgressEvent {
  cand: number;
  gate: string;
  info: string;
  passed: boolean;
  metrics?: GateMetric;
}

/** One nucleotide the pipeline changed from what was specified (ADR 0029). Positions are 1-based. */
export interface ChangeRow {
  scope: "cargo" | "genomic arm";
  part: string;
  position: number;
  codon: number | null;
  old_codon: string;
  new_codon: string;
  amino_acid: string;
  old_base: string;
  new_base: string;
}

export interface ChangesPayload {
  cargo_unchanged: boolean;
  cargo_variant: number;
  cargo_length: number;
  cargo_changes: ChangeRow[];
  arm_changes: ChangeRow[];
  /** What was protected, one line each: kept exactly as specified. */
  kept_exactly: string[];
  /** Stretches of the genomic arm(s) that were copied exactly (Prime Assembly's donor arms). */
  arm_kept_exactly?: string[];
  notes: string[];
}

export type PamStatus =
  | "broken"
  | "broken-by-recode"
  | "weakened"
  | "intact-noncoding"
  | "intact-no-synonymous-change";

/** One of the two guides of a Prime Assembly design. Positions are 0-based window coordinates. */
export interface PaGuide {
  side: "A" | "B";
  strand: "top" | "bottom";
  spacer: string;
  /** The PAM in the genome, and the PAM the product ends up with (as the guide reads it). */
  pam: string;
  product_pam: string;
  pam_status: PamStatus;
  /** The first base to the right of the nick. */
  cut: number;
  flap_len: number;
  pbs_len: number;
  pbs: string;
  rtt: string;
  tail: string;
  pegRNA: string;
  pegRNA_length: number;
  offtargets: number;
  on_target: number;
  genomic_run: number;
  /** Base pairs of the RTT's end with the spacer's end that the fold predicted and did not count. */
  nick_pairing_bp: number;
}

export interface PaDonorPart {
  name: "arm5" | "cargo" | "arm3";
  start: number;
  end: number;
  sequence: string;
}

export interface PaPayload {
  site_s: number;
  /** The two nicks: the donor replaces window[span[0]:span[1]]. */
  span: [number, number];
  window_length: number;
  cargo_length: number;
  donor_length: number;
  product_length: number;
  /** How to make the donor, from the published method. */
  donor_prep: string[];
  guides: [PaGuide, PaGuide];
  donor: {
    sequence: string;
    parts: PaDonorPart[];
    /** The end of the donor that each flap matches, as [start, end) in donor coordinates. */
    flap_a: [number, number];
    flap_b: [number, number];
  };
  recode_spans: [number, number][];
  locked_spans: [number, number][];
  simulated: boolean;
}

export interface GatesPayload {
  gate_results?: [string, string][];
  /** Everything that differs from the input, exactly; null for a run that did not produce a design. */
  changes?: ChangesPayload | null;
  recode_edits?: number;
  deletion_size?: number;
  candidates_tried?: number;
  warnings?: string[];
  /** Gates served by a stand-in rather than a real model. */
  stand_in_gates?: string[];
  scorer_mode?: string;
  strategy?: Strategy;
  /** Set by a run that did not produce a design: the CLI's own record of why. */
  failed?: boolean;
  error?: string;
  diagnostics?: string[];
  pairs_considered?: number;
  /** Prime Assembly: both guides, the donor and its parts. */
  pa?: PaPayload;
}

export interface HetEdit {
  requested: boolean;
  produced: boolean;
  note: string;
  report: string;
  components: ComponentRow[];
}

export interface DesignResult {
  type?: "result";
  ok: boolean;
  run_id: string;
  mode: Mode;
  /** Absent from a server that predates Prime Assembly (which means "tjpe"). */
  strategy?: Strategy;
  returncode: number;
  error?: string;
  stderr?: string;
  stdout?: string;
  report: string;
  components: ComponentRow[];
  has_csv: boolean;
  has_fasta?: boolean;
  docking_stub: boolean;
  gates: GatesPayload;
  het_edit: HetEdit | null;
}

export interface Motif {
  name: string;
  aa: string;
  dna: string;
  description: string;
  /** Not protein (a promoter, a UTR, a stop codon): never recoded. */
  locked?: boolean;
}

export type MotifGroups = Record<string, Motif[]>;
