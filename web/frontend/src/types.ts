export type Mode = "fake" | "api";

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
  [key: string]: unknown;
}

export interface Health {
  status: string;
  version: string;
  modes: Mode[];
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
  notes: string[];
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
  returncode: number;
  error?: string;
  stderr?: string;
  stdout?: string;
  report: string;
  components: ComponentRow[];
  has_csv: boolean;
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
