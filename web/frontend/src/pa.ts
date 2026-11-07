/**
 * Prime Assembly's settings, rules and wording, as pure functions.
 *
 * The validation below states the rules the server applies (`validate_pa_request`), so the page can say what is
 * wrong before a job is started; the server refuses whatever this lets through, so the two cannot disagree for
 * long. The numbers are the engine's: the flap range is the one the paper found robust, the span floor is the
 * project's, and the caps on pairs and PBS lengths bound what one request may spend on paid docking calls.
 */
import type { DesignConfig, PaGuide, Strategy } from "./types";

export const PA_DEFAULTS = {
  flap_length: "35",
  min_span: "30",
  max_span: "100",
  max_offtargets: "0",
  max_pairs: "30",
} as const;

export const PA_LIMITS = {
  flapMin: 30,
  flapMax: 50,
  spanFloor: 23,
  spanCap: 200,
  maxOfftargets: 10,
  maxPairs: 60,
  pbsCount: 5,
  pbsMin: 6,
  pbsMax: 30,
} as const;

/** The fields Prime Assembly adds to the form, as the form holds them (text). */
export interface PaFields {
  flap_length: string;
  min_span: string;
  max_span: string;
  max_offtargets: string;
  max_pairs: string;
  pbs_lengths: string;
  max_run: string;
}

export type PaFieldName = keyof PaFields | "cargo";

/** Keys that mean something to tj-PE only. A Prime Assembly request does not carry them. */
export const TJPE_ONLY_KEYS = [
  "rtt_aa",
  "pegrna_strand",
  "arm_side",
  "ngrna_distance_min",
  "ngrna_distance_max",
  "ngrna_max_offtargets",
  "spacer",
  "het_snps",
  "het_edit",
] as const;

/** Keys that mean something to Prime Assembly only. A tj-PE request does not carry them. */
export const PA_ONLY_KEYS = [
  "flap_length",
  "flap_min",
  "flap_max",
  "min_span",
  "max_span",
  "max_offtargets",
  "max_pairs",
] as const;

function whole(text: string): number | null {
  const t = text.trim();
  if (!/^-?\d+$/.test(t)) return null;
  return Number(t);
}

/** The PBS lengths typed as "13,15,17": every entry must be a whole number, or the list is not valid. */
export function parsePbs(text: string): number[] | null {
  const parts = text
    .split(",")
    .map((s) => s.trim())
    .filter((s) => s.length > 0);
  if (parts.length === 0) return null;
  const out: number[] = [];
  for (const p of parts) {
    const n = whole(p);
    if (n === null) return null;
    out.push(n);
  }
  return out;
}

/** Why the Prime Assembly settings cannot be used, by field. Empty when they can. */
export function validatePa(f: PaFields, cargoLength: number): Partial<Record<PaFieldName, string>> {
  const e: Partial<Record<PaFieldName, string>> = {};
  const flap = whole(f.flap_length);
  if (flap === null || flap < PA_LIMITS.flapMin || flap > PA_LIMITS.flapMax)
    e.flap_length = `A whole number from ${PA_LIMITS.flapMin} to ${PA_LIMITS.flapMax}.`;
  else if (cargoLength > 0 && cargoLength < 2 * flap)
    e.cargo = `Prime Assembly needs an insert of at least ${2 * flap} nt (two ${flap} nt flaps that must not overlap). This one is ${cargoLength} nt: write a short insert inside the pegRNA instead.`;

  const lo = whole(f.min_span);
  const hi = whole(f.max_span);
  if (lo === null || lo < PA_LIMITS.spanFloor)
    e.min_span = `A whole number, at least ${PA_LIMITS.spanFloor} bp.`;
  if (hi === null || hi > PA_LIMITS.spanCap) e.max_span = `A whole number, at most ${PA_LIMITS.spanCap} bp.`;
  else if (lo !== null && hi < lo) e.max_span = "Must be at least the minimum.";

  const off = whole(f.max_offtargets);
  if (off === null || off < 0 || off > PA_LIMITS.maxOfftargets)
    e.max_offtargets = `A whole number from 0 to ${PA_LIMITS.maxOfftargets}.`;

  const pairs = whole(f.max_pairs);
  if (pairs === null || pairs < 1 || pairs > PA_LIMITS.maxPairs)
    e.max_pairs = `A whole number from 1 to ${PA_LIMITS.maxPairs}.`;

  const pbs = parsePbs(f.pbs_lengths);
  if (
    pbs === null ||
    pbs.length > PA_LIMITS.pbsCount ||
    pbs.some((n) => n < PA_LIMITS.pbsMin || n > PA_LIMITS.pbsMax)
  )
    e.pbs_lengths = `Up to ${PA_LIMITS.pbsCount} whole numbers from ${PA_LIMITS.pbsMin} to ${PA_LIMITS.pbsMax}, separated by commas.`;

  const run = whole(f.max_run);
  if (run === null || run < 1) e.max_run = "A whole number, at least 1.";
  return e;
}

/**
 * The config a Prime Assembly request carries: the shared fields as given, the Prime Assembly settings named as the
 * engine's `design_pa` takes them, and nothing that only tj-PE understands. The caller has validated the fields.
 */
export function toPaConfig(cfg: DesignConfig, f: PaFields): DesignConfig {
  const out: DesignConfig = { ...cfg };
  for (const k of TJPE_ONLY_KEYS) delete out[k];
  out.flap_length = Number(f.flap_length);
  out.min_span = Number(f.min_span);
  out.max_span = Number(f.max_span);
  out.max_offtargets = Number(f.max_offtargets);
  out.max_pairs = Number(f.max_pairs);
  out.max_run = Number(f.max_run);
  const pbs = parsePbs(f.pbs_lengths);
  if (pbs) out.pbs_lengths = pbs;
  return out;
}

/** The config a tj-PE request carries: nothing that only Prime Assembly understands. */
export function toTjpeConfig(cfg: DesignConfig): DesignConfig {
  const out: DesignConfig = { ...cfg };
  for (const k of PA_ONLY_KEYS) delete out[k];
  return out;
}

export interface Advice {
  tone: "info" | "warn";
  text: string;
  /** The strategy the advice is for, when following it means switching. */
  switchTo?: Strategy;
  action?: string;
}

/**
 * Which strategy fits an insert of this length. It only advises: the choice is the user's. The numbers: the two
 * flaps of a Prime Assembly design must not overlap, so the insert must be at least twice the flap; a pegRNA that
 * carries the insert itself is about the insert plus 150 nt, and past about 230 nt the Cas9-binding check can no
 * longer accept or reject it.
 */
export function strategyAdvice(strategy: Strategy, cargoLength: number, flapLength: number): Advice | null {
  if (cargoLength <= 0) return null;
  const minimum = 2 * flapLength;
  if (strategy === "pa") {
    if (cargoLength >= minimum) return null;
    return {
      tone: "warn",
      text: `This insert is ${cargoLength} nt. Prime Assembly needs at least ${minimum} nt (two ${flapLength} nt flaps that must not overlap). For an insert this short, write it inside the pegRNA.`,
      switchTo: "tjpe",
      action: "Write it inside the pegRNA",
    };
  }
  if (cargoLength >= 100) {
    return {
      tone: "info",
      text: `A ${cargoLength} nt insert makes a very long pegRNA, roughly ${cargoLength + 150} nt: past about 230 nt the Cas9-binding check can no longer accept or reject it. Prime Assembly writes only two short flaps and supplies the rest as a donor DNA, so each pegRNA stays about 150 nt.`,
      switchTo: "pa",
      action: "Use Prime Assembly",
    };
  }
  if (cargoLength >= minimum) {
    return {
      tone: "info",
      text: `At ${cargoLength} nt either works. Inside the pegRNA is one RNA to order (roughly ${cargoLength + 150} nt); Prime Assembly is two short pegRNAs plus a donor DNA to make.`,
      switchTo: "pa",
      action: "Use Prime Assembly",
    };
  }
  return null;
}

/** The one warning that is not a finding about this design: it is why the whole combination is untested. */
export function isStandingCaveat(warning: string): boolean {
  return warning.startsWith("A donor that restores the DNA between the nicks");
}

export interface Fate {
  /** "ok" the guide cannot cut the product; "warn" it still can, fully or at a fraction of its rate. */
  tone: "ok" | "warn";
  text: string;
}

/** What happens to a guide's PAM in the product, in words: whether the guide can cut the product again. */
export function pamFate(g: PaGuide): Fate {
  switch (g.pam_status) {
    case "broken":
      return { tone: "ok", text: "Broken in the donor by one synonymous change: this guide cannot cut the product." };
    case "broken-by-recode":
      return { tone: "ok", text: "Broken by the recoding of the arm: this guide cannot cut the product." };
    case "weakened":
      return {
        tone: "warn",
        text: `Weakened to ${g.product_pam}: SpCas9 still cuts that at a fraction of its normal rate (about a fifth for NAG), so this guide may still nick the product.`,
      };
    case "intact-noncoding":
      return {
        tone: "warn",
        text: "Regenerated: it lies in non-coding DNA, which is kept exactly, so this guide can cut the product again.",
      };
    case "intact-no-synonymous-change":
      return {
        tone: "warn",
        text: "Regenerated: it lies in coding sequence with no synonymous change that breaks it, so this guide can cut the product again.",
      };
  }
}
