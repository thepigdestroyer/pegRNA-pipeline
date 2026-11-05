/**
 * The cargo as a list of parts, and which of them may be recoded (ADR 0029).
 *
 * Pure functions only, so the rules can be tested without a browser. The server applies the same rules and
 * refuses what this lets through, so the page and the pipeline cannot disagree.
 */

export interface BuilderPart {
  id: number;
  name: string;
  dna: string;
  /** Never recoded: a promoter, a UTR, a stop codon, or DNA the user wants kept as typed. */
  locked: boolean;
  /** True for a part that has a protein (an amino-acid sequence), so unlocking it can mean something. */
  coding: boolean;
  source: "library" | "custom";
}

let nextId = 1;
export function newPartId(): number {
  return nextId++;
}

export function cargoOf(parts: BuilderPart[]): string {
  return parts.map((p) => p.dna).join("");
}

/** The ``cargo_parts`` a design config takes. */
export function toCargoParts(parts: BuilderPart[]): { name: string; dna: string; locked: boolean }[] {
  return parts.map((p) => ({ name: p.name, dna: p.dna, locked: p.locked }));
}

/**
 * Parts that move the reading frame of coding parts after them: not a multiple of 3 long, with an unlocked part
 * somewhere later. Returns the shift (1 or 2) by part index.
 */
export function frameShifts(parts: BuilderPart[]): Map<number, number> {
  const out = new Map<number, number>();
  parts.forEach((p, i) => {
    const shift = p.dna.length % 3;
    if (shift !== 0 && parts.slice(i + 1).some((q) => !q.locked)) out.set(i, shift);
  });
  return out;
}

export interface ParsedRanges {
  ranges: [number, number][];
  error: string;
}

/**
 * "478-501, 12-20, 7" (1-based, inclusive) -> [[478,501],[12,20],[7,7]]. An empty string is no ranges.
 * Anything that is not a number or a range, or lies outside the cargo, is an error with the reason.
 */
export function parseExcluded(text: string, length: number): ParsedRanges {
  const ranges: [number, number][] = [];
  for (const raw of text.split(",")) {
    const token = raw.trim();
    if (!token) continue;
    const m = /^(\d+)\s*(?:[-–]\s*(\d+))?$/.exec(token);
    if (!m) return { ranges: [], error: `“${token}” is not a nucleotide number or a range like 12-20` };
    const first = Number(m[1]);
    const last = m[2] === undefined ? first : Number(m[2]);
    if (first < 1 || last < first) return { ranges: [], error: `“${token}” is not a valid range` };
    if (length > 0 && last > length)
      return { ranges: [], error: `“${token}” runs past the end of the cargo (${length} nt)` };
    ranges.push([first, last]);
  }
  return { ranges, error: "" };
}

export function rangesToText(ranges: unknown): string {
  if (!Array.isArray(ranges)) return "";
  return ranges
    .filter((r): r is [number, number] => Array.isArray(r) && r.length === 2)
    .map(([a, b]) => (a === b ? String(a) : `${a}-${b}`))
    .join(", ");
}

export function cleanDna(text: string): string {
  return text.replace(/\s/g, "").toUpperCase();
}

export function isDna(text: string): boolean {
  return text.length > 0 && /^[ACGT]+$/.test(text);
}
