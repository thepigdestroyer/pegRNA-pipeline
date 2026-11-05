import { useMemo, useState } from "react";
import type { Motif, MotifGroups } from "../types";
import {
  type BuilderPart,
  cargoOf,
  cleanDna,
  frameShifts,
  isDna,
  newPartId,
  parseExcluded,
} from "../cargo";
import { Badge, CopyButton, Field, cx } from "./Primitives";

/**
 * Build the cargo from the motif library and your own sections, or paste DNA.
 *
 * Every part is either protein-coding (the recoder may swap synonymous codons) or locked (kept exactly as it
 * is, always). A promoter, a UTR or a stop codon is locked by default and a library protein is recodable by
 * default; either can be changed per part. Your own sections are locked unless you say they are protein.
 * Pasted DNA can have stretches excluded by nucleotide number.
 *
 * Parts are assembled 5'->3' and the running frame check catches the mistake that actually bites: a part whose
 * length is not a multiple of three, with coding sequence after it, silently shifts everything downstream.
 */

interface Props {
  motifs: MotifGroups;
  parts: BuilderPart[];
  onPartsChange: (parts: BuilderPart[]) => void;
  rawCargo: string;
  onRawCargoChange: (dna: string) => void;
  rawExcluded: string;
  onRawExcludedChange: (text: string) => void;
  mode: "builder" | "raw";
  onModeChange: (mode: "builder" | "raw") => void;
}

function libraryPart(m: Motif): BuilderPart {
  const protein = Boolean(m.aa) && m.aa !== "*";
  return {
    id: newPartId(),
    name: m.name,
    dna: m.dna,
    locked: m.locked ?? !protein,
    coding: protein,
    source: "library",
  };
}

export function CargoBuilder({
  motifs,
  parts,
  onPartsChange,
  rawCargo,
  onRawCargoChange,
  rawExcluded,
  onRawExcludedChange,
  mode,
  onModeChange,
}: Props) {
  const [pick, setPick] = useState("");
  const [custom, setCustom] = useState({ name: "", dna: "", coding: false });
  const [customError, setCustomError] = useState("");

  const byName = useMemo(() => {
    const m = new Map<string, Motif>();
    Object.values(motifs).forEach((list) => list.forEach((x) => m.set(x.name, x)));
    return m;
  }, [motifs]);

  const cargo = mode === "raw" ? rawCargo : cargoOf(parts);
  const frameOffset = cargo.length % 3;
  const shifts = useMemo(() => frameShifts(parts), [parts]);
  const excluded = useMemo(() => parseExcluded(rawExcluded, rawCargo.length), [rawExcluded, rawCargo]);
  const lockedCount = parts.filter((p) => p.locked).length;

  function add() {
    const m = byName.get(pick);
    if (m) onPartsChange([...parts, libraryPart(m)]);
  }
  function addCustom() {
    const dna = cleanDna(custom.dna);
    if (!isDna(dna)) return setCustomError("Enter the DNA as A, C, G and T only.");
    if (custom.coding && dna.length % 3 !== 0)
      return setCustomError(
        "A protein-coding section must be a whole number of codons (a multiple of 3). Untick it to keep it exactly as typed.",
      );
    setCustomError("");
    onPartsChange([
      ...parts,
      {
        id: newPartId(),
        name: custom.name.trim() || `my section ${parts.filter((p) => p.source === "custom").length + 1}`,
        dna,
        locked: !custom.coding,
        coding: true,
        source: "custom",
      },
    ]);
    setCustom({ name: "", dna: "", coding: false });
  }
  function remove(id: number) {
    onPartsChange(parts.filter((p) => p.id !== id));
  }
  function toggleLock(id: number) {
    onPartsChange(parts.map((p) => (p.id === id ? { ...p, locked: !p.locked } : p)));
  }
  function move(i: number, delta: number) {
    const j = i + delta;
    if (j < 0 || j >= parts.length) return;
    const next = [...parts];
    [next[i], next[j]] = [next[j], next[i]];
    onPartsChange(next);
  }

  return (
    <div className="space-y-4">
      <div className="flex flex-wrap items-center gap-2">
        <div
          className="inline-flex rounded-lg border border-ink-300 p-0.5 dark:border-ink-700"
          role="tablist"
          aria-label="How to specify the cargo"
        >
          {(["builder", "raw"] as const).map((m) => (
            <button
              key={m}
              type="button"
              role="tab"
              aria-selected={mode === m}
              onClick={() => onModeChange(m)}
              className={cx(
                "rounded-md px-3 py-1 text-xs font-medium transition-colors",
                mode === m
                  ? "bg-brand-600 text-white"
                  : "text-ink-600 hover:text-ink-900 dark:text-ink-400 dark:hover:text-ink-100",
              )}
            >
              {m === "builder" ? "Build from parts" : "Paste DNA"}
            </button>
          ))}
        </div>
        <span className="text-xs text-ink-500 dark:text-ink-400">
          {mode === "builder"
            ? "Parts are joined 5′→3′ in the order listed."
            : "Enter the cargo exactly as it should be inserted."}
        </span>
      </div>

      {mode === "builder" && (
        <>
          <div className="flex flex-wrap gap-2">
            <label className="sr-only" htmlFor="motif-pick">
              Choose a part to add
            </label>
            <select
              id="motif-pick"
              className="input flex-1"
              value={pick}
              onChange={(e) => setPick(e.target.value)}
            >
              <option value="">Choose a part from the library…</option>
              {Object.keys(motifs)
                .sort()
                .map((cat) => (
                  <optgroup key={cat} label={cat.replace(/_/g, " ")}>
                    {motifs[cat].map((m) => (
                      <option key={m.name} value={m.name}>
                        {m.name}
                        {m.aa ? ` (${m.aa.length > 16 ? `${m.aa.slice(0, 16)}…` : m.aa})` : ""}
                        {m.locked ? " · locked" : ""}
                      </option>
                    ))}
                  </optgroup>
                ))}
            </select>
            <button type="button" className="btn btn-secondary" onClick={add} disabled={!pick}>
              Add part
            </button>
          </div>

          <details className="rounded-lg border border-ink-200 px-3 py-2 dark:border-ink-800">
            <summary className="cursor-pointer text-sm font-medium text-ink-700 dark:text-ink-200">
              Add your own section
            </summary>
            <div className="mt-3 space-y-3">
              <div className="grid gap-3 sm:grid-cols-[12rem_1fr]">
                <Field label="Name" htmlFor="custom-name">
                  <input
                    id="custom-name"
                    className="input text-sm"
                    value={custom.name}
                    onChange={(e) => setCustom({ ...custom, name: e.target.value })}
                    placeholder="my promoter"
                  />
                </Field>
                <Field label="DNA" htmlFor="custom-dna">
                  <input
                    id="custom-dna"
                    className="input seq text-xs"
                    value={custom.dna}
                    onChange={(e) => setCustom({ ...custom, dna: cleanDna(e.target.value) })}
                    placeholder="ATGGCT…"
                  />
                </Field>
              </div>
              <label className="flex items-start gap-2 text-sm text-ink-700 dark:text-ink-200">
                <input
                  type="checkbox"
                  className="mt-1"
                  checked={custom.coding}
                  onChange={(e) => setCustom({ ...custom, coding: e.target.checked })}
                />
                <span>
                  This is protein-coding. The recoder may swap synonymous codons in it.
                  <span className="block text-xs text-ink-500 dark:text-ink-400">
                    Leave unticked to keep it exactly as typed. That is the safe choice for anything that is
                    not protein.
                  </span>
                </span>
              </label>
              {customError && (
                <p role="alert" className="text-xs font-medium text-rose-600">
                  {customError}
                </p>
              )}
              <button type="button" className="btn btn-secondary" onClick={addCustom}>
                Add section
              </button>
            </div>
          </details>

          {parts.length === 0 ? (
            <p className="rounded-lg border border-dashed border-ink-300 px-4 py-6 text-center text-sm text-ink-500 dark:border-ink-700 dark:text-ink-400">
              No parts yet. Add a tag, linker or marker above, or your own section.
            </p>
          ) : (
            <ul className="space-y-1.5">
              {parts.map((p, i) => {
                const m = byName.get(p.name);
                const shift = shifts.get(i);
                return (
                  <li
                    key={p.id}
                    className="flex flex-wrap items-center gap-x-3 gap-y-1 rounded-lg border border-ink-200 bg-ink-50
                               px-3 py-2 dark:border-ink-800 dark:bg-ink-950/50"
                  >
                    <span className="text-xs font-semibold text-ink-400">{i + 1}</span>
                    <span className="shrink-0 text-sm font-medium">{p.name}</span>
                    <span
                      className="seq min-w-0 flex-1 truncate text-xs text-ink-500 dark:text-ink-400"
                      title={m?.aa || m?.description || p.dna}
                    >
                      {m?.aa || (p.source === "custom" ? p.dna : "—")}
                    </span>
                    <span className="shrink-0 text-xs text-ink-400">{p.dna.length} nt</span>
                    <button
                      type="button"
                      onClick={() => toggleLock(p.id)}
                      disabled={!p.coding}
                      aria-pressed={p.locked}
                      title={
                        !p.coding
                          ? "Not protein: always kept exactly as it is"
                          : p.locked
                            ? "Kept exactly as it is. Click to allow synonymous recoding."
                            : "Protein-coding: synonymous codons may change. Click to keep it exactly."
                      }
                      className="shrink-0 disabled:cursor-default"
                    >
                      <Badge tone={p.locked ? "neutral" : "info"}>
                        {p.locked ? "🔒 kept exactly" : "recodable"}
                      </Badge>
                    </button>
                    <span className="flex shrink-0 gap-0.5">
                      <button
                        type="button"
                        className="btn btn-ghost px-1.5 py-0.5 text-xs"
                        onClick={() => move(i, -1)}
                        disabled={i === 0}
                        aria-label={`Move ${p.name} earlier`}
                      >
                        ↑
                      </button>
                      <button
                        type="button"
                        className="btn btn-ghost px-1.5 py-0.5 text-xs"
                        onClick={() => move(i, 1)}
                        disabled={i === parts.length - 1}
                        aria-label={`Move ${p.name} later`}
                      >
                        ↓
                      </button>
                      <button
                        type="button"
                        className="btn btn-ghost px-1.5 py-0.5 text-xs hover:text-rose-600"
                        onClick={() => remove(p.id)}
                        aria-label={`Remove ${p.name}`}
                      >
                        ✕
                      </button>
                    </span>
                    {shift !== undefined && (
                      <span className="w-full pl-6">
                        <Badge tone="warn">
                          {p.dna.length} nt, not a multiple of 3: protein-coding parts after it are read {shift}{" "}
                          base{shift === 1 ? "" : "s"} off-frame
                        </Badge>
                      </span>
                    )}
                  </li>
                );
              })}
            </ul>
          )}
        </>
      )}

      <Field
        label="Cargo DNA"
        htmlFor="cargo"
        hint={
          mode === "builder"
            ? "Assembled from the parts above. Switch to “Paste DNA” to edit directly."
            : undefined
        }
      >
        <textarea
          id="cargo"
          rows={3}
          className="input seq text-xs"
          readOnly={mode === "builder"}
          value={cargo}
          onChange={(e) => onRawCargoChange(cleanDna(e.target.value))}
          placeholder="GACTACAAGGACGACGATGACAAG"
        />
      </Field>

      {mode === "raw" && (
        <Field
          label="Never change these nucleotides (optional)"
          htmlFor="cargo-excluded"
          hint="1-based, inclusive, comma-separated: 478-501, 12-20, 7. Without this, only the ORF up to the first stop codon is ever recoded, and a T7 promoter is never touched."
        >
          <input
            id="cargo-excluded"
            className="input seq text-xs"
            value={rawExcluded}
            onChange={(e) => onRawExcludedChange(e.target.value)}
            placeholder="478-501"
            aria-invalid={Boolean(excluded.error)}
          />
          {excluded.error && (
            <p role="alert" className="mt-1 text-xs font-medium text-rose-600">
              {excluded.error}
            </p>
          )}
        </Field>
      )}

      {cargo.length > 0 && (
        <div className="flex flex-wrap items-center gap-2">
          <Badge tone="neutral">{cargo.length} nt</Badge>
          {mode === "builder" && parts.length > 0 && (
            <Badge tone="neutral">
              {lockedCount} of {parts.length} part{parts.length === 1 ? "" : "s"} kept exactly
            </Badge>
          )}
          {frameOffset === 0 ? (
            <Badge tone="ok">Multiple of 3 — stays in frame</Badge>
          ) : (
            // Not automatically wrong: a cassette that ends in its own stop codon or a promoter does not
            // need to be a multiple of three. Flag it, but do not assert a frameshift.
            <Badge tone="warn">
              Not a multiple of 3 ({frameOffset} over) — fine at the end of the cargo or after a stop
              codon, otherwise it shifts the reading frame
            </Badge>
          )}
          <span className="ml-auto">
            <CopyButton value={cargo} label="Copy cargo" />
          </span>
        </div>
      )}
    </div>
  );
}
