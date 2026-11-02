import { useMemo, useState } from "react";
import type { Motif, MotifGroups } from "../types";
import { Badge, CopyButton, Field, cx } from "./Primitives";

/**
 * Build the cargo from the motif library, or type raw DNA.
 *
 * The library is the safe path: parts are assembled 5'->3' and the running
 * frame check catches the mistake that actually bites — a cargo whose length is
 * not a multiple of three, which silently frameshifts everything downstream.
 */

interface Props {
  motifs: MotifGroups;
  parts: string[];
  onPartsChange: (parts: string[]) => void;
  rawCargo: string;
  onRawCargoChange: (dna: string) => void;
  mode: "builder" | "raw";
  onModeChange: (mode: "builder" | "raw") => void;
}

export function CargoBuilder({
  motifs,
  parts,
  onPartsChange,
  rawCargo,
  onRawCargoChange,
  mode,
  onModeChange,
}: Props) {
  const [pick, setPick] = useState("");

  const byName = useMemo(() => {
    const m = new Map<string, Motif>();
    Object.values(motifs).forEach((list) => list.forEach((x) => m.set(x.name, x)));
    return m;
  }, [motifs]);

  const assembled = useMemo(
    () => parts.map((n) => byName.get(n)?.dna ?? "").join(""),
    [parts, byName],
  );

  const cargo = mode === "raw" ? rawCargo : assembled;
  const frameOffset = cargo.length % 3;

  function add() {
    if (pick) onPartsChange([...parts, pick]);
  }
  function remove(i: number) {
    onPartsChange(parts.filter((_, idx) => idx !== i));
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
              {m === "builder" ? "Pick from library" : "Paste DNA"}
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
              <option value="">Choose a part…</option>
              {Object.keys(motifs)
                .sort()
                .map((cat) => (
                  <optgroup key={cat} label={cat.replace(/_/g, " ")}>
                    {motifs[cat].map((m) => (
                      <option key={m.name} value={m.name}>
                        {m.name}
                        {m.aa ? ` (${m.aa.length > 16 ? `${m.aa.slice(0, 16)}…` : m.aa})` : ""}
                      </option>
                    ))}
                  </optgroup>
                ))}
            </select>
            <button type="button" className="btn btn-secondary" onClick={add} disabled={!pick}>
              Add part
            </button>
          </div>

          {parts.length === 0 ? (
            <p className="rounded-lg border border-dashed border-ink-300 px-4 py-6 text-center text-sm text-ink-500 dark:border-ink-700 dark:text-ink-400">
              No parts yet. Add a tag, linker or marker above.
            </p>
          ) : (
            <ul className="space-y-1.5">
              {parts.map((name, i) => {
                const m = byName.get(name);
                return (
                  <li
                    key={`${name}-${i}`}
                    className="flex items-center gap-3 rounded-lg border border-ink-200 bg-ink-50
                               px-3 py-2 dark:border-ink-800 dark:bg-ink-950/50"
                  >
                    <span className="text-xs font-semibold text-ink-400">{i + 1}</span>
                    <span className="shrink-0 text-sm font-medium">{name}</span>
                    <span
                      className="seq min-w-0 flex-1 truncate text-xs text-ink-500 dark:text-ink-400"
                      title={m?.aa || m?.description}
                    >
                      {m?.aa || "—"}
                    </span>
                    <span className="shrink-0 text-xs text-ink-400">{m?.dna.length ?? 0} nt</span>
                    <span className="flex shrink-0 gap-0.5">
                      <button
                        type="button"
                        className="btn btn-ghost px-1.5 py-0.5 text-xs"
                        onClick={() => move(i, -1)}
                        disabled={i === 0}
                        aria-label={`Move ${name} earlier`}
                      >
                        ↑
                      </button>
                      <button
                        type="button"
                        className="btn btn-ghost px-1.5 py-0.5 text-xs"
                        onClick={() => move(i, 1)}
                        disabled={i === parts.length - 1}
                        aria-label={`Move ${name} later`}
                      >
                        ↓
                      </button>
                      <button
                        type="button"
                        className="btn btn-ghost px-1.5 py-0.5 text-xs hover:text-rose-600"
                        onClick={() => remove(i)}
                        aria-label={`Remove ${name}`}
                      >
                        ✕
                      </button>
                    </span>
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
          onChange={(e) => onRawCargoChange(e.target.value.replace(/\s/g, "").toUpperCase())}
          placeholder="GACTACAAGGACGACGATGACAAG"
        />
      </Field>

      {cargo.length > 0 && (
        <div className="flex flex-wrap items-center gap-2">
          <Badge tone="neutral">{cargo.length} nt</Badge>
          {frameOffset === 0 ? (
            <Badge tone="ok">Multiple of 3 — stays in frame</Badge>
          ) : (
            // Not automatically wrong: a cassette that ends in its own stop codon
            // (the shipped FLAG-P2A-ZeoR example is 503 nt) does not need to be a
            // multiple of three. Flag it, but do not assert a frameshift.
            <Badge tone="warn">
              Not a multiple of 3 ({frameOffset} over) — fine if the cargo ends in a stop codon,
              otherwise it will shift the reading frame
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
