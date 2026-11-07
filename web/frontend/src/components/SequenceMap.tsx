import { useState } from "react";
import type { ComponentRow } from "../types";
import { CopyButton, cx } from "./Primitives";

/**
 * The pegRNA drawn as a proportional, colour-coded map of its parts.
 *
 * The report gives you a 232-character string and a table of components; working
 * out which stretch is which meant counting bases by hand. This renders the parts
 * to scale, in construct order, and lets you click one to read just that
 * sequence. The colours match the component table below it.
 */

const PART_STYLES: Record<string, { bar: string; dot: string; text: string }> = {
  spacer: {
    bar: "bg-[var(--color-part-spacer)]",
    dot: "bg-[var(--color-part-spacer)]",
    text: "text-[var(--color-part-spacer)]",
  },
  scaffold: {
    bar: "bg-[var(--color-part-scaffold)]",
    dot: "bg-[var(--color-part-scaffold)]",
    text: "text-[var(--color-part-scaffold)]",
  },
  rcpbs2: {
    bar: "bg-[var(--color-part-pbs2)]",
    dot: "bg-[var(--color-part-pbs2)]",
    text: "text-[var(--color-part-pbs2)]",
  },
  insert: {
    bar: "bg-[var(--color-part-insert)]",
    dot: "bg-[var(--color-part-insert)]",
    text: "text-[var(--color-part-insert)]",
  },
  replacement_gdna: {
    bar: "bg-[var(--color-part-recode)]",
    dot: "bg-[var(--color-part-recode)]",
    text: "text-[var(--color-part-recode)]",
  },
  pbs1: {
    bar: "bg-[var(--color-part-pbs1)]",
    dot: "bg-[var(--color-part-pbs1)]",
    text: "text-[var(--color-part-pbs1)]",
  },
  polyu: {
    bar: "bg-[var(--color-part-polyu)]",
    dot: "bg-[var(--color-part-polyu)]",
    text: "text-[var(--color-part-polyu)]",
  },
  // Prime Assembly. The RTT is the flap template, the PBS is tj-PE's PBS1, and the donor is two arms around the
  // cargo: the arms are restored genomic DNA (the same colour as tj-PE's recoded stretch), the cargo is the insert.
  rtt: {
    bar: "bg-[var(--color-part-pbs2)]",
    dot: "bg-[var(--color-part-pbs2)]",
    text: "text-[var(--color-part-pbs2)]",
  },
  pbs: {
    bar: "bg-[var(--color-part-pbs1)]",
    dot: "bg-[var(--color-part-pbs1)]",
    text: "text-[var(--color-part-pbs1)]",
  },
  arm5: {
    bar: "bg-[var(--color-part-recode)]",
    dot: "bg-[var(--color-part-recode)]",
    text: "text-[var(--color-part-recode)]",
  },
  arm3: {
    bar: "bg-[var(--color-part-recode)]",
    dot: "bg-[var(--color-part-recode)]",
    text: "text-[var(--color-part-recode)]",
  },
  cargo: {
    bar: "bg-[var(--color-part-insert)]",
    dot: "bg-[var(--color-part-insert)]",
    text: "text-[var(--color-part-insert)]",
  },
};

const FALLBACK = {
  bar: "bg-ink-400",
  dot: "bg-ink-400",
  text: "text-ink-500",
};

/**
 * The key into the part tables for a component name. tj-PE's names are the keys; Prime Assembly's name the guide
 * ("spacer A", "PBS B", "RTT A (flap A)"), which this strips.
 */
export function partKey(component: string | undefined): string {
  const name = (component ?? "").toLowerCase().trim();
  if (name in PART_STYLES) return name;
  return name.replace(/\s*\(.*\)$/, "").replace(/\s+[ab]$/, "");
}

export function partStyle(component: string | undefined) {
  return PART_STYLES[partKey(component)] ?? FALLBACK;
}

/** Plain-language one-liners for the construct parts. */
const PART_PLAIN: Record<string, string> = {
  spacer: "Targets Cas9 to your site and makes the nick.",
  scaffold: "The constant sgRNA body that nCas9 grips.",
  rcpbs2: "Captures the second strand so the edit becomes permanent.",
  insert: "Your cargo — the tag being knocked in.",
  replacement_gdna: "The recoded stretch that breaks homology with the genome.",
  pbs1: "Binds the nicked strand so the template can be copied.",
  polyu: "A short U tail that stabilises the pegRNA.",
  rtt: "The template for the flap: the reverse transcriptase copies it onto the nicked strand, and the donor's end must match what it writes.",
  pbs: "Binds the nicked strand so the template can be copied.",
  arm5: "Restored genomic DNA left of your insertion point, recoded where it codes for protein.",
  arm3: "Restored genomic DNA right of your insertion point, copied exactly where it does not code for protein.",
  cargo: "Your insert, used exactly as you gave it.",
};

/** The plain-language line for a part, for the Prime Assembly maps. */
export function partPlain(component: string | undefined): string | undefined {
  return PART_PLAIN[partKey(component)];
}

export function SequenceMap({ components }: { components: ComponentRow[] }) {
  const [selected, setSelected] = useState<string | null>(null);

  const parts = components.filter((c) => c.sequence && c.component);
  if (!parts.length) return null;

  const total = parts.reduce((sum, p) => sum + (p.sequence?.length ?? 0), 0);
  if (!total) return null;

  const active = parts.find((p) => p.component === selected) ?? null;

  return (
    <div className="space-y-3">
      <div className="flex flex-wrap items-baseline justify-between gap-2">
        <h4 className="text-sm font-semibold text-ink-900 dark:text-ink-50">
          pegRNA construct map
        </h4>
        <p className="text-xs text-ink-500 dark:text-ink-400">
          {total} nt, drawn to scale 5&prime;&rarr;3&prime;. Select a part to read it.
        </p>
      </div>

      {/* The map itself. Each part is a button so it is keyboard reachable. */}
      <div
        className="flex h-9 w-full overflow-hidden rounded-lg border border-ink-200 dark:border-ink-700"
        role="group"
        aria-label="pegRNA parts, to scale"
      >
        {parts.map((p) => {
          const len = p.sequence?.length ?? 0;
          const pct = (len / total) * 100;
          const style = partStyle(p.component);
          const isActive = p.component === selected;
          return (
            <button
              key={p.index ?? p.component}
              type="button"
              onClick={() => setSelected(isActive ? null : (p.component ?? null))}
              style={{ width: `${pct}%` }}
              title={`${p.component} — ${len} nt`}
              aria-label={`${p.component}, ${len} nucleotides`}
              aria-pressed={isActive}
              className={cx(
                "group relative min-w-[2px] cursor-pointer border-r border-white/25 transition-opacity last:border-r-0 dark:border-black/25",
                style.bar,
                selected && !isActive ? "opacity-35" : "opacity-100",
              )}
            >
              {pct > 7 && (
                <span className="pointer-events-none absolute inset-0 flex items-center justify-center truncate px-1 text-[10px] font-semibold text-white drop-shadow-sm">
                  {p.component}
                </span>
              )}
            </button>
          );
        })}
      </div>

      {/* Legend doubles as the selector on narrow screens. */}
      <ul className="flex flex-wrap gap-x-4 gap-y-1.5">
        {parts.map((p) => {
          const style = partStyle(p.component);
          const isActive = p.component === selected;
          return (
            <li key={`legend-${p.index ?? p.component}`}>
              <button
                type="button"
                onClick={() => setSelected(isActive ? null : (p.component ?? null))}
                className={cx(
                  "flex items-center gap-1.5 text-xs transition-colors",
                  isActive
                    ? "font-semibold text-ink-900 dark:text-ink-50"
                    : "text-ink-600 hover:text-ink-900 dark:text-ink-400 dark:hover:text-ink-100",
                )}
              >
                <span aria-hidden className={cx("h-2.5 w-2.5 rounded-sm", style.dot)} />
                {p.component}
                <span className="text-ink-400">{p.sequence?.length}</span>
              </button>
            </li>
          );
        })}
      </ul>

      {active && (
        <div className="rounded-lg border border-ink-200 bg-ink-50 p-3.5 dark:border-ink-700 dark:bg-ink-950/60">
          <div className="mb-1.5 flex flex-wrap items-center justify-between gap-2">
            <div className="flex items-center gap-2">
              <span
                aria-hidden
                className={cx("h-2.5 w-2.5 rounded-sm", partStyle(active.component).dot)}
              />
              <span className="text-sm font-semibold">{active.component}</span>
              <span className="text-xs text-ink-500">{active.sequence?.length} nt</span>
            </div>
            <CopyButton value={active.sequence ?? ""} label="Copy sequence" />
          </div>
          <p className="mb-2 text-xs leading-relaxed text-ink-600 dark:text-ink-400">
            {partPlain(active.component) ?? active.role}
          </p>
          <p className="seq rounded border border-ink-200 bg-white p-2.5 text-xs text-ink-800 dark:border-ink-700 dark:bg-ink-900 dark:text-ink-100">
            {active.sequence}
          </p>
        </div>
      )}
    </div>
  );
}
