import { Fragment, type ReactNode, useState } from "react";
import { runFileURL } from "../api";
import { pamFate } from "../pa";
import { PA_PIPELINE, failureCounts, findBlockingStep } from "../pipeline";
import type { ComponentRow, DesignResult, GatesPayload, PaGuide, PaPayload } from "../types";
import { Badge, Callout, CopyButton, Disclosure, cx } from "./Primitives";
import { partPlain, partStyle } from "./SequenceMap";

const asRna = (dna: string) => dna.replace(/T/g, "U");

interface MapPart {
  /** What the colour and the one-line explanation are looked up by. */
  key: string;
  /** What the reader sees. */
  name: string;
  sequence: string;
}

interface Bracket {
  from: number;
  to: number;
  label: string;
  side: "left" | "right";
}

/**
 * A construct drawn to scale as a bar of its parts, in order, 5'->3', with an optional bracket under it for a
 * stretch that means something (the ends of the donor that the flaps match). Select a part to read and copy it.
 */
function PartsMap({
  label,
  parts,
  brackets = [],
}: {
  label: string;
  parts: MapPart[];
  brackets?: Bracket[];
}) {
  const [selected, setSelected] = useState<string | null>(null);
  const total = parts.reduce((n, p) => n + p.sequence.length, 0);
  if (!total) return null;
  const active = parts.find((p) => p.key === selected) ?? null;
  const toggle = (key: string) => setSelected(selected === key ? null : key);

  return (
    <div className="space-y-2">
      <div
        className="flex h-8 w-full overflow-hidden rounded-lg border border-ink-200 dark:border-ink-700"
        role="group"
        aria-label={label}
      >
        {parts.map((p) => {
          const pct = (p.sequence.length / total) * 100;
          const isActive = p.key === selected;
          return (
            <button
              key={p.key}
              type="button"
              onClick={() => toggle(p.key)}
              style={{ width: `${pct}%` }}
              title={`${p.name}, ${p.sequence.length} nt`}
              aria-label={`${p.name}, ${p.sequence.length} nucleotides`}
              aria-pressed={isActive}
              className={cx(
                "relative min-w-[2px] cursor-pointer border-r border-white/25 transition-opacity last:border-r-0 dark:border-black/25",
                partStyle(p.key).bar,
                selected && !isActive ? "opacity-35" : "opacity-100",
              )}
            >
              {pct > 6 && (
                <span className="pointer-events-none absolute inset-0 block truncate px-1 text-center text-[10px] leading-8 font-semibold text-white drop-shadow-sm">
                  {p.name}
                </span>
              )}
            </button>
          );
        })}
      </div>

      {brackets.length > 0 && (
        <div className="relative h-7" aria-hidden>
          {brackets.map((b) => (
            <Fragment key={b.label}>
              <div
                className="absolute top-0 h-2 rounded-b border-x-2 border-b-2 border-ink-500 dark:border-ink-400"
                style={{ left: `${(b.from / total) * 100}%`, width: `${((b.to - b.from) / total) * 100}%` }}
              />
              <span
                className="absolute top-3 text-[11px] whitespace-nowrap text-ink-600 dark:text-ink-300"
                style={
                  b.side === "left"
                    ? { left: `${(b.from / total) * 100}%` }
                    : { right: `${100 - (b.to / total) * 100}%` }
                }
              >
                {b.label}
              </span>
            </Fragment>
          ))}
        </div>
      )}

      <ul className="flex flex-wrap gap-x-4 gap-y-1.5">
        {parts.map((p) => (
          <li key={p.key}>
            <button
              type="button"
              onClick={() => toggle(p.key)}
              className={cx(
                "flex items-center gap-1.5 text-xs transition-colors",
                p.key === selected
                  ? "font-semibold text-ink-900 dark:text-ink-50"
                  : "text-ink-600 hover:text-ink-900 dark:text-ink-400 dark:hover:text-ink-100",
              )}
            >
              <span aria-hidden className={cx("h-2.5 w-2.5 rounded-sm", partStyle(p.key).dot)} />
              {p.name}
              <span className="text-ink-400">{p.sequence.length}</span>
            </button>
          </li>
        ))}
      </ul>

      {active && (
        <div className="rounded-lg border border-ink-200 bg-ink-50 p-3 dark:border-ink-700 dark:bg-ink-950/60">
          <div className="mb-1.5 flex flex-wrap items-center justify-between gap-2">
            <span className="text-sm font-semibold">
              {active.name} <span className="text-xs font-normal text-ink-500">{active.sequence.length} nt</span>
            </span>
            <CopyButton value={active.sequence} label="Copy sequence" />
          </div>
          {partPlain(active.key) && (
            <p className="mb-2 text-xs leading-relaxed text-ink-600 dark:text-ink-400">{partPlain(active.key)}</p>
          )}
          <p className="seq rounded border border-ink-200 bg-white p-2.5 text-xs break-all text-ink-800 dark:border-ink-700 dark:bg-ink-900 dark:text-ink-100">
            {active.sequence}
          </p>
        </div>
      )}
    </div>
  );
}

function Fact({ label, children }: { label: string; children: ReactNode }) {
  return (
    <div className="min-w-0">
      <dt className="text-xs text-ink-500 dark:text-ink-400">{label}</dt>
      <dd className="mt-0.5 text-sm break-words text-ink-800 dark:text-ink-100">{children}</dd>
    </div>
  );
}

function GuideCard({ guide, components }: { guide: PaGuide; components: ComponentRow[] }) {
  const parts: MapPart[] = components
    .filter((c) => c.part_of === `pegRNA ${guide.side}` && c.sequence)
    .map((c, i) => ({ key: `${c.component}`, name: (c.component ?? `part ${i + 1}`).replace(/\s*\(.*\)$/, ""), sequence: c.sequence ?? "" }));
  const fate = pamFate(guide);
  const end = guide.side === "A" ? "5′" : "3′";
  return (
    <div className="space-y-3 rounded-xl border border-ink-200 p-4 dark:border-ink-800">
      <header className="flex flex-wrap items-baseline gap-x-3 gap-y-1">
        <h4 className="text-sm font-semibold text-ink-900 dark:text-ink-50">pegRNA {guide.side}</h4>
        <span className="text-xs text-ink-500 dark:text-ink-400">
          {guide.strand} strand · nick after nt {guide.cut} · {guide.pegRNA_length} nt
        </span>
        <span className="ml-auto flex gap-1">
          <CopyButton value={guide.pegRNA} label="Copy" />
          <CopyButton value={asRna(guide.pegRNA)} label="Copy as RNA" />
        </span>
      </header>
      <p className="seq rounded border border-ink-200 bg-ink-50 p-2.5 text-xs break-all text-ink-800 dark:border-ink-700 dark:bg-ink-950/60 dark:text-ink-100">
        {guide.pegRNA}
      </p>
      <p className="text-xs text-ink-500 dark:text-ink-400">
        A pegRNA is RNA: read each T as U (the second button copies it that way).
      </p>
      <PartsMap label={`pegRNA ${guide.side}, parts to scale`} parts={parts} />

      <dl className="grid gap-x-6 gap-y-3 sm:grid-cols-2">
        <Fact label="Spacer">
          <span className="seq text-xs">{guide.spacer}</span>
        </Fact>
        <Fact label="Template and primer-binding site">
          {guide.flap_len} nt template writing the donor&rsquo;s {end} end · {guide.pbs_len} nt PBS
        </Fact>
        <Fact label="PAM">
          <span className="seq text-xs">{guide.pam}</span> in the genome,{" "}
          <span className="seq text-xs">{guide.product_pam}</span> in the product
        </Fact>
        <Fact label="Off-target hits">
          {guide.offtargets === 0
            ? "none within the screened tier"
            : `${guide.offtargets} (tolerated by your setting: check before ordering)`}
        </Fact>
      </dl>

      <div
        className={cx(
          "rounded-lg border px-3 py-2 text-xs leading-relaxed",
          fate.tone === "ok"
            ? "border-emerald-200 bg-emerald-50 text-emerald-800 dark:border-emerald-900 dark:bg-emerald-950 dark:text-emerald-200"
            : "border-amber-200 bg-amber-50 text-amber-900 dark:border-amber-900 dark:bg-amber-950 dark:text-amber-200",
        )}
      >
        <span className="font-semibold">Guide {guide.side}&rsquo;s PAM. </span>
        {fate.text}
      </div>
      {guide.nick_pairing_bp > 0 && (
        <p className="text-xs leading-relaxed text-ink-500 dark:text-ink-400">
          The end of the template pairs with the end of the spacer for {guide.nick_pairing_bp} bp, because the donor
          copies the genome at this nick. That is expected and was not counted as a misfold.
        </p>
      )}
    </div>
  );
}

function DonorCard({ pa }: { pa: PaPayload }) {
  const label: Record<string, string> = { arm5: "5′ arm", cargo: "cargo", arm3: "3′ arm" };
  const parts: MapPart[] = pa.donor.parts.map((p) => ({
    key: p.name,
    name: label[p.name] ?? p.name,
    sequence: p.sequence,
  }));
  const [a0, a1] = pa.donor.flap_a;
  const [b0, b1] = pa.donor.flap_b;
  return (
    <div className="space-y-3 rounded-xl border border-ink-200 p-4 dark:border-ink-800">
      <header className="flex flex-wrap items-baseline gap-x-3 gap-y-1">
        <h4 className="text-sm font-semibold text-ink-900 dark:text-ink-50">Donor DNA</h4>
        <span className="text-xs text-ink-500 dark:text-ink-400">
          {pa.donor_length} nt, double-stranded · 5′ arm + cargo + 3′ arm
        </span>
        <span className="ml-auto">
          <CopyButton value={pa.donor.sequence} label="Copy" />
        </span>
      </header>
      <p className="seq max-h-40 overflow-auto rounded border border-ink-200 bg-ink-50 p-2.5 text-xs break-all text-ink-800 dark:border-ink-700 dark:bg-ink-950/60 dark:text-ink-100">
        {pa.donor.sequence}
      </p>
      <PartsMap
        label="Donor DNA, parts to scale"
        parts={parts}
        brackets={[
          { from: a0, to: a1, label: `5′ ${a1 - a0} nt match flap A`, side: "left" },
          { from: b0, to: b1, label: `3′ ${b1 - b0} nt match flap B`, side: "right" },
        ]}
      />
      <p className="text-xs leading-relaxed text-ink-600 dark:text-ink-400">
        The two arms put back the DNA between the nicks, so nothing is deleted: where it codes for protein it is
        recoded with synonymous codons, and where it does not it is copied exactly. Your cargo, in the middle, is used
        exactly as you gave it.
      </p>
    </div>
  );
}

/** The genome window with the two nicks, the insertion point and the stretch the donor replaces and restores. */
function LocusMap({ pa }: { pa: PaPayload }) {
  const width = pa.window_length;
  const pct = (x: number) => `${(x / width) * 100}%`;
  const [n1, n2] = pa.span;
  const [guideA, guideB] = pa.guides;
  const mark = (x: number, text: string, where: "above" | "below") => (
    <div
      className="absolute w-0.5 -translate-x-1/2 bg-ink-800 dark:bg-ink-100"
      style={{ left: pct(x), top: where === "above" ? "-12px" : "0", height: where === "above" ? "28px" : "28px" }}
    >
      <span
        className={cx(
          "absolute left-1/2 -translate-x-1/2 text-[11px] font-semibold whitespace-nowrap text-ink-800 dark:text-ink-100",
          where === "above" ? "-top-4" : "top-full",
        )}
      >
        {text}
      </span>
    </div>
  );
  return (
    <div className="space-y-2">
      <div
        className="relative mx-3 my-8 h-4 rounded bg-ink-200 dark:bg-ink-800"
        role="img"
        aria-label={`The ${width} nt window. Nick A after nt ${n1}, the insertion point at ${pa.site_s}, nick B after nt ${n2}. The donor replaces nt ${n1 + 1} to ${n2} and puts it back with your insert at ${pa.site_s}.`}
      >
        <div
          className="absolute inset-y-0 bg-brand-300/70 ring-2 ring-brand-600 dark:bg-brand-500/50 dark:ring-brand-300"
          style={{ left: pct(n1), width: pct(n2 - n1) }}
        />
        {pa.recode_spans.map(([a, b]) => (
          <div
            key={`recode-${a}-${b}`}
            className="absolute inset-y-0 bg-[var(--color-part-recode)]"
            style={{ left: pct(a), width: pct(b - a) }}
          />
        ))}
        {mark(n1, "A", "above")}
        {mark(n2, "B", "above")}
        {mark(pa.site_s, "S", "below")}
      </div>
      <ul className="flex flex-wrap gap-x-5 gap-y-1 text-xs text-ink-600 dark:text-ink-400">
        <li>
          <span className="font-semibold">A</span> nick, {guideA.strand} strand, after nt {n1}
        </li>
        <li>
          <span className="font-semibold">B</span> nick, {guideB.strand} strand, after nt {n2}
        </li>
        <li>
          <span className="font-semibold">S</span> insertion point, {pa.site_s}
        </li>
        <li className="flex items-center gap-1.5">
          <span aria-hidden className="h-2.5 w-2.5 rounded-sm ring-2 ring-brand-600 dark:ring-brand-300" />
          between the nicks: replaced, then restored, by the donor ({n2 - n1} bp)
        </li>
        <li className="flex items-center gap-1.5">
          <span aria-hidden className="h-2.5 w-2.5 rounded-sm bg-[var(--color-part-recode)]" />
          recoded (synonymous), where it codes for protein
        </li>
        <li className="flex items-center gap-1.5">
          <span aria-hidden className="h-2.5 w-2.5 rounded-sm bg-brand-300/70 dark:bg-brand-500/50" />
          copied exactly, where it does not
        </li>
      </ul>
    </div>
  );
}

/** Everything specific to a Prime Assembly design: the order sheet, where it goes, and how to make the donor. */
export function PaResult({ result, pa }: { result: DesignResult; pa: PaPayload }) {
  return (
    <div className="space-y-5">
      <div className="space-y-3">
        <div className="flex flex-wrap items-center gap-3">
          <h3 className="text-sm font-semibold text-ink-900 dark:text-ink-50">What to order</h3>
          <Badge tone="neutral">2 pegRNAs + 1 donor</Badge>
          {result.has_fasta && (
            <a
              className="btn btn-secondary ml-auto px-2.5 py-1 text-xs"
              href={runFileURL(result.run_id, "recoded_pa_donor.fasta")}
            >
              All three (.fasta)
            </a>
          )}
        </div>
        <div className="grid gap-4 xl:grid-cols-2">
          {pa.guides.map((g) => (
            <GuideCard key={g.side} guide={g} components={result.components} />
          ))}
        </div>
        <DonorCard pa={pa} />
      </div>

      <div className="space-y-2">
        <h3 className="text-sm font-semibold text-ink-900 dark:text-ink-50">Where it goes</h3>
        <LocusMap pa={pa} />
      </div>

      {pa.donor_prep.length > 0 && (
        <Disclosure title="How to prepare the donor" description="from the published method, not yet from this design">
          <ul className="list-disc space-y-1.5 pl-5 text-sm leading-relaxed text-ink-700 dark:text-ink-300">
            {pa.donor_prep.map((line, i) => (
              <li key={i}>{line}</li>
            ))}
          </ul>
          <p className="text-xs text-ink-500 dark:text-ink-400">
            Source:{" "}
            <a
              className="underline underline-offset-2"
              href="https://pmc.ncbi.nlm.nih.gov/articles/PMC12262537/"
              rel="noreferrer"
              target="_blank"
            >
              Liu et al., Prime Assembly with linear DNA donors (preprint)
            </a>
            .
          </p>
        </Disclosure>
      )}
    </div>
  );
}


/** What to try next, by the step that stopped the furthest candidate. */
const NEXT_STEP: Record<string, string> = {
  importance:
    "The insertion site is in a part of the protein that an insert may disturb. Choose a surface, loop or terminal site.",
  offtarget:
    "Most guides in the window match elsewhere in the genome. Allowing one off-target hit per guide (Off-target hits allowed) admits more guide pairs; each tolerated hit is listed on the result.",
  pairs:
    "No two guides face each other across the insertion point within the allowed distance. Widen the distance between the nicks, or give a wider window around the insertion point.",
  rules:
    "Most guides in the window break a sequence rule: a run of four T, GC content outside 30 to 75%, or a primer-binding site that would run off the end of the window. A wider window around the insertion point offers more guides.",
  recode:
    "The donor's arms could not be recoded without changing the protein. Try other guide pairs (more pairs to try) or a different distance between the nicks.",
  inv0: "The donor's arms could not be built without a deletion or a protein change for the pairs tried. Try more pairs.",
  flaps:
    "No usable template for the pairs tried: either a run of four T, which stops RNA polymerase III, lies in the end of the donor at every flap length from 30 to 50 nt, or the donor's two ends share a long identical stretch, so the two flaps could pair with each other instead of with the donor. Another pair, a wider window, or a different end to your insert may work.",
  strand: "The strand check found the ordered pegRNAs would not rebuild the intended allele. This should not happen: please report it.",
  structure_a: "The first guide's pegRNAs fold back on themselves for the pairs tried. Try more PBS lengths, a different flap length, or more pairs.",
  structure_b: "The second guide's pegRNAs fold back on themselves for the pairs tried. Try more PBS lengths, a different flap length, or more pairs.",
  cryptic: "The edit would create a new splice site in the pairs tried. Try more pairs, or another insertion site.",
  splice: "SpliceAI judged the edit to disturb splicing in the pairs tried. Try more pairs, or another insertion site.",
  docking_a: "The first guide's pegRNAs bound Cas9 too weakly for the pairs tried. Try more pairs, other PBS lengths or another flap length.",
  docking_b: "The second guide's pegRNAs bound Cas9 too weakly for the pairs tried. Try more pairs, other PBS lengths or another flap length.",
};

/**
 * Why a Prime Assembly run found no design. Many guide pairs are tried, so there is no single reason: this counts the
 * rejections by step, says where the furthest candidate stopped and what could be tried, and keeps every line of the
 * pipeline's own record behind a disclosure.
 */
export function PaFailure({ gates }: { gates: GatesPayload }) {
  const diagnostics = gates.diagnostics ?? [];
  const counts = failureCounts(PA_PIPELINE, diagnostics);
  const blockIdx = findBlockingStep("", "", PA_PIPELINE, diagnostics);
  const next = blockIdx >= 0 ? NEXT_STEP[PA_PIPELINE[blockIdx].id] : undefined;
  return (
    <Callout tone="neutral" title="Why there is no design">
      {counts.length > 0 ? (
        <ul className="mt-1 list-disc space-y-0.5 pl-5">
          {counts.map((c) => (
            <li key={c.id}>
              <span className="font-medium">{c.title}:</span> {c.count} {c.unit}
              {c.count === 1 ? "" : "s"} rejected
            </li>
          ))}
        </ul>
      ) : (
        gates.error && <p className="seq text-xs break-words">{gates.error}</p>
      )}
      {next && <p className="mt-2">{next}</p>}
      {diagnostics.length > 0 && (
        <div className="mt-3">
          <Disclosure title="Every rejection, in the order the pipeline recorded it" description={`${diagnostics.length} lines`}>
            <pre className="seq max-h-72 overflow-auto rounded-lg bg-ink-950 p-3 text-xs leading-relaxed whitespace-pre-wrap text-ink-100">
              {diagnostics.join("\n")}
            </pre>
          </Disclosure>
        </div>
      )}
    </Callout>
  );
}
