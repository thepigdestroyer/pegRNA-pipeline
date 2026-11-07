import { runFileURL } from "../api";
import type { ChangeRow, ChangesPayload, ComponentRow, DesignResult, Strategy } from "../types";
import { isStandingCaveat } from "../pa";
import { PaFailure, PaResult } from "./PaResult";
import { Badge, Callout, CopyButton, Disclosure, cx } from "./Primitives";
import { SequenceMap, partStyle } from "./SequenceMap";

function ComponentTable({ components }: { components: ComponentRow[] }) {
  if (!components.length) return null;
  const cols = Object.keys(components[0]);
  const numeric = new Set(["index", "length"]);
  return (
    <div className="overflow-x-auto">
      <table className="w-full table-fixed border-collapse text-xs">
        <colgroup>
          {cols.map((c) => (
            <col
              key={c}
              style={{
                width:
                  { index: "48px", component: "132px", role: "150px", part_of: "118px", length: "56px" }[
                    c
                  ] ?? "auto",
              }}
            />
          ))}
        </colgroup>
        <thead>
          <tr className="border-b border-ink-200 dark:border-ink-700">
            {cols.map((c) => (
              <th
                key={c}
                scope="col"
                className={cx(
                  "px-2 py-2 text-left font-semibold text-ink-600 dark:text-ink-300",
                  numeric.has(c) && "text-right",
                )}
              >
                {c.replace(/_/g, " ")}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {components.map((row, i) => (
            <tr
              key={row.index ?? i}
              className="border-b border-ink-100 last:border-0 dark:border-ink-800"
            >
              {cols.map((c) => (
                <td
                  key={c}
                  className={cx(
                    "px-2 py-1.5 align-top break-words",
                    c === "sequence" && "seq text-ink-700 dark:text-ink-300",
                    numeric.has(c) && "text-right tabular-nums",
                  )}
                >
                  {c === "component" ? (
                    <span className="flex items-center gap-1.5">
                      <span
                        aria-hidden
                        className={cx("h-2 w-2 shrink-0 rounded-sm", partStyle(row[c]).dot)}
                      />
                      <span className="font-medium">{row[c]}</span>
                    </span>
                  ) : (
                    row[c]
                  )}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

const CHANGE_LABEL: Record<ChangeRow["scope"], string> = {
  cargo: "Cargo",
  "genomic arm": "Genomic arm",
};

/**
 * Exactly what differs from what was specified (ADR 0029): every changed nucleotide, which part it is in, and
 * what was kept exactly. The amino acid never changes; it is shown so that can be checked by eye.
 */
function ChangesPanel({
  changes,
  runId,
  strategy,
}: {
  changes: ChangesPayload;
  runId: string;
  strategy: Strategy;
}) {
  const isPa = strategy === "pa";
  const armKept = changes.arm_kept_exactly ?? [];
  const rows = [...changes.cargo_changes, ...changes.arm_changes];
  const shown = rows.slice(0, 250);
  const summary = changes.cargo_unchanged
    ? `cargo as given · ${changes.arm_changes.length} arm nt`
    : `${changes.cargo_changes.length} cargo nt · ${changes.arm_changes.length} arm nt`;
  return (
    <Disclosure title="Changes from your input" description={summary} defaultOpen>
      <div className="space-y-3 text-sm">
        <p>
          <span className="font-semibold">Cargo: </span>
          {changes.cargo_unchanged ? (
            <>
              used exactly as you gave it ({changes.cargo_length} nt, <strong>0 changed</strong>).
            </>
          ) : (
            <>
              <strong>recoded, not as you gave it</strong>: {changes.cargo_changes.length} nucleotide
              {changes.cargo_changes.length === 1 ? "" : "s"} changed (variant {changes.cargo_variant}).
              Every amino acid is identical.
            </>
          )}
        </p>
        <p>
          <span className="font-semibold">{isPa ? "Donor arms: " : "Genomic arm: "}</span>
          {changes.arm_changes.length === 0
            ? "copied exactly (no nucleotide changed)."
            : `${changes.arm_changes.length} nucleotide${changes.arm_changes.length === 1 ? "" : "s"} of the design window changed by the silent recode. Every amino acid is identical.`}
        </p>

        {rows.length > 0 && (
          <div className="overflow-x-auto">
            <table className="w-full border-collapse text-xs">
              <thead>
                <tr className="border-b border-ink-200 dark:border-ink-700">
                  {["Where", "Part", "nt", "Codon", "Change", "aa"].map((h) => (
                    <th
                      key={h}
                      scope="col"
                      className="px-2 py-2 text-left font-semibold text-ink-600 dark:text-ink-300"
                    >
                      {h}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {shown.map((c, i) => (
                  <tr
                    key={`${c.scope}-${c.position}-${i}`}
                    className="border-b border-ink-100 last:border-0 dark:border-ink-800"
                  >
                    <td className="px-2 py-1.5">{CHANGE_LABEL[c.scope]}</td>
                    <td className="px-2 py-1.5">{c.part}</td>
                    <td className="px-2 py-1.5 tabular-nums">{c.position}</td>
                    <td className="px-2 py-1.5 tabular-nums">{c.codon ?? "—"}</td>
                    <td className="seq px-2 py-1.5">
                      {c.old_codon && c.new_codon
                        ? `${c.old_codon} → ${c.new_codon}`
                        : `${c.old_base} → ${c.new_base}`}
                    </td>
                    <td className="px-2 py-1.5">{c.amino_acid || "—"}</td>
                  </tr>
                ))}
              </tbody>
            </table>
            {rows.length > shown.length && (
              <p className="mt-2 text-xs text-ink-500 dark:text-ink-400">
                Showing {shown.length} of {rows.length}. The whole list is in the CSV.
              </p>
            )}
          </div>
        )}

        {changes.kept_exactly.length + armKept.length > 0 && (
          <div>
            <p className="font-semibold">Kept exactly as specified</p>
            <ul className="mt-1 list-disc space-y-0.5 pl-5 text-ink-700 dark:text-ink-300">
              {[...changes.kept_exactly, ...armKept].map((line, i) => (
                <li key={i}>{line}</li>
              ))}
            </ul>
          </div>
        )}
        {changes.notes.length > 0 && (
          <Callout tone="warn">
            <ul className="list-disc space-y-0.5 pl-5">
              {changes.notes.map((n, i) => (
                <li key={i}>{n}</li>
              ))}
            </ul>
          </Callout>
        )}

        <div className="flex flex-wrap items-center gap-3">
          <p className="text-xs text-ink-500 dark:text-ink-400">
            Everything not listed here is exactly as you specified.
          </p>
          <a
            className="btn btn-secondary ml-auto px-2.5 py-1 text-xs"
            href={runFileURL(runId, isPa ? "recoded_pa_changes.csv" : "recoded_tjpe_changes.csv")}
          >
            Changes (.csv)
          </a>
        </div>
      </div>
    </Disclosure>
  );
}

export function ResultPanel({ result }: { result: DesignResult }) {
  const standIn = result.gates?.stand_in_gates ?? [];
  const isDryRun = result.mode === "fake";
  // A server that predates Prime Assembly sends no strategy, and that means tj-PE.
  const strategy: Strategy = result.strategy ?? "tjpe";
  const isPa = strategy === "pa";
  const prefix = isPa ? "recoded_pa" : "recoded_tjpe";
  const pa = result.gates?.pa;
  const warnings = result.gates?.warnings ?? [];
  // The one standing caveat (the combination is untested) is not a finding about this design: it gets its own place.
  const caveat = isPa ? warnings.find(isStandingCaveat) : undefined;
  const findings = isPa ? warnings.filter((w) => !isStandingCaveat(w)) : warnings;

  return (
    <div className="card overflow-hidden">
      <header className="flex flex-wrap items-center gap-3 border-b border-ink-200 px-5 py-4 dark:border-ink-800">
        <h2 className="text-base font-semibold text-ink-900 dark:text-ink-50">Result</h2>
        {result.ok ? <Badge tone="ok">Design found</Badge> : <Badge tone="error">No design</Badge>}
        <Badge tone="neutral">{isPa ? "Prime Assembly" : "Inside the pegRNA"}</Badge>
        <span className="text-xs text-ink-500 dark:text-ink-400">
          run {result.run_id} · {result.mode} mode
        </span>
        {result.has_csv && (
          <span className="ml-auto flex flex-wrap gap-2">
            <a
              className="btn btn-secondary px-2.5 py-1 text-xs"
              href={runFileURL(result.run_id, `${prefix}_report.txt`)}
            >
              Report (.txt)
            </a>
            <a
              className="btn btn-secondary px-2.5 py-1 text-xs"
              href={runFileURL(result.run_id, `${prefix}_components.csv`)}
            >
              Components (.csv)
            </a>
            {isPa && result.has_fasta && (
              <a
                className="btn btn-secondary px-2.5 py-1 text-xs"
                href={runFileURL(result.run_id, "recoded_pa_donor.fasta")}
              >
                Sequences (.fasta)
              </a>
            )}
            <CopyButton value={result.report} label="Copy report" />
          </span>
        )}
      </header>

      <div className="space-y-5 p-5">
        {/*
          The headline safety message. A stand-in gate means part of this result
          is a placeholder, and that has to be the first thing read — not a line
          buried in the provenance block.
        */}
        {standIn.length > 0 && (
          <Callout
            tone={isDryRun ? "info" : "warn"}
            title={
              isDryRun
                ? "Dry run — this is not a validated design"
                : "Some checks did not run a real model"
            }
          >
            {isDryRun ? (
              <>
                <span className="font-medium">fake</span> mode exercises the design logic with
                stand-in scorers, so every verdict here is a placeholder. Use it to check your
                inputs flow through, then re-run in <span className="font-medium">api</span> mode
                before trusting the result.
              </>
            ) : (
              <>
                These steps returned a placeholder rather than a real verdict:{" "}
                <span className="font-semibold">{standIn.join(", ")}</span>.
                {standIn.includes("docking") && (
                  <> Docking needs an NVIDIA API key on the server.</>
                )}
                {standIn.includes("offtarget") && (
                  <> Off-target scoring needs a genome FASTA for cas-offinder.</>
                )}{" "}
                Everything else ran for real. Do not order DNA off the placeholder steps.
              </>
            )}
          </Callout>
        )}

        {!result.ok && result.mode === "api" && (
          <Callout tone="neutral">
            In <span className="font-medium">api</span> mode, “no candidate survived” is a genuine
            verdict from the checks — not a crash. The blocking step is marked in the panel above.
            Results are reproducible: the same inputs give the same verdict every time.
          </Callout>
        )}

        {!result.ok && standIn.length > 0 && (
          <Callout tone="warn" title="Part of this rejection came from a placeholder">
            These checks did not run a real model:{" "}
            <span className="font-semibold">{standIn.join(", ")}</span>. A rejection by a stand-in
            is not evidence of anything — fix the configuration above before concluding the design
            is bad.
          </Callout>
        )}

        {result.error && (
          <Callout tone="error" title="The run could not complete">
            {result.error}
          </Callout>
        )}

        {!result.ok && isPa && result.gates?.failed && <PaFailure gates={result.gates} />}

        {result.ok && caveat && (
          <Callout tone="warn" title="This combination is untested">
            {caveat}
          </Callout>
        )}

        {/*
          Advisories the design carries (a check that could only report, a tolerated
          off-target, a guide sequence worth a second look). They were only in the collapsed
          text report, so a reader ordering DNA could miss them.
        */}
        {result.ok && findings.length > 0 && (
          <Callout tone="warn" title="Check these before you order">
            <ul className="mt-1 list-disc space-y-1 pl-5">
              {findings.map((w, i) => (
                <li key={i}>{w}</li>
              ))}
            </ul>
          </Callout>
        )}

        {result.ok && isPa && pa && <PaResult result={result} pa={pa} />}

        {result.ok && result.gates?.changes && (
          <ChangesPanel changes={result.gates.changes} runId={result.run_id} strategy={strategy} />
        )}

        {result.ok && !isPa && result.components.length > 0 && (
          <SequenceMap components={result.components} />
        )}

        {result.gates?.gate_results && result.gates.gate_results.length > 0 && (
          <dl className="grid grid-cols-2 gap-3 sm:grid-cols-4">
            {(isPa
              ? [
                  ["Guide pairs considered", result.gates.pairs_considered],
                  ["Pairs tried", result.gates.candidates_tried],
                  ["Recode edits", result.gates.recode_edits],
                  ["Deletion size", `${result.gates.deletion_size ?? 0} bp`],
                ]
              : [
                  ["Candidates tried", result.gates.candidates_tried],
                  ["Recode edits", result.gates.recode_edits],
                  ["Deletion size", `${result.gates.deletion_size ?? 0} bp`],
                  ["Checks passed", `${result.gates.gate_results.length}`],
                ]
            ).map(([label, value]) => (
              <div
                key={String(label)}
                className="rounded-lg border border-ink-200 bg-ink-50 px-3 py-2 dark:border-ink-800 dark:bg-ink-950/50"
              >
                <dt className="text-xs text-ink-500 dark:text-ink-400">{label}</dt>
                <dd className="mt-0.5 text-lg font-semibold tabular-nums">{String(value)}</dd>
              </div>
            ))}
          </dl>
        )}

        {result.components.length > 0 && (
          <Disclosure title="All components" description={`${result.components.length} rows`}>
            <ComponentTable components={result.components} />
          </Disclosure>
        )}

        {result.report && (
          <Disclosure title="Full text report">
            <pre className="seq max-h-[28rem] overflow-auto rounded-lg bg-ink-950 p-4 text-xs leading-relaxed whitespace-pre-wrap text-ink-100">
              {result.report}
            </pre>
          </Disclosure>
        )}

        {result.het_edit?.requested && (
          <Disclosure
            title="Allele-specific pass"
            defaultOpen={result.het_edit.produced}
            description={result.het_edit.produced ? "pair produced" : "primary design only"}
          >
            <p className="mb-3 text-sm text-ink-600 dark:text-ink-400">{result.het_edit.note}</p>
            {result.het_edit.produced && result.het_edit.components.length > 0 && (
              <ComponentTable components={result.het_edit.components} />
            )}
          </Disclosure>
        )}

        {result.stderr && (
          <Disclosure title="Diagnostics" defaultOpen={!result.ok}>
            <pre className="seq max-h-80 overflow-auto rounded-lg bg-ink-950 p-4 text-xs leading-relaxed whitespace-pre-wrap text-ink-100">
              {result.stderr}
            </pre>
          </Disclosure>
        )}
      </div>
    </div>
  );
}
