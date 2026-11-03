import { runFileURL } from "../api";
import type { ComponentRow, DesignResult } from "../types";
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

export function ResultPanel({ result }: { result: DesignResult }) {
  const standIn = result.gates?.stand_in_gates ?? [];
  const isDryRun = result.mode === "fake";

  return (
    <div className="card overflow-hidden">
      <header className="flex flex-wrap items-center gap-3 border-b border-ink-200 px-5 py-4 dark:border-ink-800">
        <h2 className="text-base font-semibold text-ink-900 dark:text-ink-50">Result</h2>
        {result.ok ? <Badge tone="ok">Design found</Badge> : <Badge tone="error">No design</Badge>}
        <span className="text-xs text-ink-500 dark:text-ink-400">
          run {result.run_id} · {result.mode} mode
        </span>
        {result.has_csv && (
          <span className="ml-auto flex flex-wrap gap-2">
            <a
              className="btn btn-secondary px-2.5 py-1 text-xs"
              href={runFileURL(result.run_id, "recoded_tjpe_report.txt")}
            >
              Report (.txt)
            </a>
            <a
              className="btn btn-secondary px-2.5 py-1 text-xs"
              href={runFileURL(result.run_id, "recoded_tjpe_components.csv")}
            >
              Components (.csv)
            </a>
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
            Docking in particular is stochastic, so the same inputs can pass on another run.
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

        {result.ok && result.components.length > 0 && (
          <SequenceMap components={result.components} />
        )}

        {result.gates?.gate_results && result.gates.gate_results.length > 0 && (
          <dl className="grid grid-cols-2 gap-3 sm:grid-cols-4">
            {[
              ["Candidates tried", result.gates.candidates_tried],
              ["Recode edits", result.gates.recode_edits],
              ["Deletion size", `${result.gates.deletion_size ?? 0} bp`],
              ["Checks passed", `${result.gates.gate_results.length}`],
            ].map(([label, value]) => (
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
