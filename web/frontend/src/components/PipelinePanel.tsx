import type { DesignResult, GateMetric, ProgressEvent } from "../types";
import { COMPARISON_SYMBOL, PIPELINE, blockReason, findBlockingStep } from "../pipeline";
import { InfoHint, cx } from "./Primitives";

export type StepState = "idle" | "running" | "passed" | "blocked" | "visited" | "unreached";

const DOT: Record<StepState, string> = {
  idle: "border-ink-300 bg-white dark:border-ink-700 dark:bg-ink-900",
  running: "border-brand-500 bg-brand-500",
  passed: "border-emerald-500 bg-emerald-500",
  blocked: "border-rose-500 bg-rose-500",
  visited: "border-ink-300 bg-ink-300 dark:border-ink-600 dark:bg-ink-600",
  unreached: "border-ink-200 bg-white dark:border-ink-800 dark:bg-ink-900",
};

/** Render a gate's raw value against the cutoff it was judged on. */
function MetricReadout({ metric, info, passed }: { metric?: GateMetric; info: string; passed: boolean }) {
  if (metric && metric.value != null && metric.cutoff != null) {
    const symbol = COMPARISON_SYMBOL[metric.cmp ?? ""] ?? metric.cmp ?? "";
    return (
      <span className="flex flex-wrap items-baseline gap-x-1.5 gap-y-0.5">
        <span className="text-ink-500 dark:text-ink-400">{metric.metric}</span>
        <span className="seq font-semibold text-ink-900 dark:text-ink-50">
          {String(metric.value)}
        </span>
        <span className="text-ink-400">
          ({symbol}
          {String(metric.cutoff)})
        </span>
        <span
          className={cx(
            "rounded-full px-1.5 text-[10px] font-semibold",
            passed
              ? "bg-emerald-100 text-emerald-700 dark:bg-emerald-950 dark:text-emerald-300"
              : "bg-rose-100 text-rose-700 dark:bg-rose-950 dark:text-rose-300",
          )}
        >
          {passed ? "within cutoff" : "outside cutoff"}
        </span>
        {metric.extra && <span className="block w-full text-ink-400">{metric.extra}</span>}
      </span>
    );
  }
  if (metric?.note) return <span className="text-ink-500 dark:text-ink-400">{metric.note}</span>;
  return <span className="text-ink-500 dark:text-ink-400">{info}</span>;
}

interface Props {
  running: boolean;
  progress: Record<string, ProgressEvent>;
  result: DesignResult | null;
}

export function PipelinePanel({ running, progress, result }: Props) {
  const gateProse: Record<string, string> = {};
  (result?.gates?.gate_results ?? []).forEach(([id, info]) => {
    gateProse[id] = info;
  });

  const standIn = new Set(result?.gates?.stand_in_gates ?? []);
  const blockIdx = result && !result.ok ? findBlockingStep(result.stderr ?? "", result.report) : -1;
  const activeGate = running
    ? Object.values(progress).reduce<ProgressEvent | null>(
        (latest, ev) => (latest === null || ev.cand >= latest.cand ? ev : latest),
        null,
      )?.gate
    : undefined;

  function stateFor(index: number, id: string): StepState {
    if (running) {
      if (id === activeGate) return "running";
      return progress[id] ? (progress[id].passed ? "passed" : "blocked") : "idle";
    }
    if (!result) return "idle";
    if (result.ok) return "passed";
    if (blockIdx < 0) return "idle";
    if (index < blockIdx) return "passed";
    if (index === blockIdx) return "blocked";
    return "unreached";
  }

  return (
    <div className="card overflow-hidden">
      <header className="border-b border-ink-200 px-5 py-4 dark:border-ink-800">
        <h2 className="text-base font-semibold text-ink-900 dark:text-ink-50">
          How your design is checked
        </h2>
        <p className="mt-1 text-sm leading-relaxed text-ink-600 dark:text-ink-400">
          {running
            ? "Running now — each step shows its measured value as it completes."
            : result
              ? result.ok
                ? "Every check passed. The winning candidate is shown below."
                : "A check rejected every candidate. The blocking step is marked in red."
              : "Candidates run cheapest check first. The first one to pass every step wins."}
        </p>
      </header>

      <ol className="divide-y divide-ink-100 dark:divide-ink-800">
        {PIPELINE.map((step, i) => {
          const state = stateFor(i, step.id);
          const live = progress[step.id];
          const isStandIn = standIn.has(step.id);
          return (
            <li
              key={step.id}
              className={cx(
                "flex gap-3 px-5 py-3",
                state === "running" && "bg-brand-50/60 dark:bg-brand-700/10",
                state === "blocked" && "bg-rose-50/60 dark:bg-rose-950/20",
              )}
            >
              <div className="flex flex-col items-center">
                <span
                  aria-hidden
                  className={cx(
                    "mt-1 flex h-3.5 w-3.5 shrink-0 items-center justify-center rounded-full border-2",
                    DOT[state],
                    state === "running" && "animate-pulse",
                  )}
                />
                {i < PIPELINE.length - 1 && (
                  <span
                    aria-hidden
                    className="mt-1 w-px flex-1 bg-ink-200 dark:bg-ink-800"
                  />
                )}
              </div>

              <div className="min-w-0 flex-1">
                <div className="flex flex-wrap items-baseline gap-x-2 gap-y-1">
                  <span
                    className={cx(
                      "text-sm font-medium",
                      state === "blocked"
                        ? "text-rose-700 dark:text-rose-300"
                        : "text-ink-900 dark:text-ink-50",
                    )}
                  >
                    {step.title}
                  </span>
                  <span className="text-xs text-ink-400">{step.engine}</span>
                  <InfoHint label={step.title}>{step.why}</InfoHint>
                  {isStandIn && (
                    <span
                      className="rounded-full bg-amber-100 px-1.5 text-[10px] font-semibold text-amber-800
                                 dark:bg-amber-950 dark:text-amber-300"
                      title="This step did not run a real model — its verdict is a placeholder."
                    >
                      stand-in
                    </span>
                  )}
                  <span className="ml-auto text-xs font-semibold">
                    {state === "passed" && <span className="text-emerald-600">passed</span>}
                    {state === "blocked" && <span className="text-rose-600">blocked</span>}
                    {state === "running" && <span className="text-brand-600">running</span>}
                    {state === "unreached" && <span className="text-ink-400">not reached</span>}
                  </span>
                </div>

                <div className="mt-1 text-xs leading-relaxed">
                  {live ? (
                    <MetricReadout metric={live.metrics} info={live.info} passed={live.passed} />
                  ) : gateProse[step.id] ? (
                    <span className="text-ink-500 dark:text-ink-400">{gateProse[step.id]}</span>
                  ) : state === "blocked" && result ? (
                    <span className="text-rose-700 dark:text-rose-300">
                      {blockReason(result.stderr ?? "")}
                    </span>
                  ) : null}
                </div>
              </div>
            </li>
          );
        })}
      </ol>
    </div>
  );
}
