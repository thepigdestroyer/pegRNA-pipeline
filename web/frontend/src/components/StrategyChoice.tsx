import type { Advice } from "../pa";
import type { Strategy } from "../types";
import { Callout, cx } from "./Primitives";

interface Option {
  value: Strategy;
  title: string;
  fits: string;
  body: string;
}

const OPTIONS: Option[] = [
  {
    value: "tjpe",
    title: "Inside the pegRNA",
    fits: "tags and short inserts",
    body: "The insert is carried in the pegRNA's own extension (tj-PE), so there is one RNA to order. Best for a tag or a short insert: the pegRNA grows with the insert, and past a few dozen nucleotides the Cas9-binding check can no longer judge it.",
  },
  {
    value: "pa",
    title: "As a separate donor DNA",
    fits: "cassettes and genes",
    body: "Prime Assembly. Two pegRNAs write short flaps at two nicks, and a donor DNA you make by PCR, whose ends match the flaps, supplies the rest. Two short pegRNAs plus a donor to order, and the insert can be a whole cassette (at least twice the flap length, 70 nt by default). This version restores the DNA between the nicks, which has not been tested at the bench.",
  },
];

interface Props {
  value: Strategy;
  onChange: (value: Strategy) => void;
  /** Locked while a design runs: the pipeline being drawn must stay the one that is running. */
  disabled: boolean;
  /** False when the server does not offer Prime Assembly (an older deployment). */
  paAvailable: boolean;
  advice: Advice | null;
}

export function StrategyChoice({ value, onChange, disabled, paAvailable, advice }: Props) {
  return (
    <div className="space-y-3">
      <div role="radiogroup" aria-label="How the insert is delivered" className="grid gap-3 sm:grid-cols-2">
        {OPTIONS.map((o) => {
          const selected = value === o.value;
          const unavailable = o.value === "pa" && !paAvailable;
          return (
            <button
              key={o.value}
              type="button"
              role="radio"
              aria-checked={selected}
              disabled={disabled || unavailable}
              onClick={() => onChange(o.value)}
              className={cx(
                "rounded-xl border p-4 text-left transition-colors focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-brand-500",
                selected
                  ? "border-brand-500 bg-brand-50 ring-1 ring-brand-500 dark:bg-brand-700/15"
                  : "border-ink-200 bg-white hover:border-ink-300 dark:border-ink-700 dark:bg-ink-900 dark:hover:border-ink-600",
                (disabled || unavailable) && "cursor-not-allowed opacity-60",
              )}
            >
              <span className="flex items-center gap-2">
                <span
                  aria-hidden
                  className={cx(
                    "flex h-4 w-4 shrink-0 items-center justify-center rounded-full border-2",
                    selected ? "border-brand-600" : "border-ink-300 dark:border-ink-600",
                  )}
                >
                  {selected && <span className="h-2 w-2 rounded-full bg-brand-600" />}
                </span>
                <span className="text-sm font-semibold text-ink-900 dark:text-ink-50">{o.title}</span>
                <span className="ml-auto text-xs text-ink-500 dark:text-ink-400">{o.fits}</span>
              </span>
              <span className="mt-2 block text-xs leading-relaxed text-ink-600 dark:text-ink-400">
                {unavailable ? "This server has not been updated for Prime Assembly yet." : o.body}
              </span>
            </button>
          );
        })}
      </div>

      {advice && (
        <Callout tone={advice.tone}>
          {advice.text}
          {advice.switchTo && advice.action && (
            <span className="mt-2 flex">
              <button
                type="button"
                className="btn btn-secondary px-2.5 py-1 text-xs"
                disabled={disabled || (advice.switchTo === "pa" && !paAvailable)}
                onClick={() => onChange(advice.switchTo as Strategy)}
              >
                {advice.action}
              </button>
            </span>
          )}
        </Callout>
      )}
    </div>
  );
}
