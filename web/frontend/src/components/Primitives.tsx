import { useId, useState, type ReactNode } from "react";

export function cx(...parts: (string | false | null | undefined)[]): string {
  return parts.filter(Boolean).join(" ");
}

/**
 * An inline "what is this?" affordance.
 *
 * Every piece of jargon on this page gets one. It is a real <button> with the
 * explanation rendered as text (not a title attribute) so it works on touch and
 * is reachable by keyboard and screen readers.
 */
export function InfoHint({ label, children }: { label: string; children: ReactNode }) {
  const [open, setOpen] = useState(false);
  const id = useId();
  return (
    <span className="relative inline-block align-middle">
      <button
        type="button"
        aria-expanded={open}
        aria-controls={id}
        aria-label={`What is ${label}?`}
        onClick={() => setOpen((v) => !v)}
        className="ml-1 inline-flex h-4 w-4 items-center justify-center rounded-full border
                   border-ink-300 text-[10px] leading-none font-semibold text-ink-500
                   hover:border-brand-500 hover:text-brand-600 dark:border-ink-600 dark:text-ink-400"
      >
        ?
      </button>
      {open && (
        <span
          id={id}
          role="note"
          className="absolute top-6 left-0 z-30 w-72 rounded-lg border border-ink-200 bg-white
                     p-3 text-xs leading-relaxed font-normal text-ink-700 shadow-lg
                     dark:border-ink-700 dark:bg-ink-800 dark:text-ink-200"
        >
          {children}
        </span>
      )}
    </span>
  );
}

interface FieldProps {
  label: string;
  htmlFor?: string;
  hint?: ReactNode;
  help?: ReactNode;
  error?: string;
  children: ReactNode;
  className?: string;
}

export function Field({ label, htmlFor, hint, help, error, children, className }: FieldProps) {
  return (
    <div className={className}>
      <label className="field-label" htmlFor={htmlFor}>
        {label}
        {help && <InfoHint label={label}>{help}</InfoHint>}
      </label>
      {children}
      {error ? (
        <p role="alert" className="mt-1.5 text-xs font-medium text-rose-600 dark:text-rose-400">
          {error}
        </p>
      ) : (
        hint && <p className="field-hint">{hint}</p>
      )}
    </div>
  );
}

type Tone = "neutral" | "ok" | "warn" | "error" | "info";

const TONES: Record<Tone, string> = {
  neutral: "bg-ink-100 text-ink-700 border-ink-200 dark:bg-ink-800 dark:text-ink-200 dark:border-ink-700",
  ok: "bg-emerald-50 text-emerald-700 border-emerald-200 dark:bg-emerald-950 dark:text-emerald-300 dark:border-emerald-900",
  warn: "bg-amber-50 text-amber-800 border-amber-200 dark:bg-amber-950 dark:text-amber-300 dark:border-amber-900",
  error: "bg-rose-50 text-rose-700 border-rose-200 dark:bg-rose-950 dark:text-rose-300 dark:border-rose-900",
  info: "bg-brand-50 text-brand-700 border-brand-100 dark:bg-brand-700/20 dark:text-brand-300 dark:border-brand-700/40",
};

export function Badge({
  tone = "neutral",
  children,
  className,
}: {
  tone?: Tone;
  children: ReactNode;
  className?: string;
}) {
  return (
    <span
      className={cx(
        "inline-flex items-center gap-1.5 rounded-full border px-2.5 py-0.5 text-xs font-medium",
        TONES[tone],
        className,
      )}
    >
      {children}
    </span>
  );
}

export function Callout({
  tone = "info",
  title,
  children,
}: {
  tone?: Tone;
  title?: ReactNode;
  children: ReactNode;
}) {
  return (
    <div className={cx("rounded-lg border px-3.5 py-3 text-sm", TONES[tone])}>
      {title && <p className="mb-1 font-semibold">{title}</p>}
      <div className="leading-relaxed [&_code]:font-mono [&_code]:text-xs">{children}</div>
    </div>
  );
}

/** A numbered form section, so the page reads as a sequence of decisions. */
export function Section({
  step,
  title,
  description,
  children,
  aside,
}: {
  step: number;
  title: string;
  description?: ReactNode;
  children: ReactNode;
  aside?: ReactNode;
}) {
  return (
    <section className="card p-5 sm:p-6">
      <header className="mb-5 flex items-start gap-3">
        <span
          aria-hidden
          className="mt-0.5 flex h-7 w-7 shrink-0 items-center justify-center rounded-full
                     bg-brand-600 text-sm font-semibold text-white"
        >
          {step}
        </span>
        <div className="min-w-0 flex-1">
          <h2 className="text-base font-semibold text-ink-900 dark:text-ink-50">{title}</h2>
          {description && (
            <p className="mt-1 text-sm leading-relaxed text-ink-600 dark:text-ink-400">
              {description}
            </p>
          )}
        </div>
        {aside}
      </header>
      {children}
    </section>
  );
}

/** Progressive disclosure for the expert-only inputs. */
export function Disclosure({
  title,
  description,
  children,
  defaultOpen = false,
}: {
  title: string;
  description?: string;
  children: ReactNode;
  defaultOpen?: boolean;
}) {
  return (
    <details
      className="group rounded-lg border border-ink-200 bg-ink-50/60 dark:border-ink-800 dark:bg-ink-950/40"
      open={defaultOpen}
    >
      <summary
        className="flex cursor-pointer list-none items-center gap-2 px-3.5 py-2.5 text-sm
                   font-medium text-ink-700 select-none dark:text-ink-200"
      >
        <svg
          aria-hidden
          viewBox="0 0 20 20"
          className="h-4 w-4 shrink-0 text-ink-400 transition-transform group-open:rotate-90"
          fill="currentColor"
        >
          <path d="M7.5 5l5 5-5 5V5z" />
        </svg>
        <span className="flex-1">{title}</span>
        {description && (
          <span className="hidden text-xs font-normal text-ink-500 sm:inline">{description}</span>
        )}
      </summary>
      <div className="space-y-4 border-t border-ink-200 px-3.5 py-4 dark:border-ink-800">
        {children}
      </div>
    </details>
  );
}

export function CopyButton({ value, label = "Copy" }: { value: string; label?: string }) {
  const [done, setDone] = useState(false);
  return (
    <button
      type="button"
      className="btn btn-ghost px-2 py-1 text-xs"
      onClick={() => {
        navigator.clipboard
          .writeText(value)
          .then(() => {
            setDone(true);
            setTimeout(() => setDone(false), 1500);
          })
          .catch(() => {
            /* clipboard can be blocked; the text stays selectable either way */
          });
      }}
    >
      {done ? "Copied" : label}
    </button>
  );
}
