import type { Health } from "../types";
import { cx } from "./Primitives";

function Chip({
  tone,
  hidden,
  children,
}: {
  tone: "ok" | "warn";
  hidden?: boolean;
  children: React.ReactNode;
}) {
  return (
    <span
      className={cx(
        "inline-flex items-center gap-1.5 text-xs text-ink-600 dark:text-ink-300",
        hidden && "invisible",
      )}
    >
      <span
        aria-hidden
        className={cx(
          "h-1.5 w-1.5 rounded-full",
          tone === "ok" ? "bg-emerald-500" : "bg-amber-500",
        )}
      />
      {children}
    </span>
  );
}

/** What the chips say when everything is up. Used to hold their space until the backend answers. */
const EXPECTED_CHIPS = ["service ready", "real models installed", "docking live"];

export function Header({
  health,
  theme,
  onToggleTheme,
}: {
  health: Health | null;
  theme: "light" | "dark";
  onToggleTheme: () => void;
}) {
  return (
    <header className="border-b border-ink-200 bg-white dark:border-ink-800 dark:bg-ink-900">
      <div className="mx-auto flex max-w-7xl flex-wrap items-center gap-x-6 gap-y-3 px-4 py-4 sm:px-6">
        <div className="order-1 min-w-0 flex-1">
          <h1 className="text-lg font-semibold tracking-tight text-ink-900 dark:text-ink-50">
            pegRNA designer
          </h1>
          <p className="mt-0.5 text-sm text-ink-600 dark:text-ink-400">
            Knock a tag into the C-terminus of a human gene, by prime editing.
          </p>
        </div>

        {/*
          The status chips sit under the title on a phone and beside it from lg up. Beside it on a
          narrow screen they crushed the subtitle into a column one word wide. Their space is held
          (hidden, same words) until the backend answers, so nothing moves when they appear.
        */}
        <div
          className="order-3 flex w-full flex-wrap items-center gap-x-4 gap-y-1.5 lg:order-2 lg:w-auto"
          aria-live="polite"
          aria-label="Service status"
        >
          {health ? (
            <>
              <Chip tone={health.repo_ok ? "ok" : "warn"}>
                {health.repo_ok ? "service ready" : "pipeline not found"}
              </Chip>
              <Chip tone={health.api_python_ok ? "ok" : "warn"}>
                {health.api_python_ok ? "real models installed" : "real models unavailable"}
              </Chip>
              <Chip tone={health.nvidia_key ? "ok" : "warn"}>
                {health.nvidia_key ? "docking live" : "docking stand-in"}
              </Chip>
            </>
          ) : (
            EXPECTED_CHIPS.map((text) => (
              <Chip key={text} tone="ok" hidden>
                {text}
              </Chip>
            ))
          )}
        </div>

        <button
          type="button"
          onClick={onToggleTheme}
          className="btn btn-secondary order-2 px-2.5 py-1.5 lg:order-3"
          aria-label={`Switch to ${theme === "dark" ? "light" : "dark"} theme`}
        >
          {theme === "dark" ? "☀" : "☾"}
        </button>
      </div>
    </header>
  );
}
