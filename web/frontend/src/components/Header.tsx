import type { Health } from "../types";
import { cx } from "./Primitives";

function Chip({ tone, children }: { tone: "ok" | "warn"; children: React.ReactNode }) {
  return (
    <span className="inline-flex items-center gap-1.5 text-xs text-ink-600 dark:text-ink-300">
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
        <div className="min-w-0 flex-1">
          <h1 className="text-lg font-semibold tracking-tight text-ink-900 dark:text-ink-50">
            pegRNA designer
          </h1>
          <p className="mt-0.5 text-sm text-ink-600 dark:text-ink-400">
            Knock a tag into the C-terminus of a human gene, by prime editing.
          </p>
        </div>

        {health && (
          <div
            className="flex flex-wrap items-center gap-x-4 gap-y-1.5"
            aria-live="polite"
            aria-label="Service status"
          >
            <Chip tone={health.repo_ok ? "ok" : "warn"}>
              {health.repo_ok ? "service ready" : "pipeline not found"}
            </Chip>
            <Chip tone={health.api_python_ok ? "ok" : "warn"}>
              {health.api_python_ok ? "real models installed" : "real models unavailable"}
            </Chip>
            <Chip tone={health.nvidia_key ? "ok" : "warn"}>
              {health.nvidia_key ? "docking live" : "docking stand-in"}
            </Chip>
          </div>
        )}

        <button
          type="button"
          onClick={onToggleTheme}
          className="btn btn-secondary px-2.5 py-1.5"
          aria-label={`Switch to ${theme === "dark" ? "light" : "dark"} theme`}
        >
          {theme === "dark" ? "☀" : "☾"}
        </button>
      </div>
    </header>
  );
}
