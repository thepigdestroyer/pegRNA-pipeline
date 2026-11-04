import { useEffect, useState } from "react";

export type Theme = "light" | "dark";

const KEY = "pegrna.theme";

/** Browser chrome colours (address bar, title bar), matching the page background. */
const THEME_COLOR: Record<Theme, string> = { light: "#ffffff", dark: "#0b0f14" };

/** Theme preference, persisted, defaulting to the OS setting. */
export function useTheme(): [Theme, () => void] {
  const [theme, setTheme] = useState<Theme>(() => {
    const saved = localStorage.getItem(KEY);
    if (saved === "light" || saved === "dark") return saved;
    return window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  });

  useEffect(() => {
    document.documentElement.classList.toggle("dark", theme === "dark");
    // index.html ships one theme-color tag per OS setting; once the toggle decides, they must agree.
    document.querySelectorAll('meta[name="theme-color"]').forEach((meta) => {
      meta.removeAttribute("media");
      meta.setAttribute("content", THEME_COLOR[theme]);
    });
    localStorage.setItem(KEY, theme);
  }, [theme]);

  return [theme, () => setTheme((t) => (t === "dark" ? "light" : "dark"))];
}
