import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import {
  getExample,
  getHealth,
  getModes,
  getMotifs,
  getRealExample,
  getToken,
  setToken,
  streamDesign,
} from "./api";
import type {
  DesignConfig,
  DesignResult,
  Health,
  Mode,
  MotifGroups,
  ProgressEvent,
} from "./types";
import { Header } from "./components/Header";
import { CargoBuilder } from "./components/CargoBuilder";
import { PipelinePanel } from "./components/PipelinePanel";
import { ResultPanel } from "./components/ResultPanel";
import { Badge, Callout, Disclosure, Field, Section, cx } from "./components/Primitives";
import { useTheme } from "./useTheme";

/** PRIDICT2 needs this much sequence either side of the insertion site. */
const MIN_FLANK = 100;

interface FormState {
  gene: string;
  transcript: string;
  uniprot_acc: string;
  rtt_aa: string;
  genome: string;
  site_s: string;
  pegrna_strand: "Fw" | "Rv";
  arm_side: "5prime" | "3prime";
  pbs_lengths: string;
  ngrna_distance_min: string;
  ngrna_distance_max: string;
  ngrna_max_offtargets: string;
  max_run: string;
  spacer: string;
  het_snps: string;
  het_edit: boolean;
}

const EMPTY_FORM: FormState = {
  gene: "",
  transcript: "",
  uniprot_acc: "",
  rtt_aa: "",
  genome: "",
  site_s: "",
  pegrna_strand: "Rv",
  arm_side: "5prime",
  pbs_lengths: "13,15,17",
  ngrna_distance_min: "40",
  ngrna_distance_max: "90",
  ngrna_max_offtargets: "0",
  max_run: "8",
  spacer: "",
  het_snps: "",
  het_edit: false,
};

export default function App() {
  const [theme, toggleTheme] = useTheme();
  const [health, setHealth] = useState<Health | null>(null);
  const [backendDown, setBackendDown] = useState(false);
  const [modes, setModes] = useState<Mode[]>(["fake"]);
  const [mode, setMode] = useState<Mode>("fake");
  const [motifs, setMotifs] = useState<MotifGroups>({});

  const [form, setForm] = useState<FormState>(EMPTY_FORM);
  const [cargoMode, setCargoMode] = useState<"builder" | "raw">("builder");
  const [cargoParts, setCargoParts] = useState<string[]>([]);
  const [rawCargo, setRawCargo] = useState("");

  /** Keys from a loaded example that the form does not surface, preserved on submit. */
  const [passthrough, setPassthrough] = useState<DesignConfig>({});

  /**
   * Which fields the user has actually interacted with, plus whether they have
   * tried to submit. An empty form is incomplete, not wrong -- showing it in red
   * before anyone has typed is just noise, so errors stay hidden until a field is
   * touched or the design is submitted.
   */
  const [touched, setTouched] = useState<Set<string>>(new Set());
  const [submitAttempted, setSubmitAttempted] = useState(false);
  const touch = useCallback(
    (name: string) => setTouched((t) => (t.has(name) ? t : new Set(t).add(name))),
    [],
  );

  const [token, setTokenState] = useState<string>(() => getToken());
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState<Record<string, ProgressEvent>>({});
  const [result, setResult] = useState<DesignResult | null>(null);
  const [submitError, setSubmitError] = useState<string | null>(null);
  const abortRef = useRef<AbortController | null>(null);
  const resultRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => {
    getHealth()
      .then((h) => {
        setHealth(h);
        setBackendDown(false);
      })
      .catch(() => setBackendDown(true));
    getModes()
      .then((m) => {
        setModes(m.modes);
        setMode(m.default);
      })
      .catch(() => undefined);
    getMotifs()
      .then((m) => setMotifs(m.groups))
      .catch(() => undefined);
  }, []);

  const set = useCallback(
    <K extends keyof FormState>(key: K, value: FormState[K]) =>
      setForm((f) => ({ ...f, [key]: value })),
    [],
  );

  const cargo = useMemo(() => {
    if (cargoMode === "raw") return rawCargo;
    const byName = new Map<string, string>();
    Object.values(motifs).forEach((list) => list.forEach((m) => byName.set(m.name, m.dna)));
    return cargoParts.map((n) => byName.get(n) ?? "").join("");
  }, [cargoMode, rawCargo, cargoParts, motifs]);

  function applyConfig(cfg: DesignConfig) {
    setForm({
      gene: cfg.gene ?? "",
      transcript: cfg.transcript ?? "",
      uniprot_acc: cfg.uniprot_acc ?? "",
      rtt_aa: cfg.rtt_aa ?? "",
      genome: cfg.genome ?? "",
      site_s: cfg.site_s != null ? String(cfg.site_s) : "",
      pegrna_strand: cfg.pegrna_strand ?? "Rv",
      arm_side: cfg.arm_side ?? "5prime",
      pbs_lengths: (cfg.pbs_lengths ?? []).join(",") || "13,15,17",
      ngrna_distance_min: cfg.ngrna_distance_min != null ? String(cfg.ngrna_distance_min) : "40",
      ngrna_distance_max: cfg.ngrna_distance_max != null ? String(cfg.ngrna_distance_max) : "90",
      ngrna_max_offtargets:
        cfg.ngrna_max_offtargets != null ? String(cfg.ngrna_max_offtargets) : "0",
      max_run: cfg.max_run != null ? String(cfg.max_run) : "8",
      spacer: cfg.spacer ?? "",
      het_snps: cfg.het_snps ? JSON.stringify(cfg.het_snps) : "",
      het_edit: Boolean(cfg.het_edit),
    });
    // A loaded config carries raw cargo DNA, not motif names.
    setCargoMode("raw");
    setRawCargo(cfg.cargo ?? "");
    setPassthrough(cfg);
    // A freshly loaded example is valid by construction; clear any error state
    // left over from whatever was in the form before.
    setTouched(new Set());
    setSubmitAttempted(false);
    setSubmitError(null);
  }

  const genomeLength = form.genome.length;
  const siteS = Number(form.site_s || 0);
  const leftFlank = siteS;
  const rightFlank = genomeLength - siteS;

  const errors = useMemo(() => {
    const e: Partial<Record<keyof FormState | "cargo", string>> = {};
    if (!form.genome.trim()) e.genome = "A target sequence window is required.";
    else if (/[^ACGTNacgtn]/.test(form.genome.trim()))
      e.genome = "The window should contain only A, C, G, T or N.";
    if (!form.site_s.trim()) e.site_s = "Tell the pipeline where the tag goes.";
    else if (siteS < 0 || siteS > genomeLength)
      e.site_s = `Must be between 0 and ${genomeLength}, the length of the window.`;
    if (!cargo) e.cargo = "Add the cargo you want to insert.";
    if (form.spacer && form.spacer.trim().length !== 20)
      e.spacer = `A spacer is 20 nt (this is ${form.spacer.trim().length}).`;
    const tolerance = Number(form.ngrna_max_offtargets);
    if (!Number.isInteger(tolerance) || tolerance < 0 || tolerance > 10)
      e.ngrna_max_offtargets = "A whole number from 0 to 10.";
    if (form.het_snps.trim()) {
      try {
        const parsed = JSON.parse(form.het_snps);
        if (!Array.isArray(parsed)) e.het_snps = "This must be a JSON array.";
      } catch (err) {
        e.het_snps = `Not valid JSON — ${(err as Error).message}`;
      }
    }
    return e;
  }, [form, cargo, siteS, genomeLength]);

  const hasErrors = Object.keys(errors).length > 0;

  /** The error to display for a field, or undefined while it is still untouched. */
  const shownError = useCallback(
    (name: keyof FormState | "cargo") =>
      submitAttempted || touched.has(name) ? errors[name] : undefined,
    [errors, touched, submitAttempted],
  );

  function buildConfig(): DesignConfig {
    const cfg: DesignConfig = { ...passthrough };
    const assign = (k: string, v: string) => {
      if (v.trim()) cfg[k] = v.trim();
      else delete cfg[k];
    };
    assign("gene", form.gene);
    assign("transcript", form.transcript);
    assign("uniprot_acc", form.uniprot_acc);
    assign("rtt_aa", form.rtt_aa);
    cfg.genome = form.genome.trim().toUpperCase();
    cfg.cargo = cargo;
    cfg.site_s = siteS;
    cfg.pegrna_strand = form.pegrna_strand;
    cfg.arm_side = form.arm_side;
    cfg.max_run = Number(form.max_run) || 8;
    const pbs = form.pbs_lengths
      .split(",")
      .map((s) => Number(s.trim()))
      .filter((n) => Number.isFinite(n) && n > 0);
    if (pbs.length) cfg.pbs_lengths = pbs;
    cfg.ngrna_distance_min = Number(form.ngrna_distance_min) || 40;
    cfg.ngrna_distance_max = Number(form.ngrna_distance_max) || 90;
    cfg.ngrna_max_offtargets = Number(form.ngrna_max_offtargets) || 0;

    // Spacer is an override: blank means PRIDICT2 designs it, and must clear any
    // spacer inherited from a loaded example.
    if (form.spacer.trim()) cfg.spacer = form.spacer.trim().toUpperCase();
    else delete cfg.spacer;

    if (form.het_snps.trim()) cfg.het_snps = JSON.parse(form.het_snps);
    else delete cfg.het_snps;
    if (form.het_edit) cfg.het_edit = true;
    else delete cfg.het_edit;

    // Coding up to the insertion site, non-coding after it.
    cfg.cds_mask = Array.from({ length: cfg.genome.length }, (_, i) => i < siteS);
    return cfg;
  }

  async function run() {
    setSubmitAttempted(true);
    if (hasErrors) {
      setSubmitError("Fix the highlighted fields first.");
      return;
    }
    setSubmitError(null);
    setResult(null);
    setProgress({});
    setRunning(true);
    abortRef.current = new AbortController();
    try {
      const res = await streamDesign(
        buildConfig(),
        mode,
        (ev) => setProgress((p) => ({ ...p, [ev.gate]: ev })),
        abortRef.current.signal,
      );
      setResult(res);
      requestAnimationFrame(() =>
        resultRef.current?.scrollIntoView({ behavior: "smooth", block: "start" }),
      );
    } catch (err) {
      if ((err as Error).name !== "AbortError") setSubmitError((err as Error).message);
    } finally {
      setRunning(false);
      abortRef.current = null;
    }
  }

  return (
    <div className="min-h-screen">
      <Header health={health} theme={theme} onToggleTheme={toggleTheme} />

      <main className="mx-auto max-w-7xl px-4 py-6 sm:px-6">
        <div className="grid gap-6 lg:grid-cols-[minmax(0,1fr)_380px]">
          {/* ---- form column ---- */}
          <div className="space-y-5">
            {backendDown && (
              <Callout tone="error" title="The design server is not responding">
                The page has loaded, but it cannot reach the machine that runs the pipeline.
                It may be offline, or its address may have changed. Designs will not run until
                it is back.
              </Callout>
            )}

            <Callout tone="info">
              New here? Load the worked <strong>CTCF</strong> example below to see every field
              filled in correctly, then change one thing at a time.
              <span className="mt-2 flex flex-wrap gap-2">
                <button
                  type="button"
                  className="btn btn-secondary px-2.5 py-1 text-xs"
                  onClick={() => getExample().then(applyConfig).catch(() => undefined)}
                >
                  Load CTCF example
                </button>
                <button
                  type="button"
                  className="btn btn-secondary px-2.5 py-1 text-xs"
                  onClick={() =>
                    getRealExample()
                      .then((cfg) => {
                        applyConfig(cfg);
                        if (modes.includes("api")) setMode("api");
                      })
                      .catch(() => undefined)
                  }
                >
                  Load example for real models
                </button>
              </span>
            </Callout>

            <Section
              step={1}
              title="What are you tagging?"
              description="The gene and the stretch of sequence around the spot where the tag goes."
            >
              <div className="space-y-4">
                <div className="grid gap-4 sm:grid-cols-3">
                  <Field label="Gene" htmlFor="gene" hint="For the report only.">
                    <input
                      id="gene"
                      className="input"
                      value={form.gene}
                      onChange={(e) => set("gene", e.target.value)}
                      placeholder="CTCF"
                    />
                  </Field>
                  <Field
                    label="Transcript"
                    htmlFor="transcript"
                    help="The Ensembl transcript ID. If your window is too short, the pipeline uses this to fetch more sequence from Ensembl."
                  >
                    <input
                      id="transcript"
                      className="input"
                      value={form.transcript}
                      onChange={(e) => set("transcript", e.target.value)}
                      placeholder="ENST00000264010"
                    />
                  </Field>
                  <Field
                    label="UniProt"
                    htmlFor="uniprot_acc"
                    help="Used to pull the AlphaFold structure, so the pipeline can check the tag goes somewhere the protein tolerates."
                  >
                    <input
                      id="uniprot_acc"
                      className="input"
                      value={form.uniprot_acc}
                      onChange={(e) => set("uniprot_acc", e.target.value)}
                      placeholder="P49711"
                    />
                  </Field>
                </div>

                <Field
                  label="Target sequence window"
                  htmlFor="genome"
                  error={shownError("genome")}
                  help="Genomic sequence (5′→3′) spanning your insertion point, with at least 100 bases either side."
                  hint={
                    genomeLength > 0
                      ? `${genomeLength} bases pasted.`
                      : "Paste the genomic sequence around your insertion point."
                  }
                >
                  <textarea
                    id="genome"
                    onBlur={() => touch("genome")}
                    rows={4}
                    className={cx("input seq text-xs", shownError("genome") && "input-invalid")}
                    value={form.genome}
                    onChange={(e) => set("genome", e.target.value.replace(/\s/g, ""))}
                    placeholder="GGTGCAATTGAGAACATTATAGTT…"
                  />
                </Field>

                <div className="grid gap-4 sm:grid-cols-2">
                  <Field
                    label="Insertion point"
                    htmlFor="site_s"
                    error={shownError("site_s")}
                    help="Position in the window where the tag is inserted, counting from 0. For a C-terminal tag this is the stop codon."
                  >
                    <input
                      id="site_s"
                    onBlur={() => touch("site_s")}
                      type="number"
                      className={cx("input", shownError("site_s") && "input-invalid")}
                      value={form.site_s}
                      onChange={(e) => set("site_s", e.target.value)}
                      placeholder="147"
                    />
                  </Field>
                  <Field
                    label="Anchor amino acid (optional)"
                    htmlFor="rtt_aa"
                    help="Does not change the design: the recoded arm and the RTT are built from the genome and your cargo, so any residue, or none, gives the same pegRNA. Kept so older configs still load."
                  >
                    <input
                      id="rtt_aa"
                      className="input"
                      value={form.rtt_aa}
                      onChange={(e) => set("rtt_aa", e.target.value.toUpperCase().slice(0, 1))}
                      placeholder="not needed"
                    />
                  </Field>
                </div>

                {genomeLength > 0 && form.site_s !== "" && (
                  <div className="flex flex-wrap items-center gap-2">
                    <Badge tone={leftFlank >= MIN_FLANK ? "ok" : "warn"}>
                      {leftFlank} bases before
                    </Badge>
                    <Badge tone={rightFlank >= MIN_FLANK ? "ok" : "warn"}>
                      {rightFlank} bases after
                    </Badge>
                    {(leftFlank < MIN_FLANK || rightFlank < MIN_FLANK) && (
                      <span className="text-xs text-amber-700 dark:text-amber-400">
                        {form.transcript
                          ? `Under ${MIN_FLANK} — the pipeline will fetch more sequence from Ensembl for ${form.transcript}.`
                          : `Under ${MIN_FLANK}. Paste a wider window, or add a transcript so more can be fetched.`}
                      </span>
                    )}
                  </div>
                )}
              </div>
            </Section>

            <Section
              step={2}
              title="What are you inserting?"
              description="The cargo knocked in at the insertion point — a tag, linker, marker, or any combination."
            >
              {shownError("cargo") && (
                <p role="alert" className="mb-3 text-xs font-medium text-rose-600">
                  {shownError("cargo")}
                </p>
              )}
              <CargoBuilder
                motifs={motifs}
                parts={cargoParts}
                onPartsChange={setCargoParts}
                rawCargo={rawCargo}
                onRawCargoChange={setRawCargo}
                mode={cargoMode}
                onModeChange={setCargoMode}
              />
            </Section>

            <Section
              step={3}
              title="How should it be designed?"
              description="Sensible defaults are filled in. Change these only if you have a reason."
            >
              <div className="space-y-4">
                <div className="grid gap-4 sm:grid-cols-2">
                  <Field
                    label="pegRNA strand"
                    htmlFor="pegrna_strand"
                    help="Which strand the edit is written on. This names the flap strand, not where the guide reads."
                  >
                    <select
                      id="pegrna_strand"
                      className="input"
                      value={form.pegrna_strand}
                      onChange={(e) => set("pegrna_strand", e.target.value as "Fw" | "Rv")}
                    >
                      <option value="Fw">Forward</option>
                      <option value="Rv">Reverse</option>
                    </select>
                  </Field>
                  <Field
                    label="Recoded arm side"
                    htmlFor="arm_side"
                    help="Which side of the insertion point gets rewritten with synonymous codons to break homology."
                  >
                    <select
                      id="arm_side"
                      className="input"
                      value={form.arm_side}
                      onChange={(e) => set("arm_side", e.target.value as "5prime" | "3prime")}
                    >
                      <option value="5prime">5′ side</option>
                      <option value="3prime">3′ side</option>
                    </select>
                  </Field>
                </div>

                <Disclosure
                  title="Fine-tuning"
                  description="PBS lengths, nick distance, homopolymer limit, off-target tolerance"
                >
                  <div className="grid gap-4 sm:grid-cols-2">
                    <Field
                      label="PBS lengths to try"
                      htmlFor="pbs_lengths"
                      help="Primer-binding site lengths, in nucleotides. The pipeline tries each and keeps the first that passes everything."
                    >
                      <input
                        id="pbs_lengths"
                        className="input"
                        value={form.pbs_lengths}
                        onChange={(e) => set("pbs_lengths", e.target.value)}
                        placeholder="13,15,17"
                      />
                    </Field>
                    <Field
                      label="Longest allowed repeat"
                      htmlFor="max_run"
                      help="Maximum run of one base. Long homopolymers are hard to synthesise and encourage misrepair."
                    >
                      <input
                        id="max_run"
                        type="number"
                        className="input"
                        value={form.max_run}
                        onChange={(e) => set("max_run", e.target.value)}
                      />
                    </Field>
                    <Field label="Nick distance, minimum" htmlFor="ngrna_distance_min">
                      <input
                        id="ngrna_distance_min"
                        type="number"
                        className="input"
                        value={form.ngrna_distance_min}
                        onChange={(e) => set("ngrna_distance_min", e.target.value)}
                      />
                    </Field>
                    <Field
                      label="Nick distance, maximum"
                      htmlFor="ngrna_distance_max"
                      help="How far the second nick may sit from the insertion point, in bases."
                    >
                      <input
                        id="ngrna_distance_max"
                        type="number"
                        className="input"
                        value={form.ngrna_distance_max}
                        onChange={(e) => set("ngrna_distance_max", e.target.value)}
                      />
                    </Field>
                    <Field
                      label="Off-target hits allowed, second nick"
                      htmlFor="ngrna_max_offtargets"
                      error={shownError("ngrna_max_offtargets")}
                      help="How many near-matches elsewhere in the genome the second nick's guide may have. The second nick cuts one strand only, so a near-match costs less than for a cutting guide. 0 is strict and recommended; raise it only if every good guide is being rejected. Whatever the chosen guide has is shown on the result."
                    >
                      <input
                        id="ngrna_max_offtargets"
                        type="number"
                        min={0}
                        max={10}
                        className={cx(
                          "input",
                          shownError("ngrna_max_offtargets") && "input-invalid",
                        )}
                        value={form.ngrna_max_offtargets}
                        onBlur={() => touch("ngrna_max_offtargets")}
                        onChange={(e) => set("ngrna_max_offtargets", e.target.value)}
                      />
                    </Field>
                  </div>
                </Disclosure>

                <Disclosure title="Use my own guide" description="optional override">
                  <Field
                    label="Spacer"
                    htmlFor="spacer"
                    error={shownError("spacer")}
                    hint="Leave blank and PRIDICT2 designs the guide for you. This is the recommended path."
                  >
                    <input
                      id="spacer"
                    onBlur={() => touch("spacer")}
                      className={cx("input seq", shownError("spacer") && "input-invalid")}
                      value={form.spacer}
                      onChange={(e) => set("spacer", e.target.value.toUpperCase())}
                      placeholder="20-nt spacer, e.g. TCAGCATGATGGACCGGTGA"
                    />
                  </Field>
                </Disclosure>

                <Disclosure
                  title="Allele-specific design"
                  description="for heterozygous sites"
                >
                  <Field
                    label="Phased heterozygous SNPs"
                    htmlFor="het_snps"
                    error={shownError("het_snps")}
                    help="A JSON array locating each SNP by the guide element it falls in. Discrimination is strongest in the PAM, then the seed region."
                    hint='e.g. [{"element":"pegrna_pam","offset":1,"haplotype":1,"phased":true}]'
                  >
                    <textarea
                      id="het_snps"
                    onBlur={() => touch("het_snps")}
                      rows={3}
                      className={cx("input seq text-xs", shownError("het_snps") && "input-invalid")}
                      value={form.het_snps}
                      onChange={(e) => set("het_snps", e.target.value)}
                    />
                  </Field>
                  <label className="flex items-center gap-2 text-sm">
                    <input
                      type="checkbox"
                      checked={form.het_edit}
                      onChange={(e) => set("het_edit", e.target.checked)}
                      className="h-4 w-4 rounded border-ink-300 text-brand-600 focus:ring-brand-500"
                    />
                    Also design an allele-specific pair
                  </label>
                  <p className="field-hint">
                    Allele specificity is limited by the natural variation at your locus — it is
                    never manufactured, so a pair appears only when a discriminating SNP exists.
                  </p>
                </Disclosure>
              </div>
            </Section>
          </div>

          {/* ---- pipeline column ---- */}
          <div className="space-y-5 lg:sticky lg:top-6 lg:self-start">
            <PipelinePanel running={running} progress={progress} result={result} />
          </div>
        </div>

        {/*
          Sticky action bar. It spans the viewport (not just the grid) and carries
          an upward shadow, so it reads as a toolbar pinned above the page rather
          than a panel floating on top of the form.
        */}
        <div
          className="sticky bottom-0 z-20 mt-6 -mx-4 border-t border-ink-200 bg-white/95
                     shadow-[0_-4px_16px_-6px_rgba(0,0,0,0.18)] backdrop-blur
                     sm:-mx-6 dark:border-ink-800 dark:bg-ink-900/95"
        >
          <div className="flex flex-wrap items-end gap-x-4 gap-y-3 px-4 py-3 sm:px-6">
            <div className="w-44">
              <label className="field-label" htmlFor="mode">
                Mode
              </label>
              <select
                id="mode"
                className="input"
                value={mode}
                onChange={(e) => setMode(e.target.value as Mode)}
              >
                {modes.map((m) => (
                  <option key={m} value={m}>
                    {m === "fake" ? "Dry run (fast)" : "Real models"}
                  </option>
                ))}
              </select>
            </div>
            {/*
              The token is only asked for when the server says api mode needs
              one, and never for fake mode -- a visitor should be able to try the
              tool without credentials.
            */}
            {mode === "api" && health?.api_mode_token_required && (
              <div className="w-56">
                <label className="field-label" htmlFor="access-token">
                  Access token
                </label>
                <input
                  id="access-token"
                  type="password"
                  className="input"
                  value={token}
                  onChange={(e) => {
                    setTokenState(e.target.value);
                    setToken(e.target.value);
                  }}
                  placeholder="required for real models"
                />
              </div>
            )}
            <p className="max-w-sm flex-1 pb-2 text-xs leading-relaxed text-ink-500 dark:text-ink-400">
              {mode === "fake"
                ? "Checks your inputs flow through the design logic in about a second. Every verdict is a placeholder."
                : health?.api_mode_token_required && !token
                  ? "Running the real models needs an access token. Dry run is open to everyone."
                  : "Runs the real models — this can take 15 minutes or more. Progress appears as each check completes."}
            </p>
            {running && (
              <button
                type="button"
                className="btn btn-secondary"
                onClick={() => abortRef.current?.abort()}
              >
                Cancel
              </button>
            )}
            <button
              type="button"
              className="btn btn-primary px-6 py-2.5"
              onClick={run}
              disabled={running}
            >
              {running && (
                <span
                  aria-hidden
                  className="h-3.5 w-3.5 animate-spin rounded-full border-2 border-white/40 border-t-white"
                />
              )}
              {running ? "Designing…" : "Design pegRNA"}
            </button>
          </div>
          {submitError && (
            <p
              role="alert"
              className="border-t border-rose-200 bg-rose-50 px-4 py-2 text-sm font-medium
                         text-rose-700 sm:px-6 dark:border-rose-900 dark:bg-rose-950 dark:text-rose-300"
            >
              {submitError}
            </p>
          )}
        </div>

        <div ref={resultRef} className="mt-6">
          {result && <ResultPanel result={result} />}
        </div>
      </main>

      <footer className="mx-auto max-w-7xl px-4 py-8 text-xs text-ink-500 sm:px-6 dark:text-ink-400">
        Recoded TJ-PE pegRNA designer · designs are engineering previews — validate at the bench
        before ordering.
      </footer>
    </div>
  );
}
