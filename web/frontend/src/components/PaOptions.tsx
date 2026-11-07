import { PA_LIMITS, type PaFieldName, type PaFields } from "../pa";
import { Callout, Disclosure, Field, cx } from "./Primitives";

interface Props {
  form: PaFields;
  set: (key: keyof PaFields, value: string) => void;
  error: (name: PaFieldName) => string | undefined;
  touch: (name: string) => void;
}

/** The settings of a Prime Assembly design. The guides are chosen by the pipeline, so there is no guide input. */
export function PaOptions({ form, set, error, touch }: Props) {
  const input = (name: keyof PaFields, extra?: { min?: number; max?: number; type?: string }) => ({
    id: name,
    type: extra?.type ?? "number",
    min: extra?.min,
    max: extra?.max,
    className: cx("input", error(name) && "input-invalid"),
    value: form[name],
    onBlur: () => touch(name),
    onChange: (e: { target: { value: string } }) => set(name, e.target.value),
  });

  return (
    <div className="space-y-4">
      <Callout tone="neutral">
        Prime Assembly chooses the guides for you. It screens every pair of NGG guides that face each other across
        your insertion point, tries the pairs best first, and keeps the first that passes every check.
      </Callout>

      <div className="grid gap-4 sm:grid-cols-2">
        <Field
          label="Flap length (nt)"
          htmlFor="flap_length"
          error={error("flap_length")}
          help="How many nucleotides each pegRNA writes, and so how much of each end of the donor must match it. The paper used 35 and found 30 to 50 robust (13 was very weak). The insert must be at least twice this."
        >
          <input {...input("flap_length", { min: PA_LIMITS.flapMin, max: PA_LIMITS.flapMax })} />
        </Field>
        <Field
          label="Off-target hits allowed, either guide"
          htmlFor="max_offtargets"
          error={error("max_offtargets")}
          help="Both guides are screened against the whole genome. Each cuts one strand only, so a near-match costs less than for a cutting guide. 0 is strict and recommended; raise it only if no pair survives. In the CTCF example, 0 leaves one usable guide pair and 1 leaves four. Whatever each chosen guide has is shown on the result, and an exact duplicate elsewhere is never allowed."
        >
          <input {...input("max_offtargets", { min: 0, max: PA_LIMITS.maxOfftargets })} />
        </Field>
      </div>

      <Disclosure
        title="Fine-tuning"
        description="PBS lengths, distance between the nicks, repeat limit, pairs to try"
      >
        <div className="grid gap-4 sm:grid-cols-2">
          <Field
            label="PBS lengths to try"
            htmlFor="pbs_lengths"
            error={error("pbs_lengths")}
            help={`Primer-binding site lengths to try for each pegRNA, in nucleotides, up to ${PA_LIMITS.pbsCount} values. Every one is folded and docked.`}
          >
            <input {...input("pbs_lengths", { type: "text" })} placeholder="13,15,17" />
          </Field>
          <Field
            label="Longest allowed repeat"
            htmlFor="max_run"
            error={error("max_run")}
            help="The longest run of identical bases the donor's recoded arm may keep from the genome. A long identical stretch lets a flap re-anneal where it came from instead of on the donor."
          >
            <input {...input("max_run", { min: 1 })} />
          </Field>
          <Field
            label="Distance between the nicks, minimum (bp)"
            htmlFor="min_span"
            error={error("min_span")}
            help={`The donor restores everything between the two nicks, so a wider span means more DNA to recode. ${PA_LIMITS.spanFloor} is a hard floor.`}
          >
            <input {...input("min_span", { min: PA_LIMITS.spanFloor, max: PA_LIMITS.spanCap })} />
          </Field>
          <Field
            label="Distance between the nicks, maximum (bp)"
            htmlFor="max_span"
            error={error("max_span")}
          >
            <input {...input("max_span", { min: PA_LIMITS.spanFloor, max: PA_LIMITS.spanCap })} />
          </Field>
          <Field
            label="Guide pairs to try"
            htmlFor="max_pairs"
            error={error("max_pairs")}
            help={`The walk tries guide pairs best first and stops at the first that passes everything. This caps how many it tries, at most ${PA_LIMITS.maxPairs}; each pair can need a Cas9-binding (docking) call.`}
          >
            <input {...input("max_pairs", { min: 1, max: PA_LIMITS.maxPairs })} />
          </Field>
        </div>
      </Disclosure>
    </div>
  );
}
