#!/usr/bin/env python3
"""On-CX3 smoke validation for every REAL gate model.

The real adapter bodies (Boltz-2, SpliceAI, AlphaFold-importance, RiNALMo,
RNAformer) were written against documented APIs but could NOT be executed on the
development box (no Ampere GPU / no weights). This harness is the gate that
proves each one actually runs and returns sane output **on CX3 with staged
weights**, BEFORE any production design run trusts it.

Run after staging:

    python scripts/hpc/smoke_models.py --data $PEGRNA_DATA

Exit code 0 only if every model that is expected to be present runs and returns
a value in range. Any failure prints the model + error and exits non-zero — so
``run_recoded_pipeline.pbs`` aborts before designing against a broken gate.
"""

from __future__ import annotations

import argparse
import sys
import traceback

_SEQ = "GGGAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC"
_RTT_RNA = "ACGUACGUACGUACGUACGU"


def _check(name, fn):
    try:
        ok, detail = fn()
        status = "PASS" if ok else "FAIL"
        print(f"[{status}] {name}: {detail}")
        return ok
    except Exception as exc:  # noqa: BLE001
        print(f"[ERROR] {name}: {type(exc).__name__}: {exc}")
        traceback.print_exc()
        return False


def smoke_splice():
    from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer

    d = SpliceAIScorer().delta("A" * 11000, "A" * 5500 + "T" + "A" * 5499)
    return (0.0 <= d <= 1.0, f"delta={d:.4f}")


def smoke_importance(data):
    # Calmodulin P0DP23: small, single-fragment, verified on CX3 login node.
    # The cif must be pre-staged into $PEGRNA_DATA/afdb (compute nodes are offline);
    # stage_models.sh fetches it. Override accession via $PEGRNA_SMOKE_UNIPROT.
    import os

    from pegrna_pipeline.hpc.structure_importance import AFDBImportanceScorer

    acc = os.environ.get("PEGRNA_SMOKE_UNIPROT", "P0DP23")
    m = AFDBImportanceScorer(acc, afdb_dir=f"{data}/afdb").metrics(10)
    return (0.0 <= m.rsa <= 2.0 and 0.0 <= m.plddt <= 100.0, f"plddt={m.plddt:.1f} rsa={m.rsa:.2f}")


def smoke_rinalmo():
    import os

    from pegrna_pipeline.hpc.full_structure import RiNALMoScorer

    pp = RiNALMoScorer(weights=os.environ.get("RINALMO_SS_WEIGHTS")).pair_probs(_SEQ)
    return (isinstance(pp, dict), f"{len(pp)} pairs >= threshold")


def smoke_rnaformer():
    import os

    from pegrna_pipeline.hpc.full_structure import RNAformerScorer

    model_dir = os.path.dirname(os.environ.get("RNAFORMER_CKPT", "")) or None
    pp = RNAformerScorer(model_dir=model_dir).pair_probs(_SEQ)
    return (isinstance(pp, dict), f"{len(pp)} pairs >= 0.5")


def smoke_boltz(data):
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer

    # Point --cache at the pre-staged Boltz weights so the offline compute node
    # is self-contained (the fix for the "no confidence JSON" failure).
    iptm = Boltz2DockScorer(cache_dir=f"{data}/boltz_cache").interface_iptm(_RTT_RNA, "MKTAYIAKQR")
    return (0.0 <= iptm <= 1.0, f"iptm={iptm:.3f}")


def smoke_structure_auto(data):
    """The REAL design-path structure gate: RiNALMo (Ampere+ GPU) with a
    transparent RNAformer (CPU) fallback. Passes if EITHER backend runs, so a node
    without FlashAttention still validates via RNAformer."""
    import os

    from pegrna_pipeline.hpc.full_structure import AutoStructureScorer

    model_dir = os.path.dirname(os.environ.get("RNAFORMER_CKPT", "")) or f"{data}/rnaformer_weights"
    sc = AutoStructureScorer(
        rinalmo_weights=os.environ.get("RINALMO_SS_WEIGHTS"), rnaformer_model_dir=model_dir
    )
    pp = sc.pair_probs(_SEQ)
    return (isinstance(pp, dict), f"{len(pp)} pairs via backend={sc.backend_name}")


def smoke_casoffinder(data):
    import os

    from pegrna_pipeline.hpc.offtarget_scorer import CasOffinderScorer

    genome = os.environ.get("PEGRNA_GENOME") or f"{data}/GRCh38.fa"
    if not os.path.exists(genome):
        return (False, f"genome not found at {genome} (set PEGRNA_GENOME or stage GRCh38.fa)")
    sc = CasOffinderScorer(genome=genome, max_mismatches=2)
    res = sc(["GACGTTGACGATCGATCGAT"])  # arbitrary 20-mer; verifies the tool runs+parses
    on, off, _has = next(iter(res.values()))
    return (isinstance(res, dict) and len(res) == 1, f"ran; ontarget={on} offtargets={off}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=".", help="staging dir ($PEGRNA_DATA)")
    ap.add_argument(
        "--only",
        nargs="*",
        default=None,
        help="subset: splice importance rinalmo rnaformer boltz",
    )
    args = ap.parse_args(argv)

    # Default smoke = the 5 gates a real design run actually uses (structure via
    # the auto RiNALMo->RNAformer path). RiNALMo/RNAformer are individually
    # addressable via --only for debugging a specific backend.
    checks = {
        "splice (SpliceAI)": smoke_splice,
        "importance (AlphaFold/biotite)": lambda: smoke_importance(args.data),
        "structure (auto RiNALMo->RNAformer)": lambda: smoke_structure_auto(args.data),
        "docking (Boltz-2)": lambda: smoke_boltz(args.data),
        "offtarget (Cas-OFFinder)": lambda: smoke_casoffinder(args.data),
    }
    debug_checks = {
        "full_structure (RiNALMo)": smoke_rinalmo,
        "full_structure (RNAformer)": smoke_rnaformer,
    }
    keys = {
        "splice": "splice (SpliceAI)",
        "importance": "importance (AlphaFold/biotite)",
        "structure": "structure (auto RiNALMo->RNAformer)",
        "boltz": "docking (Boltz-2)",
        "offtarget": "offtarget (Cas-OFFinder)",
        "rinalmo": "full_structure (RiNALMo)",
        "rnaformer": "full_structure (RNAformer)",
    }
    all_checks = {**checks, **debug_checks}
    selected = [keys[k] for k in args.only] if args.only else list(checks)
    results = [_check(name, all_checks[name]) for name in selected]
    n_ok = sum(results)
    print(f"\n{n_ok}/{len(results)} gate models passed smoke validation")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
