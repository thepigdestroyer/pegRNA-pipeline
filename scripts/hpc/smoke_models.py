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
    # Use a small well-characterised UniProt (ubiquitin P0CG48 fragment) if online,
    # else a locally staged AFDB cif under $PEGRNA_DATA/afdb.
    from pegrna_pipeline.hpc.structure_importance import AFDBImportanceScorer

    m = AFDBImportanceScorer("P0DTC2", afdb_dir=f"{data}/afdb").metrics(0)
    return (0.0 <= m.rsa <= 2.0 and 0.0 <= m.plddt <= 100.0, f"plddt={m.plddt:.1f} rsa={m.rsa:.2f}")


def smoke_rinalmo():
    from pegrna_pipeline.hpc.full_structure import RiNALMoScorer

    pp = RiNALMoScorer().pair_probs(_SEQ)
    return (isinstance(pp, dict), f"{len(pp)} pairs >= threshold")


def smoke_rnaformer():
    from pegrna_pipeline.hpc.full_structure import RNAformerScorer

    pp = RNAformerScorer().pair_probs(_SEQ)
    return (isinstance(pp, dict), f"{len(pp)} pairs >= 0.5")


def smoke_boltz():
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer

    iptm = Boltz2DockScorer().interface_iptm(_RTT_RNA, "MKTAYIAKQR")
    return (0.0 <= iptm <= 1.0, f"iptm={iptm:.3f}")


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

    checks = {
        "splice (SpliceAI)": smoke_splice,
        "importance (AlphaFold/biotite)": lambda: smoke_importance(args.data),
        "full_structure (RiNALMo)": smoke_rinalmo,
        "full_structure (RNAformer)": smoke_rnaformer,
        "docking (Boltz-2)": smoke_boltz,
    }
    keys = {
        "splice": "splice (SpliceAI)",
        "importance": "importance (AlphaFold/biotite)",
        "rinalmo": "full_structure (RiNALMo)",
        "rnaformer": "full_structure (RNAformer)",
        "boltz": "docking (Boltz-2)",
    }
    selected = [keys[k] for k in args.only] if args.only else list(checks)
    results = [_check(name, checks[name]) for name in selected]
    n_ok = sum(results)
    print(f"\n{n_ok}/{len(results)} gate models passed smoke validation")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
