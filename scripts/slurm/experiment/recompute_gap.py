#!/usr/bin/env python3
"""
recompute_gap.py - recompute a literature-comparable primal gap for existing
results.json files, using each instance's known optimal/best-known objective
value from a MIPLIB .solu file, instead of SCIP's own (weak, early) dual bound.

Formula (Definition 1.2, Mexi 2021 - "New Ideas for the Feasibility Pump:
Using Multiple Reference Vectors"):

           / 0,                                  if c^T x* = c^T x~ = 0
    gap = |  |c^T x* - c^T x~| / |c^T x*|,        if (c^T x*)*(c^T x~) > 0
           \\ inf,                                 otherwise

where x* is the known optimal/best-known solution (from the .solu file) and
x~ is the solution this run actually found (primalBound in results.json).

Writes two new fields into each results.json:
  "knownOptimal": the reference value looked up from the .solu file
  "gap": the primal gap per the formula above (a float, or the string "Inf"
         for not-found instances / the formula's own infinity case, matching
         the string-sentinel convention already used elsewhere in this data)

Does NOT touch dualBound (removed earlier, not part of this methodology).

Usage:
  ./recompute_gap.py <folder> [<folder> ...]
  e.g. ./recompute_gap.py run1/fpfw_baseline_s1 run1/fpfw_baseline_s2 ...
  or:  ./recompute_gap.py run1   (recurses into every *_s*/instance/ under it)
"""
import json
import os
import sys

PROJECT_DIR = os.environ.get("PROJECT_DIR", "/home/htc/aleoputra/project")
COMP_RESULT = os.path.join(PROJECT_DIR, "compResult")
FPFW_DIR = os.path.join(PROJECT_DIR, "FPmeetsFW")
# Bundled inside the repo (trimmed to the 240-instance MIPLIB 2017 Benchmark Set),
# same file the Julia source reads (src/constants.jl's SOLU_FILE) - single source
# of truth, and works for anyone who clones the repo, not just this cluster layout.
SOLU_FILE = os.path.join(FPFW_DIR, "data/miplib2017-v36.solu")


def load_solu(path):
    ref = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            status = parts[0].strip("=")
            name = parts[1]
            if status in ("opt", "best"):
                ref[name] = float(parts[2])
    return ref


def instance_name(fname):
    for suf in (".mps.gz", ".mps"):
        if fname.endswith(suf):
            return fname[: -len(suf)]
    return fname


def compute_gap(x_star, x_found):
    if x_star == 0.0 and x_found == 0.0:
        return 0.0
    if x_star * x_found > 0:
        return abs(x_star - x_found) / abs(x_star)
    return None  # infinite


def process_file(rj_path, ref):
    with open(rj_path) as f:
        d = json.load(f)

    inst = instance_name(d.get("instance", ""))
    x_star = ref.get(inst)
    if x_star is None:
        print(f"  warning: no reference value for instance '{inst}' ({rj_path})", file=sys.stderr)
        return False

    d["knownOptimal"] = x_star

    if not d.get("solutionFound", False):
        d["gap"] = "Inf"
    else:
        primal = d.get("primalBound")
        if primal is None or primal >= 1e19:
            d["gap"] = "Inf"
        else:
            g = compute_gap(x_star, primal)
            d["gap"] = g if g is not None else "Inf"

    with open(rj_path, "w") as f:
        json.dump(d, f, indent=2)
    return True


def find_results_json(target):
    full = os.path.join(COMP_RESULT, target)
    out = []
    for root, dirs, files in os.walk(full):
        if "results.json" in files:
            out.append(os.path.join(root, "results.json"))
    return out


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <folder> [<folder> ...]", file=sys.stderr)
        sys.exit(1)

    ref = load_solu(SOLU_FILE)
    print(f"Loaded {len(ref)} reference values from {SOLU_FILE}", file=sys.stderr)

    total = 0
    ok = 0
    for target in sys.argv[1:]:
        files = find_results_json(target)
        for rj in files:
            total += 1
            if process_file(rj, ref):
                ok += 1

    print(f"Processed {ok}/{total} results.json files", file=sys.stderr)


if __name__ == "__main__":
    main()
