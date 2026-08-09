"""Run probe.py against two S4 builds and diff the numbers.

Usage:

    python compare.py --a <path-to-build-a> --b <path-to-build-b> [test ...]

Each path is prepended to sys.path of a fresh child process, so the two builds
never share an interpreter.  Passing the empty string for a path means "use
whatever S4 the environment already has installed", which is how the old oracle
is reached.

Every test runs in its own subprocess.  A test that segfaults is reported as a
crash with its exit status instead of taking the run down.
"""

import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROBE = os.path.join(HERE, "probe.py")

# probe.py must not be run from inside the source tree: the repository has an
# S4/ directory that shadows the extension module as a namespace package.
RUNDIR = os.path.join(os.path.sep, "tmp")


def run_one(build, name, timeout=1800, python=None):
    env = dict(os.environ)
    pp = [HERE]
    if build:
        pp.insert(0, build)
    old = env.get("PYTHONPATH")
    env["PYTHONPATH"] = os.pathsep.join(pp + ([old] if old else []))
    try:
        p = subprocess.run([python or sys.executable, PROBE, "run", name],
                           cwd=RUNDIR, env=env, timeout=timeout,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.TimeoutExpired:
        return {"status": "timeout"}
    for line in p.stdout.decode("utf-8", "replace").splitlines():
        if line.startswith("@@JSON@@"):
            return json.loads(line[len("@@JSON@@"):])
    return {"status": "crash", "returncode": p.returncode,
            "stderr": p.stderr.decode("utf-8", "replace")[-2000:],
            "stdout": p.stdout.decode("utf-8", "replace")[-2000:]}


def decode(v):
    if isinstance(v, dict) and "__c__" in v:
        return complex(v["__c__"][0], v["__c__"][1])
    return v


def fmt(v):
    if isinstance(v, complex):
        return "%.12g%+.12gj" % (v.real, v.imag)
    return "%.12g" % v


def reldiff(a, b):
    d = abs(a - b)
    s = max(abs(a), abs(b))
    return d, (d / s if s > 0 else 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", default="", help="sys.path entry for build A (old oracle if empty)")
    ap.add_argument("--b", required=True, help="sys.path entry for build B")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--python-a", default=None,
                    help="interpreter to run build A under; defaults to this one. "
                         "Needed when the two builds target different Python "
                         "versions, or link different BLAS/FFT backends.")
    ap.add_argument("--python-b", default=None)
    ap.add_argument("--json", default=None, help="write the raw results here")
    ap.add_argument("--tol", type=float, default=1e-9,
                    help="relative tolerance for calling a value 'same'")
    ap.add_argument("--atol", type=float, default=1e-12,
                    help="absolute floor; below it a value counts as agreeing "
                         "regardless of relative difference (residuals that are "
                         "already at the rounding noise level)")
    ap.add_argument("tests", nargs="*")
    args = ap.parse_args()

    names = args.tests
    if not names:
        env = dict(os.environ)
        env["PYTHONPATH"] = HERE + os.pathsep + env.get("PYTHONPATH", "")
        out = subprocess.check_output([args.python_a or sys.executable,
                                       PROBE, "list"], cwd=RUNDIR, env=env)
        names = out.decode().split()

    raw = {}
    nsame = ndiff = nfail = 0
    for name in names:
        ra = run_one(args.a, name, python=args.python_a)
        rb = run_one(args.b, name, python=args.python_b)
        raw[name] = {"a": ra, "b": rb}

        print("=" * 78)
        print(name)
        print("-" * 78)
        if ra["status"] != "ok" or rb["status"] != "ok":
            nfail += 1
            for lbl, r in ((args.label_a, ra), (args.label_b, rb)):
                if r["status"] != "ok":
                    print("  %-6s %s" % (lbl, r["status"].upper()))
                    for k in ("returncode", "traceback", "stderr"):
                        if k in r and r[k]:
                            txt = str(r[k]).strip().splitlines()
                            for ln in txt[-6:]:
                                print("         | " + ln)
                else:
                    print("  %-6s ok" % lbl)
            # a build that ran gets its numbers printed anyway
            for lbl, r in ((args.label_a, ra), (args.label_b, rb)):
                if r["status"] == "ok":
                    for k, v in r["values"].items():
                        print("    %-28s %s" % (k, fmt(decode(v))))
            continue

        va = dict((k, decode(v)) for k, v in ra["values"].items())
        vb = dict((k, decode(v)) for k, v in rb["values"].items())
        keys = sorted(set(va) | set(vb))
        worst = 0.0
        print("  %-28s %-24s %-24s %s" % ("key", args.label_a, args.label_b, "rel"))
        for k in keys:
            if k not in va or k not in vb:
                print("  %-28s %-24s %-24s %s" % (
                    k, fmt(va[k]) if k in va else "-",
                    fmt(vb[k]) if k in vb else "-", "MISSING"))
                worst = float("inf")
                continue
            d, rd = reldiff(va[k], vb[k])
            eff = 0.0 if d <= args.atol else rd
            worst = max(worst, eff)
            flag = "" if eff <= args.tol else "  <<<"
            print("  %-28s %-24s %-24s %.2e%s" % (k, fmt(va[k]), fmt(vb[k]), rd, flag))
        if worst <= args.tol:
            nsame += 1
            print("  -> SAME (worst rel %.2e)" % worst)
        else:
            ndiff += 1
            print("  -> DIFFERS (worst rel %.2e)" % worst)

    print("=" * 78)
    print("same: %d   differs: %d   failed: %d" % (nsame, ndiff, nfail))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(raw, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
