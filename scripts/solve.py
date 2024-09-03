"""Solve one structure file and print what the solver found.

python scripts/solve.py cases/rib_slab90.wg
"""

import sys

from modegrid import engine

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    r = engine.run(sys.argv[1], "check")
    m = r["metrics"]
    print(r["title"])
    print(
        f"  grid {int(m['nx'])} x {int(m['ny'])} cells of {1e3 * m['dx']:g} x {1e3 * m['dy']:g} nm, {int(m['iterations'])} iterations"
    )
    print(
        f"  n_eff {m['neff']:.5f}, group index {m['group_index']:.4f}, cladding {m['n_cladding']:.4f}, {'guided' if m['guided'] else 'not guided'}"
    )
    print(
        f"  confinement {m['confinement']:.3f}, area {m['area']:.4f} um2, rms size {m['rms_x']:.3f} x {m['rms_y']:.3f} um"
    )
    if "coupling_length" in m:
        print(f"  supermodes {m['n_even']:.5f} and {m['n_odd']:.5f}, coupling length {m['coupling_length']:.2f} um")
    for c in r["checks"]:
        print(f"  {c['status']}  {c['check']}  ({c['detail']})")
