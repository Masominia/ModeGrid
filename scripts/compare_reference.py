"""Solve one structure again with SciPy and print the differences.

python scripts/compare_reference.py cases/coupler_te_gap200.wg
"""

import sys

from modegrid import engine
from modegrid.compare import compare

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    result = engine.run(sys.argv[1])
    cmp = compare(sys.argv[1], result)
    print(f"{'quantity':<18}{'solver':>16}{'difference':>14}")
    for name, diff in cmp["differences"].items():
        print(f"{name:<18}{result['metrics'][name]:>16.9f}{diff:>14.1e}")
    print(f"field differs by at most {cmp['field']:.1e} of its peak")
    sys.exit(1 if cmp["outside"] else 0)
