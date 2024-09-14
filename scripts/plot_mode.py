"""Draw the mode profile of one structure.

python scripts/plot_mode.py cases/slot_gap100.wg mode.png
"""

import sys

from modegrid import engine, plotting, reference

if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    result = engine.run(sys.argv[1])
    odd = reference.read_structure(sys.argv[1])["symmetry"] == "odd"
    plotting.case_figure(result, result["title"], sys.argv[2], odd=odd)
    print(f"wrote {sys.argv[2]}")
