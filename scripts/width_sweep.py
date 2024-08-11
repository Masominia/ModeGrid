"""Sweep the width of a buried strip and print index, group index and confinement.

python scripts/width_sweep.py [thickness_um] [te|tm]
"""

import sys
import tempfile
from pathlib import Path

from modegrid import engine

if __name__ == "__main__":
    thickness = float(sys.argv[1]) if len(sys.argv) > 1 else 0.22
    field = sys.argv[2] if len(sys.argv) > 2 else "te"
    print(f"{'width, nm':>10}{'n_eff':>10}{'group':>10}{'confined':>10}{'guided':>8}")
    with tempfile.TemporaryDirectory() as tmp:
        for width in range(300, 1001, 50):
            half = 0.5e-3 * width
            path = Path(tmp) / "strip.wg"
            path.write_text(
                f"grid dx=0.025 dy=0.02\nwindow width=3.0 ymin=-1.0 ymax={1.0 + thickness:g}\n"
                f"rect si {-half:g} {half:g} 0 {thickness:g}\nsymmetry even\nfield {field}\n"
            )
            m = engine.run(path, "check")["metrics"]
            guided = "yes" if m["guided"] else "no"
            print(f"{width:>10}{m['neff']:>10.4f}{m['group_index']:>10.4f}{m['confinement']:>10.3f}{guided:>8}")
