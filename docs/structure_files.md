# Structure files

One statement per line. `#` starts a comment; the comment lines at the top of a file are its description and appear in the report. Lengths are in micrometres. The window is centred on x = 0.

## Statements

| statement | meaning | default |
|---|---|---|
| `title <text>` | name used in reports | empty |
| `wavelength <um>` | vacuum wavelength | 1.55 |
| `grid dx=<um> dy=<um>` | cell size | 0.02, 0.02 |
| `window width=<um> ymin=<um> ymax=<um>` | computational window; must be a whole number of cells | 3.0, -1.0, 1.22 |
| `background <material>` | what fills the window before painting | sio2 |
| `layer <material> <ymin> <ymax>` | a slab across the whole window | |
| `rect <material> <xmin> <xmax> <ymin> <ymax>` | a rectangle; later statements paint over earlier ones | |
| `symmetry none` / `even` / `odd` | mirror at x = 0; with a mirror only x >= 0 is solved | none |
| `field te` / `tm` / `scalar` | dominant field along x, along y, or no polarisation | te |
| `coupler` | solve the even and the odd supermode of a structure symmetric in x | off |
| `expect <kind> <args>` | a check, see below | |

Materials: `si`, `sio2`, `sin`, `air` (linear dispersion around 1.55 um) or a number for a constant index of at least 1.

With a mirror the structure must be symmetric in x; only the part at x >= 0 is used, and rectangles should still be written in full. `coupler` cannot be combined with `symmetry`.

Choose cell sizes that divide the dimensions of the waveguide where possible: cells cut by a wall take a mean permittivity and the result scatters by a few thousandths from grid to grid.

## Expectations

| kind | arguments | passes when |
|---|---|---|
| `neff` | value, tolerance | the effective index is within the tolerance |
| `neff_between` | low, high | the index lies in the range |
| `group_index` | value, tolerance | |
| `confinement` | fraction, tolerance | share of the squared field in the highest-index material |
| `area` | um^2, relative tolerance | effective mode area |
| `rms_x`, `rms_y` | um, tolerance | root mean square half-widths |
| `centre_y` | um, tolerance | height of the centroid |
| `guided` | `yes` or `no` | the index is above (or not above) the cladding index |
| `cells` | count | the grid has that many cells |
| `n_even`, `n_odd` | value, tolerance | supermode indices (needs `coupler`) |
| `splitting` | value, relative tolerance | n_even - n_odd |
| `coupling_length` | um, relative tolerance | |

A malformed expectation is reported as a failed check with its line number; the other checks still run.

## Output

`modegrid check <file>` prints `key value` lines and one `check PASS|FAIL <name> :: <detail>` line per expectation. Exit status: 0 if all checks pass, 1 if one fails, 2 if the file cannot be read or is invalid. `modegrid run <file>` adds one `eps` and one `field` line per grid row, bottom row first, with one value per cell.

## Example

    # A rib waveguide: 500 nm wide, etched down to a 90 nm slab.
    # Checked by a bound: the index lies between the two slab indices.
    title Silicon rib 500 nm wide on a 90 nm slab
    grid dx=0.025 dy=0.01
    window width=4.0 ymin=-1.0 ymax=1.22
    layer si 0 0.09
    rect si -0.25 0.25 0 0.22
    symmetry even
    field te
    expect neff_between 2.105 2.848
    expect guided yes
