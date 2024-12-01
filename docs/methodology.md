# Methodology

ModeGrid finds the guided mode of highest effective index of a dielectric waveguide cross-section, for one polarisation and one symmetry class at a time, with a semi-vectorial finite-difference method.

## 1. Equation

For a field that is mostly along x (quasi-TE), neglecting the coupling to the minor component, the transverse wave equation is

    d/dx [ (1/eps) d(eps E_x)/dx ] + d2 E_x/dy2 + k^2 eps E_x = beta^2 E_x

and the same with x and y exchanged for a field mostly along y (quasi-TM). The scalar equation drops the permittivity from the derivatives. The effective index is beta / k.

The first term is what distinguishes the method from a scalar one: across an interface normal to the field, eps E is continuous and E itself jumps by the permittivity ratio. This is what makes slot waveguides work, and what lowers the quasi-TM index of a thin silicon strip far below the scalar value.

## 2. Discretisation

The window is divided into equal rectangular cells, each with one permittivity; the field is sampled at cell centres. Between a cell P and its neighbour Q along the direction of the dominant field, with the interface permittivity taken as the mean of the two,

    off-diagonal:  2 eps_Q / (eps_P + eps_Q) / h^2
    diagonal:     -2 eps_P / (eps_P + eps_Q) / h^2

and 1/h^2, -1/h^2 in the other direction (and in both directions for the scalar operator). In a uniform region the two forms coincide. The matrix is not symmetric where the permittivity changes, but its eigenvalues are real.

Rectangles are painted in order. A cell partly covered by a rectangle takes the area-weighted mean permittivity. This keeps any geometry usable, but the result then scatters from grid to grid; grids that put cell edges on the walls converge smoothly (section 6).

## 3. Boundaries

The outer walls are metallic: the field vanishes half a cell beyond the last centre. There are no absorbing layers, so only modes confined well inside the window are meaningful; the library checks that the index of a guided strip does not move by more than 5e-5 when the window is varied.

A mirror at x = 0 halves the problem. An even mirror selects the modes whose dominant field is symmetric in x (the fundamental mode, the even supermode of a coupler); an odd mirror selects the antisymmetric ones (the first higher-order mode, the odd supermode).

## 4. Eigenvalue

Every eigenvalue lies below k^2 eps_max. The solver shifts the matrix by slightly more than that, factorises it once as a banded LU (cells numbered along the shorter side of the grid) and iterates x <- (A - shift)^-1 x, which converges to the largest eigenvalue. After the estimate settles, the shift is moved once to a twentieth of its distance and the matrix refactorised; a typical solve takes 20 to 40 iterations. The factorisation has no pivoting: the shifted matrix is negative definite in the symmetric case and has behaved the same way in every structure of the library.

Only the mode of highest index in its symmetry class is found. Higher modes in the same class are out of reach of this iteration.

## 5. What is reported

- effective index, and whether the mode is guided: its index must exceed the cladding index. For uniform walls that is the highest material index on them. Where a slab reaches the side wall (rib waveguides), it is the index of that slab's own mode, solved on the same grid in one dimension. A mode below this value is a mode of the box, not of the waveguide.
- group index n - lambda dn/dlambda by a central difference over +-5 nm, with the material dispersion included through the linear material models.
- confinement: the share of the squared dominant field in the cells of highest index. This is a field measure, not a power-flow measure; for quasi-TM modes and slots it is much smaller than the power confinement usually quoted.
- effective area (integral E^2)^2 / integral E^4, root mean square half-widths and centroid.
- for couplers: the even and odd supermode indices and the coupling length lambda / (2 (n_even - n_odd)).

## 6. Accuracy

Discretisation. On grids aligned with the walls the index converges at first order in the cell size, from above. For a 480 x 240 nm silicon strip a 20 nm grid is 0.006 (quasi-TE) and 0.009 (quasi-TM) above the extrapolated value; a 10 nm grid halves that. Group index and coupling lengths are less sensitive than the index itself because the error largely cancels in differences.

Approximation. The semi-vectorial equation ignores the minor field component. For low index contrast (silicon nitride in oxide) this is a small effect. For silicon strips near the single-mode width the fields are strongly hybrid at the corners and a full-vector solver gives a quasi-TE index that is lower by a few hundredths. That difference does not shrink with the grid.

## 7. Validation

1. C++ unit tests: banded LU against a known solution, the eigen-iteration against the analytic spectrum of a second-difference matrix, the operator of a uniform box against its analytic eigenvalue for every polarisation and mirror, operator coefficients at an interface, and mode-level physics.
2. Slab limit: a film with no structure in x must give the exact slab index (lowered by the known lateral half cosine). Agreement on a 10 nm grid: 5e-5 for nitride, 1e-3 for 220 nm silicon TE, 1.9e-3 for 220 nm silicon TM.
3. The waveguide library: each file's expectations were computed by `reference.py`, which reads, paints, assembles and solves the same grid with SciPy's sparse eigen-solver, sharing no code with the C++ solver. The report pipeline repeats that comparison for every structure, including the field.
4. Consistency across structures: index rising with width, agreement with the effective index method, exponential decay of coupler splitting with gap, whole structure against mirrored half.

Items 3 and 4 verify the numerics and the implementation. Only item 2 compares with an exact result, and only in one dimension. There is no comparison with a full-vector solver or with measurements in this repository.
