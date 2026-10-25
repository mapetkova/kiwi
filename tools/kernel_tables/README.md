# kernel_tables: derivation and verification of the thesis kernel-integral coefficients

Standard-library Python only (exact `fractions` arithmetic, hand-written Gauss-Legendre quadrature).

| file | purpose |
|---|---|
| `kernels.py` | the six kernels of thesis appendix A (M4, M5, M6, Wendland C2/C4/C6): exact pieces, breakpoints, normalisation |
| `derive.py` | derives the appendix-A coefficient tables (2D `a_p`; 3D `a_i` and `B(q0)`) from the kernel definition alone |
| `parse_thesis.py` | parses the LaTeX tables of `appendix_other_kernels.tex` into `thesis_tables.json` |
| `compare_tables.py` | thesis tables vs derived tables, entry by entry |
| `reference.py` | brute-force quadrature of the vertex-pyramid integral `I_P` straight from the kernel definition |
| `verify.py` | table-based `I_P` vs the brute-force reference (run by `ctest` as `kernel_tables`) |

```
python3 parse_thesis.py /path/to/appendix_other_kernels.tex   # regenerates thesis_tables.json
python3 compare_tables.py                                      # lists disagreements with the thesis
python3 verify.py 8                                            # numerical check, 8 random points per kernel
```

## Results

The derived tables agree with brute-force quadrature to about 1e-12 (relative) for all six kernels,
in 2D and 3D. The thesis tables agree with the derived ones except for three sign errors:

1. 3D cubic spline, `B_3`, `q <= 1/2`: `b_1` printed as `-7/160`, should be `+7/160`.
2. 3D quintic spline, `B_4`, `1/3 < q <= 2/3`: `b_1` printed as `-473/(3^7 14)`, should be `+473/(3^7 14)`.
3. 2D Wendland C4, `q <= 1`: `a_10` printed as `-7/6`, should be `+7/6`.

Notation: appendix eq. for the 3D `I_P` reads `I_P = A sum_i a_i I_i q0^(1-i)`; for `i = 0` the term
is `B_j I_0` with no factor `q0` (the tabulated `b_k` already contain every power of `q0`).

Numerical note: `verify.py` integrates over phi with 20-point Gauss-Legendre per segment, which is
the limiting accuracy at small `r0` (the integrand is sharply peaked there); with 40 points the same
cases agree to about 3e-12.  The closed-form table evaluation in C++ should still be checked for
cancellation at small `r0` and for the high-order kernels (Wendland C6 reaches `I_-13`).
