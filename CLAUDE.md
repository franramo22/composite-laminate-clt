# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

A MATLAB implementation of **Classical Laminate Theory (CLT)** for composite materials analysis. The tool interactively accepts material properties, layup definition, and applied loads, then computes the full laminate response and optionally renders 21 engineering visualizations.

## Running the analysis

From the MATLAB command window:

```matlab
res = laminate_analysis();       % interactive solver, prints results to console
laminate_viz(res);               % generate all 21 figures from a previous run
```

There is no build step, test suite, or package manager. All code runs directly in MATLAB R2018b or later (`yline`, `xline`, `sgtitle` are the version-limiting functions).

## File roles

| File | Role |
|---|---|
| `laminate_analysis.m` | Entry point. Prompts user for all inputs, runs the full CLT pipeline, prints tabulated results, returns `res` struct. |
| `laminate_viz.m` | Visualization-only module. Accepts the `res` struct and generates figures 1–21 as separate named windows. |
| `lamina.m` | Placeholder / scratch file (currently empty). |

## Architecture: the `res` struct

`laminate_analysis` returns a single struct that is the contract between the solver and the visualizer. Every field is in **SI units** internally. Key fields:

```
res.mat          % material properties struct (E1, E2, G12, nu12, nu21, t, rho, Xt, Xc, Yt, Yc, S)
res.theta_deg    % 1×n ply angles [deg], bottom→top
res.z            % 1×(n+1) ply boundary z-coordinates [m], measured from mid-plane
res.Q            % 3×3 reduced stiffness matrix [Pa]
res.Qbar         % 3×3×n transformed stiffness per ply [Pa]
res.A, B, D      % 3×3 extensional / coupling / bending stiffness [N/m, N, N·m]
res.ABD          % 6×6 full laminate stiffness
res.eps0, kap    % 3×1 midplane strains [-] and curvatures [1/m]
res.sig_g        % 3×npts global stresses at sampled z positions [Pa]
res.sig_l        % 3×npts local (fiber-frame) stresses [Pa]
res.FI_TW/TH/MS  % 1×n failure indices (Tsai-Wu, Tsai-Hill, Max Stress) at ply midplanes
res.delam_risk   % 1×(n-1) normalized delamination risk per interface [0–1]
res.N_cr, SF_buckle, m_cr  % buckling results
res.ply_*_pct    % 1×n percentage contribution of each ply to Ex, Ey, Gxy, D11
```

## CLT computation order (inside `laminate_analysis`)

1. Build `Q` from E1, E2, G12, nu12 (plane-stress reduced stiffness)
2. Build `Qbar(:,:,k)` per ply using the direct trigonometric formula (Jones 1999 convention — index 1=fiber, 2=transverse, 6=shear; matrix positions map as [1,1]→11, [1,2]→12, [2,2]→22, [3,3]→66)
3. Integrate through thickness to get A, B, D (using exact z-boundary differences, not midplane approximation)
4. Solve `ABD \ NMvec` for midplane strains `eps0` and curvatures `kap`
5. Evaluate strains/stresses at top and bottom face of each ply; de-duplicate shared boundaries
6. Transform global stresses to local frame with T = [c², s², 2cs; s², c², -2cs; −cs, cs, c²−s²]
7. Compute failure indices at ply midplane
8. Delamination risk = |Δθ| / 90°, capped at 1 (interface angle-change metric)
9. Buckling: minimize `(π/b)² [D11·R² + 2(D12+2D66) + D22/R²]` over m=1…10 half-waves, where R = mb/a

## Unit conventions (critical)

- All internal calculations use **SI** (Pa, m, N, N·m)
- User inputs are entered in **GPa / MPa / mm** and converted by `prompt_scalar` via its `scale` argument
- Printed outputs are converted back to GPa / MPa / mm / MN·m as labeled

## Visualization module structure (`laminate_viz.m`)

The main `laminate_viz(res)` function calls 21 local subfunctions named `viz_01_*` through `viz_21_*`. Each creates its own figure window. Two shared helpers at the top of the file:

- `ang_color(ang, all_angles)` — maps a ply angle to a consistent RGB color using `lines()` indexed by `unique(all_angles)`
- `risk_color(v)` — traffic-light RGB: green < 0.4 ≤ yellow < 0.7 ≤ red

The failure envelope (viz 19) solves the Tsai-Wu quadratic in σ₂ for a sweep of σ₁ values with τ₁₂ = 0. The polar stiffness (viz 18) uses the exact compliance rotation formula: `ā₁₁(θ) = a11·c⁴ + (2a12+a66)·c²s² + a22·s⁴ + 2a16·c³s + 2a26·cs³`.

## Saving changes to GitHub

```bash
git add -A
git commit -m "describe what changed"
git push
```

Remote: `https://github.com/franramo22/composite-laminate-clt`
