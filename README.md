# Neutronics-Analysis-of-a-Simplified-Thermal-Reactor-Using-OpenMC-and-Python-
A small reactor-physics project built in Python, using OpenMC's Monte Carlo neutron transport engine to model a single PWR-style fuel pin and pull out three numbers that actually matter in reactor physics: the multiplication factor, the neutron flux spectrum by region, and the moderator density coefficient of reactivity.

I put this together as a second computational project to go alongside my CFD work — same underlying idea (build a model, run a solver, check the output against what theory predicts), just applied to neutron transport instead of fluid flow. It's also a fairly direct match for the kind of work PSI's heat-pipe microreactor projects list as core skills: Python, Linux, and OpenMC, more or less end to end.

## Table of Contents
- [What Was Actually Done](#what-was-actually-done)
- [Setup / Software Used](#setup--software-used)
- [Reproducing the Results](#reproducing-the-results)
- [Results](#results)
- [How the Scripts Work](#how-the-scripts-work)
- [Limitations](#limitations)
- [Repository Structure](#repository-structure)
- [Author](#author)

## What was actually done

- Built a single UO2 / Zircaloy-4 / borated-water pin-cell in OpenMC's Python API, on a standard 1.26 cm PWR pitch, with reflective boundaries (the usual way of approximating an infinite lattice with one unit cell)
- Ran a base-case criticality calculation (10,000 particles/batch, 150 batches, 30 inactive) and got k-eff = 1.30949 ± 0.00083
- Tallied flux by energy — 50 log-spaced bins from 10⁻⁵ eV to 2×10⁷ eV — separately in the fuel, clad, and moderator cells, then converted it to flux-per-unit-lethargy so the whole spectrum plots sensibly on one log-log axis
- Noticed, and actually reported rather than hid, the self-shielding gap between the fuel and moderator flux curves in the thermal/epithermal range
- Swept moderator density across nine values (0.60–0.90 g/cm³) to get the sign of the moderator coefficient, dropping the particle count to 5,000/batch for the sweep since it meant nine separate full criticality runs instead of one
- Caught a small mismatch — the sweep's k-eff at 0.74 g/cm³ (1.30772) didn't land exactly on the base case's k-eff at the same density (1.30949) — traced it to the particle-count difference between the two runs rather than just assuming the model was broken, and confirmed the two values agree within about two standard deviations

## Setup / Software Used

* OpenMC, installed via `conda create -n openmc-env -c conda-forge openmc`, inside WSL2 running Ubuntu 24.04 (OpenMC doesn't have a native Windows build)
* ENDF/B-VII.1 (NNDC) cross-section and thermal-scattering data, pulled in through the `openmc_data_downloader` package
* Python 3.13, with NumPy, Matplotlib, and Pandas doing the post-processing work

Honestly, getting the environment running ate up more time than the actual modelling did. WSL's networking had a bug where any large HTTPS download — conda packages, the nuclear data files, didn't matter which — kept dying partway through with SSL errors. Turned out to be a known WSL2 networking issue, fixed by switching to mirrored networking mode (adding `networkingMode=mirrored` under `[wsl2]` in `.wslconfig`). Worth mentioning here since anyone else setting this up on Windows will probably hit the same wall.

## Reproducing the Results

1. `conda create -n openmc-env -c conda-forge openmc`, then `conda activate openmc-env`
2. `pip install openmc_data_downloader`, then run it against an `openmc.Materials` object matching the fuel/clad/water compositions in the report (this also points `OPENMC_CROSS_SECTIONS` at the right file automatically)
3. `python pin_cell_model.py` — builds the model and runs the base-case calculation
4. `python analyze_results.py` — prints k-eff, writes `flux_spectrum.png`
5. `python moderator_sweep.py` — runs the nine-point sweep, writes `keff_vs_density.png` and `.csv` (takes a good few minutes — nine separate simulations)

## Results

**Base case, moderator density 0.74 g/cm³:** k-eff = 1.30949 ± 0.00083. Since this is a reflective-boundary pin-cell, there's zero neutron leakage by construction, so what's actually being computed is k-infinity, not k-effective. A value around 1.3 for fresh 3.2%-enriched fuel is exactly what you'd expect from this kind of infinite-lattice setup.

**Flux spectrum:** clean two-hump shape — fast peak around 10⁵–10⁷ eV, thermal peak below ~1 eV, dip in between as neutrons slow down. The fuel curve sits visibly below the moderator curve through the thermal/epithermal range — that's neutron self-shielding, where the outer layer of fuel absorbs resonance-energy neutrons before they reach the pellet interior.

**Moderator density sweep:** k-eff vs. density comes out linear — slope = +0.1478 per g/cm³ (`results/keff_vs_density.png`, `.csv`). Positive slope means losing coolant density (boiling, voiding) drops k-eff, which is a negative moderator coefficient — the safety-favourable direction for a thermal reactor.

| Moderator density (g/cm³) | k-eff | σ |
|---|---|---|
| 0.60 | 1.28145 | 0.00127 |
| 0.65 | 1.29183 | 0.00135 |
| 0.70 | 1.30040 | 0.00133 |
| 0.72 | 1.30335 | 0.00139 |
| 0.74 | 1.30772 | 0.00145 |
| 0.76 | 1.31122 | 0.00132 |
| 0.80 | 1.31689 | 0.00144 |
| 0.85 | 1.32180 | 0.00139 |
| 0.90 | 1.32505 | 0.00132 |

## How the scripts work

- **`pin_cell_model.py`** — defines the three materials (UO2 fuel, Zircaloy-4 clad, borated water with the `c_H_in_H2O` thermal scattering law attached), builds the pin-cell geometry out of concentric `ZCylinder` surfaces plus a reflective bounding box, sets the run to 10,000 particles/batch over 150 batches (30 inactive), and defines the energy-binned flux tally. Running it kicks off the base-case criticality calculation directly.
- **`analyze_results.py`** — opens the statepoint HDF5 file OpenMC writes at the end of a run, reads out `k-eff` and its uncertainty, pulls the flux tally as a pandas dataframe, and does the flux-to-per-unit-lethargy conversion before plotting — without that conversion the fast peak visually swamps the thermal region and the plot is basically useless.
- **`moderator_sweep.py`** — the same model, rebuilt from scratch nine times at nine different moderator densities, each run isolated in its own subdirectory so the XML files OpenMC writes per run don't collide. Uses a lighter particle count than the base case purely to keep total runtime sane across nine full simulations.

## Limitations

This is one infinite pin-cell, not a fuel assembly or a full core — it has nothing to say about neutron leakage at a real core boundary, spatial power peaking, or control-rod worth, and its k-eff can't be compared directly to a real reactor's operating k-eff (which sits near 1.0 because of leakage and poisons this model doesn't include). The fuel is fresh — zero burnup, no Xe-135, no depletion of U-235 — so this doesn't capture how k-eff falls across an actual fuel cycle. I also didn't run a formal source-convergence check (tracking Shannon entropy of the fission source across inactive batches, for instance) beyond trusting OpenMC's default settings — the uncertainties quoted are plain Monte Carlo statistical uncertainty, not a full uncertainty budget.

If I come back to this, the obvious next steps are sweeping enrichment or pin pitch the same way the moderator density was swept, or hooking in `openmc.deplete` to watch k-eff decline over a burnup cycle instead of only ever looking at fresh fuel.

## Repository structure

```
├── README.md
├── report/
│   └── Neutronics_Pin_Cell_Report.docx
├── pin_cell_model.py        → base-case geometry, materials, settings, tallies
├── analyze_results.py       → extracts k-eff and plots the flux spectrum
├── moderator_sweep.py       → nine-point moderator density sweep + linear fit
└── results/
    ├── flux_spectrum.png
    ├── keff_vs_density.png
    └── keff_vs_density.csv
```

## Author

**Pratyush Dash**

B.Tech Chemical Engineering, KIIT University, Bhubaneswar
