"""
analyze_results.py
Reads the statepoint file produced by pin_cell_model.py, prints k-eff with
its statistical uncertainty, and plots the neutron flux spectrum by region.
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
import openmc

statepoints = sorted(glob.glob("statepoint.*.h5"))
if not statepoints:
    raise FileNotFoundError("No statepoint.*.h5 found -- run pin_cell_model.py first.")

with openmc.StatePoint(statepoints[-1]) as sp:
    keff = sp.keff
    print(f"k-eff = {keff.nominal_value:.5f} +/- {keff.std_dev:.5f}")

    tally = sp.get_tally(name="flux_by_region_and_energy")
    df = tally.get_pandas_dataframe()

energy_low = sorted(df["energy low [eV]"].unique())
energy_high = sorted(df["energy high [eV]"].unique())
regions = df["cell"].unique()

fig, ax = plt.subplots(figsize=(7, 5))

for region in regions:
    sub = df[df["cell"] == region].sort_values("energy low [eV]")
    e_lo = sub["energy low [eV]"].values
    e_hi = sub["energy high [eV]"].values
    flux = sub["mean"].values

    lethargy_width = np.log(e_hi / e_lo)
    flux_per_lethargy = flux / lethargy_width
    e_mid = np.sqrt(e_lo * e_hi)

    ax.step(e_mid, flux_per_lethargy, where="mid", label=f"cell {region}")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Energy (eV)")
ax.set_ylabel("Flux per unit lethargy (a.u.)")
ax.set_title("Neutron flux spectrum by region")
ax.legend()
ax.grid(alpha=0.3, which="both")

plt.tight_layout()
plt.savefig("flux_spectrum.png", dpi=150)
print("Wrote flux_spectrum.png")