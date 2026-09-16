"""
moderator_sweep.py
Sweeps moderator density and reruns the pin-cell model at each value, to
study the moderator density / void coefficient of reactivity.
"""

import os
import csv
import glob
import numpy as np
import matplotlib.pyplot as plt
import openmc

DENSITIES = [0.60, 0.65, 0.70, 0.72, 0.74, 0.76, 0.80, 0.85, 0.90]   # g/cm3
BASE_DIR = "moderator_sweep_runs"


def build_and_run(water_density, run_dir):
    os.makedirs(run_dir, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(run_dir)
    try:
        fuel = openmc.Material(name="UO2 fuel")
        fuel.add_element("U", 1.0, enrichment=3.2)
        fuel.add_element("O", 2.0)
        fuel.set_density("g/cm3", 10.4)

        clad = openmc.Material(name="Zircaloy-4 clad")
        clad.add_element("Zr", 0.98)
        clad.add_element("Sn", 0.015)
        clad.add_element("Fe", 0.002)
        clad.add_element("Cr", 0.001)
        clad.set_density("g/cm3", 6.56)

        water = openmc.Material(name="Borated water moderator")
        water.add_element("H", 2.0)
        water.add_element("O", 1.0)
        water.add_element("B", 700e-6)
        water.set_density("g/cm3", water_density)
        water.add_s_alpha_beta("c_H_in_H2O")

        openmc.Materials([fuel, clad, water]).export_to_xml()

        pellet_radius, clad_ir_r, clad_or_r, pitch = 0.4095, 0.418, 0.475, 1.26
        fuel_or = openmc.ZCylinder(r=pellet_radius)
        clad_ir = openmc.ZCylinder(r=clad_ir_r)
        clad_or = openmc.ZCylinder(r=clad_or_r)

        fuel_cell = openmc.Cell(fill=fuel, region=-fuel_or)
        gap_cell = openmc.Cell(region=+fuel_or & -clad_ir)
        clad_cell = openmc.Cell(fill=clad, region=+clad_ir & -clad_or)

        box = openmc.model.rectangular_prism(width=pitch, height=pitch, boundary_type="reflective")
        water_cell = openmc.Cell(fill=water, region=+clad_or & box)

        geometry = openmc.Geometry(openmc.Universe(cells=[fuel_cell, gap_cell, clad_cell, water_cell]))
        geometry.export_to_xml()

        settings = openmc.Settings()
        settings.batches = 120
        settings.inactive = 25
        settings.particles = 5000
        bounds = [-pitch / 2, -pitch / 2, -1, pitch / 2, pitch / 2, 1]
        settings.source = openmc.Source(
            space=openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
        )
        settings.export_to_xml()
        openmc.Tallies([]).export_to_xml()

        openmc.run(output=False)

        sp_file = sorted(glob.glob("statepoint.*.h5"))[-1]
        with openmc.StatePoint(sp_file) as sp:
            keff = sp.keff
        return keff.nominal_value, keff.std_dev
    finally:
        os.chdir(cwd)


def main():
    results = []
    for rho in DENSITIES:
        run_dir = os.path.join(BASE_DIR, f"rho_{rho:.2f}")
        print(f"Running moderator density = {rho} g/cm3 ...")
        keff, std = build_and_run(rho, run_dir)
        results.append((rho, keff, std))
        print(f"  k-eff = {keff:.5f} +/- {std:.5f}")

    with open("keff_vs_density.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["water_density_g_cm3", "k_eff", "k_eff_std_dev"])
        writer.writerows(results)

    rhos = np.array([r[0] for r in results])
    keffs = np.array([r[1] for r in results])
    stds = np.array([r[2] for r in results])
    slope, intercept = np.polyfit(rhos, keffs, 1)

    plt.figure(figsize=(6, 4.5))
    plt.errorbar(rhos, keffs, yerr=stds, fmt="o", capsize=3, label="OpenMC results")
    fit_x = np.linspace(rhos.min(), rhos.max(), 50)
    plt.plot(fit_x, slope * fit_x + intercept, "--", label=f"linear fit (slope={slope:.4f})")
    plt.xlabel("Moderator density (g/cm3)")
    plt.ylabel("k-eff")
    plt.title("Moderator density coefficient (pin-cell model)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("keff_vs_density.png", dpi=150)

    print("\nWrote keff_vs_density.csv and keff_vs_density.png")
    print(f"d(k-eff)/d(density) ~= {slope:.4f} per g/cm3")


if __name__ == "__main__":
    main()