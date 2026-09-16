"""
pin_cell_model.py
A simplified thermal (PWR-style) pin-cell neutronics model in OpenMC.

Geometry: single UO2 fuel rod + gap + Zircaloy cladding, surrounded by
borated water moderator on a square pitch, with reflective boundaries
(this approximates an infinite lattice -- standard pin-cell simplification
of a full reactor core).

Computes:
    - k-eff (multiplication factor)
    - energy-binned neutron flux, tallied separately in fuel / clad / moderator

Requires: `pip install openmc` (or conda-forge) AND a cross-section data
library on disk with OPENMC_CROSS_SECTIONS pointing at its cross_sections.xml
-- see README.md for exactly how to get this set up.

Run:
    python pin_cell_model.py
"""

import numpy as np
import openmc

# ---------------------------------------------------------------------------
# 1. MATERIALS
# ---------------------------------------------------------------------------
fuel = openmc.Material(name="UO2 fuel")
fuel.add_element("U", 1.0, enrichment=3.2)   # 3.2 w/o U-235, typical PWR
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
water.add_element("B", 700e-6)      # 700 ppm natural boron -- typical PWR
water.set_density("g/cm3", 0.740)   # hot-operating density; sweep this in
                                     # moderator_sweep.py to study the
                                     # moderator density (void) coefficient
water.add_s_alpha_beta("c_H_in_H2O")

materials = openmc.Materials([fuel, clad, water])
materials.export_to_xml()

# ---------------------------------------------------------------------------
# 2. GEOMETRY (typical 17x17 PWR assembly pin-cell dimensions, cm)
# ---------------------------------------------------------------------------
pellet_radius = 0.4095
clad_inner_radius = 0.418
clad_outer_radius = 0.475
pitch = 1.26   # square pin pitch

fuel_or = openmc.ZCylinder(r=pellet_radius)
clad_ir = openmc.ZCylinder(r=clad_inner_radius)
clad_or = openmc.ZCylinder(r=clad_outer_radius)

fuel_cell = openmc.Cell(name="fuel", fill=fuel, region=-fuel_or)
gap_cell = openmc.Cell(name="gap", region=+fuel_or & -clad_ir)          # void
clad_cell = openmc.Cell(name="clad", fill=clad, region=+clad_ir & -clad_or)

box = openmc.model.rectangular_prism(width=pitch, height=pitch, boundary_type="reflective")
water_cell = openmc.Cell(name="moderator", fill=water, region=+clad_or & box)

root_universe = openmc.Universe(cells=[fuel_cell, gap_cell, clad_cell, water_cell])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# ---------------------------------------------------------------------------
# 3. SETTINGS
# ---------------------------------------------------------------------------
settings = openmc.Settings()
settings.batches = 150
settings.inactive = 30       # discard first 30 batches while fission source converges
settings.particles = 10000

bounds = [-pitch / 2, -pitch / 2, -1, pitch / 2, pitch / 2, 1]
settings.source = openmc.Source(
    space=openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
)
settings.export_to_xml()

# ---------------------------------------------------------------------------
# 4. TALLIES -- energy-binned flux, by region (fuel / clad / moderator)
# ---------------------------------------------------------------------------
# Log-spaced bins from 1e-5 eV (cold thermal) to 20 MeV (fast fission
# spectrum tail) -- 50 bins is enough to see the thermal peak, epithermal
# resonance dip, and fast peak without needing a named group structure.
energy_bins = np.logspace(-5, np.log10(2.0e7), 50)
energy_filter = openmc.EnergyFilter(energy_bins)
cell_filter = openmc.CellFilter([fuel_cell, clad_cell, water_cell])

flux_tally = openmc.Tally(name="flux_by_region_and_energy")
flux_tally.filters = [cell_filter, energy_filter]
flux_tally.scores = ["flux"]

tallies = openmc.Tallies([flux_tally])
tallies.export_to_xml()

# ---------------------------------------------------------------------------
# 5. RUN
# ---------------------------------------------------------------------------
openmc.run()

print("\nDone. k-eff is printed above and stored in the statepoint file.")
print("Run analyze_results.py next to extract k-eff + plot the flux spectrum.")
