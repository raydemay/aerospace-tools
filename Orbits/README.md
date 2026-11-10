# Orbits

A collection of orbital mechanics and trajectory tools.

## J2_Perturbations.py / J2Test.m
The MATLAB file is a derivative of two-body orbit determination that was used for Aero 5626, with J2 perturbations included in the differential equation. The Python script replicates that MATLAB script.

## plot_spk_trajectory.py
A 3D trajectory visualizer for NASA SPICE SPK ephemeris files (`.bsp` / `.spk`) using **SpiceyPy**, **Matplotlib**, and **Plotly**.

### Prerequisites
```bash
pip install spiceypy matplotlib plotly numpy
```

### Features
- Reads binary SPICE SPK files and automatically detects target bodies and coverage windows.
- Computes 3D trajectory in heliocentric ecliptic coordinates (`ECLIPJ2000`).
- Plots:
  - Spacecraft / body trajectory with labeled start and end locations
  - Earth's orbital path (approx. 1 AU reference orbit around the Sun)
  - Sun marker at the heliocentric origin
- Renders via both **Matplotlib** (static 3D window and PNG export) and **Plotly** (interactive HTML with 3D rotation, zooming, and date hover tooltips).

### Usage
```bash
# Interactive prompt
python plot_spk_trajectory.py

# Specify SPK kernel directly
python plot_spk_trajectory.py path/to/trajectory.bsp

# Specify target body and renderer
python plot_spk_trajectory.py path/to/trajectory.bsp --target -98 --renderer plotly --output-html new_horizons.html
```
