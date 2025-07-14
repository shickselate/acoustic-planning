# Focused Ultrasound Beam Simulation in MATLAB

This repository contains a MATLAB script for simulating 3D focused ultrasound beams using the **k-Wave** toolbox and GPU-accelerated time-domain wave propagation (`kspaceFirstOrder3DG`). The simulation models a circular planar transducer emitting time-delayed, apodised sinusoidal waveforms to focus energy at a specified point in space.

## 🔍 Overview

The simulation:
- Defines a 3D computational grid and homogeneous medium
- Models a circular transducer using a binary mask
- Applies spatial apodisation (Tukey window)
- Delays the phase of each source element to focus energy at a given 3D focal point
- Computes the resulting pressure field over time using a GPU-accelerated solver
- Visualises the beam profile in axial, sagittal, and coronal slices, and as a 3D isosurface

## 📂 Files

- `simulate_focused_beam.m`: Main script that runs the simulation and generates visualisations
- `myTukeywin.m`: Custom Tukey window generator used for spatial apodisation
- Output: Pressure fields, 2D slices, and a 3D isosurface of the focal region

## ⚙️ Dependencies

- [k-Wave Toolbox](https://www.k-wave.org/)
- GPU-compatible version of `kspaceFirstOrder3DG`
- MATLAB (tested with R2023a or later)
- Compatible GPU + CUDA (for GPU acceleration)

## 🧠 Applications

- Focused ultrasound neuromodulation
- Transducer array simulation
- Acoustic lens design validation
- 3D beam characterisation and targeting

## 📊 Visual Output

- Axial, sagittal, and coronal pressure maps
- Adaptive 3D isosurface of high-pressure region

## 📝 Notes

- The focus point, grid size, and transducer aperture can be adjusted in the script
- Apodisation improves beam quality and reduces sidelobes
- This script does not model skull heterogeneity or lens phase correction

## 📈 Example Visuals

_(optional: add screenshots here)_

---

### Author

Stephen Hicks | stephenlhicks@gmail.com  
If you use this code or adapt it in research, please consider citing or acknowledging this repository.
