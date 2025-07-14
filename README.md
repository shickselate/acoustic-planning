# Acoustic Simulation & Lens Design Tools (MATLAB + k-Wave)

This repository provides two simulation pipelines for focused ultrasound research using MATLAB and the [k-Wave Toolbox](https://www.k-wave.org/):

- `beam_simulation/`: Defines a circular transducer mask and simulates phase-delayed sources to produce a focused ultrasound beam.
- `lens_design_timereversal/`: Uses time-reversal to generate a phase-compensated lens design for arbitrary focal points.

Each module has its own README and is independently runnable. Visual outputs include axial pressure maps, beam profiles, and 3D isosurfaces.

## Requirements

- MATLAB R2022b or later
- k-Wave Toolbox
- Compatible GPU + CUDA (for `kspaceFirstOrder3DG`)
