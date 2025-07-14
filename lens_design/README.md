# Time-Reversal Acoustic Lens Design with k-Wave (MATLAB)

This repository contains a MATLAB script for designing custom acoustic lenses using **time-reversal wave propagation** via the **k-Wave** toolbox. The approach allows you to focus ultrasound energy at an arbitrary point by simulating a point source and recording the arrival times across a virtual transducer surface, which are then converted into a phase-compensated lens thickness map.

## 🔍 Overview

The simulation workflow:
1. Defines a 3D grid and homogeneous propagation medium
2. Places a point source at the target focal point
3. Records the time-domain pressure response at the source plane
4. Time-reverses the recorded signals and replays them into the medium
5. Visualises the focal region via axial slices and 3D isosurfaces
6. Extracts arrival times and converts them into a **lens thickness map**
7. Exports a circular acoustic lens design as both a **2D PNG heightmap** and a **3D surface mesh**

## 📂 Files

- `lens_design_from_timereversal.m`: Full pipeline for simulation, time reversal, and lens generation
- Output:  
  - `lens_heightmap.png`: Full 2D grayscale PNG of lens thickness  
  - `lens_heightmap_44mm_circle.png`: Circular 44 mm lens version (cropped and masked)  
  - Visualisation figures: axial/sagittal slices, beam profile, and 3D surfaces

## ⚙️ Dependencies

- [k-Wave Toolbox](https://www.k-wave.org/)  
- GPU-accelerated binary `kspaceFirstOrder3DG`  
- MATLAB R2023a or later recommended  
- No Image Processing Toolbox required (manual grayscale conversion used)

## 🧠 Applications

- Acoustic lens design for focused ultrasound neuromodulation
- Patient-specific transducer shaping
- Phase-compensation strategies for skull/heterogeneous media (future extension)
- Educational tool for understanding time-reversal focusing

## 📈 Visual Output

- **Axial and sagittal pressure slices**
- **3D isosurface of beam focus**
- **Axial line profile of beam**
- **Lens thickness map (visual and PNG)**
- **3D surface render of final lens**

## 📝 Notes

- The script assumes a circular 44 mm aperture with 1 mm grid resolution  
- Adjust `focus`, `f0`, `lens diameter`, and `lens material speed` to suit your setup  
- Lens material speed (`c_lens`) is assumed to be 2500 m/s by default  
- Scaling to a 5 mm or 20 mm lens height is configurable

---

### Author

Stephen Hicks | stephenlhicks@gmail.com  
If you use this code or adapt it in research, please consider citing or acknowledging this repository.
