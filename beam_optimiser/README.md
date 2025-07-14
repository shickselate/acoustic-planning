
# Iterative Phase Optimisation for Acoustic Focusing

This MATLAB script performs iterative optimisation of a phase map to maximise the focus-to-sidelobe ratio (FSLR) in a simulated 3D acoustic field using the k-Wave toolbox.

## 🧠 Purpose

The goal is to improve the precision of ultrasound energy delivery by tuning the phase profile of a planar source. The process enhances pressure at a defined focal point while suppressing surrounding sidelobes.

## 🔁 Workflow Summary

1. Initialise a planar transducer with apodisation and zero phase.
2. Define a 3×3×3 voxel focus region within a 3D grid.
3. Iteratively update the phase map using gradient-based reinforcement:
   - Run 3D simulation (CUDA-accelerated via `kspaceFirstOrder3DG`)
   - Compute pressure field and FSLR
   - Adjust phases to boost contribution to focus
4. Save best-performing phase map and pressure field.

## 📊 Outputs

- **FSLR plot** across iterations (in dB)
- **Best axial slice** through the focus
- **Phase map** yielding the highest FSLR
- `.mat` file containing:
  - `best_phases`: 2D phase map
  - `best_p_map`: 3D pressure distribution
  - simulation grid info

## ⚙️ Parameters

- Grid size: 44×44×96 (1 mm resolution)
- Ultrasound freq: 500 kHz
- Medium speed: 1500 m/s
- Focus: voxel `[22, 22, 70]`
- Iterations: 10
- Tukey window apodisation
- Learning rate: 0.1 (phase step scale)

## 📁 Files

- `best_beam_iterationX.mat`: best result
- Figures showing optimisation performance

## 📝 Notes

- The optimisation is local and may benefit from more iterations or alternative gradient methods.
- Requires k-Wave with CUDA support.
- All units in millimetres unless otherwise stated.
