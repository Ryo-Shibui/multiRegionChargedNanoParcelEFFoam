# Bare gas gap: native steady/particle sequence

Run `./Allrun` in an OpenFOAM-v2406 environment. `steady/Allrun` solves the gas flow and all-region electrostatics for 2000 steady iterations, starting with zero interior velocity and a right-to-left parabolic inlet of mean .01 m/s. Pressure relaxation .3 and U relaxation .7 follow the validated native flow case. Species, energy and parcels are disabled.

`particles/Allrun` first runs `copySteadyFields.py`: all native fluid/solid output fields from steady time2000 are copied into particle physical time0. The particle run uses Euler dt=.0001 s and freezes flow/species/energy/electrostatics while `solveParcels true`. Native `reactingCloud1Properties` injects one +1e, 10 nm particle with Cunningham-corrected Stokes drag and Coulomb force. The carrier temperature is293.15 K. The gas gap is20mm.

Every rerun archives previous generated outputs before starting. Source inputs and cloud dictionaries remain editable. Deposition/trajectory CSV files are native cloud-function outputs under `particles/postProcessing/lagrangian/reactingCloud1/chargedParticleDiagnostics/fluid/processor0/`. Check results are saved as `acceptance.json` in each phase.
