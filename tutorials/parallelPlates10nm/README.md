# Parallel plates with and without a dielectric

Run `./Allrun` in this directory after sourcing OpenFOAM. Both variants first
solve steady gas flow/electrostatics, then copy the converged carrier fields to a
separate transient particle case. The scripts compare flow, electric fields,
frozen fields and wall-deposition time to their independent expectations.

The gas gap is 20 mm; one electrode carries a 1 mm dielectric with relative
permittivity 4 in `withDielectric`. Voltage is 1000 V. Mean gas speed is
0.01 m/s from right to left. A 10 nm particle with charge +1e starts halfway
across the gas gap. Stokes-Cunningham drift uses a mean free path of 65 nm.
Brownian motion is disabled for this deterministic verification.

Theory gas fields: 50000 and 49382.7160494 V/m. Theory deposition times:
0.0962724 and 0.0974758 s. These are idealized verification cases, not a
validation of material data for a particular device.
