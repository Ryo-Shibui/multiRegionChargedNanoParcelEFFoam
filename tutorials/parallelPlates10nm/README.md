# Dielectric-coated parallel plates: charged nanoparticle deposition

This tutorial answers a physical question: how does a thin dielectric on the collecting electrode change the electric field and the deposition time of a singly charged 10 nm particle? Solve the carrier field first, then freeze it while tracking one particle. Compare the coated electrode with an otherwise identical bare electrode.

## Run

Load the OpenCFD OpenFOAM environment containing the compiled solver and its libraries, then run:

~~~bash
cd tutorials/parallelPlates10nm
./Allrun
~~~

For example, a v2406 installation may be loaded with source ~/OpenFOAM/OpenFOAM-v2406/etc/bashrc. NumPy is required for verification. Matplotlib produces standard PNG figures when installed; without it, the CSV comparisons and all numerical checks still run. No helper imports data or code from an external research workspace.

Allrun runs both comparisons, each in two phases:

1. noDielectric/steady and withDielectric/steady use native blockMesh -region ..., solve the gas flow and electrostatics for 2000 steady iterations, and write the solved fields.
2. Each particles case copies the corresponding native mesh and solved fields into time zero, disables carrier equations, and tracks the particle from 0 to 0.11 s with a 0.0001 s step.

Every phase keeps pristine inputs in 0.orig. Its Allrun restores those inputs and replaces its generated time directories and postprocessing, making repeated runs deterministic. It does not alter 0.orig. Run individual phases with ./steady/Allrun and ./particles/Allrun from either comparison directory. Rerun python3 check_results.py from the main tutorial directory to rebuild comparison CSVs, figures, and the acceptance report.

## Geometry and boundary conditions

The gas is 20 mm long in x, 20 mm high in y, and 1 mm deep in z. Gas mesh: 20 × 100 × 1 cells. Symmetry planes at front/back give a planar flow. A fully developed laminar inlet at x=20 mm drives dry nitrogen right to left at mean speed 0.01 m/s; the outlet at x=0 has pressure 101325 Pa. Both electrode surfaces are no slip.

In withDielectric, a 1 mm solid layer occupies y=0–1 mm, with relative permittivity 4 and a 20 × 5 × 1 mesh. Gas occupies y=1–21 mm. The lower electrode at y=0 is 0 V; the upper electrode at y=21 mm is 1000 V. Native mapped-wall patches couple potential and dielectric displacement at the interface. The solid electric field is solved, not imposed.

In noDielectric, gas occupies y=0–20 mm between 0 V and 1000 V electrodes. The gas gap is unchanged; removing the solid shortens the total electrode separation by 1 mm.

The steady run starts from zero velocity and zero potential internally. The inlet specifies the developed parabolic velocity profile; the interior velocity, pressure, and electrostatic potential are solved. Isothermal T=320 K and gravity zero isolate electric drift from buoyancy, thermal diffusion, and evaporation.

## Actual gas properties

constant/fluid/thermophysicalProperties selects sutherland, janaf, and perfectGas. foam.inp activates N2 and H2O; the initial mass fractions are N2=1 and H2O=0. This is dry nitrogen. The local foam.dat contains literature-linked N2/O2/H2O data and the supplied DEG entry, but unused O2/DEG entries do not add carrier species. See [transport properties and sources](../../docs/transportProperties.md).

At 320 K, the N2 coefficients used in the case give

~~~text
As = 1.4525113378005e-6
Ts = 123.027512608955 K
mu = As*T^(3/2)/(T+Ts) = 1.87678187747934e-05 Pa s
~~~

Density follows the ideal-gas equation of state; Cp follows the NASA/JANAF polynomial. The particle's material density and heat capacity remain simple fixed NaCl surrogate values because particle heat/mass transfer is disabled.

Every fluid case explicitly contains constant/fluid/speciesDiffusionProperties. It selects the independent dilute H2O-in-N2 molecular law, with Dref=2.505362314376117e-5 m²/s at Tref=298.15 K and pref=101325 Pa, temperature exponent 2.072 and pressure exponent 1. Its Marrero–Mason fit covers 282–373 K. Species transport is disabled here: no species gradient is needed to demonstrate the electric deposition mechanism.

## Particle and analytical comparison

One particle starts near the centre at x=10.13 mm, z=0.513 mm, and 10 mm above the collecting gas-side surface. Its diameter is 10 nm and charge number +1 means q=1 elementary charge. The force list contains nanoStokesDrag (fixed mean free path lambda=65 nm) and CoulombSep. Brownian forcing, gravity, evaporation, and parcel–carrier momentum coupling are disabled. All collecting walls use native stick interaction.

The steady field is uniform in each material, with

~~~text
Egas = 1000 / (Hgas + Hsolid/epsilonR)
Esolid = Egas/epsilonR
~~~

Both point downward, because potential increases upward. Use the case's actual viscosity in the terminal-drift theory:

~~~text
Cc = 1 + (2*lambda/d)*(1.257 + 0.4*exp(-1.1*d/(2*lambda)))
vElectric = q*abs(Egas)*Cc/(3*pi*mu*d)
tDeposit = 0.010/vElectric
~~~

| Quantity | Bare electrode | 1 mm dielectric, epsilonR=4 |
|---|---:|---:|
| Gas field magnitude (V/m) | 50000 | 49382.7160494 |
| Solid field magnitude (V/m) | — | 12345.6790123 |
| Slip correction Cc | 22.1191012568 | 22.1191012568 |
| Downward terminal speed (m/s) | 0.100175829259 | 0.098939090626 |
| Deposition time (s) | 0.099824479357 | 0.101072285349 |

The injected velocity is (-0.015, -vElectric, 0) m/s. The vertical component starts at terminal drift, avoiding a startup position offset. The horizontal component is the carrier centreline speed. As the particle descends, its horizontal speed follows the solved parabolic gas flow. The plotted horizontal theory integrates the parabolic velocity along the constant vertical-drift path. The charge and diameter are retained throughout tracking.

## Inspect the results

- steady/2000/fluid/ and, for the coated electrode, steady/2000/dielectric/: solved U, p, T, phiE, E and other native fields.
- particles/0/: copied steady state; steady_copy.json records the source time and copied fields.
- particles/postProcessing/**/trajectory.csv and deposition.csv: native cloud callback records, including event time, coordinates, diameter, charge number, and deposition patch.
- results/*_fieldProfiles.csv: solved and analytical potential/field versus height.
- results/*_trajectoryComparison.csv: native trajectory alongside analytical x/y.
- results/electricFieldStrength.png: side-by-side maps of the actual native gas/solid |E| on their physical meshes, with a shared 0–50 kV/m colour scale, native particle paths, and deposition events.
- results/*_comparison.png: potential profiles, particle path, and solved carrier flow.
- results/acceptance.json: actual deposition times, field/trajectory errors, frozen-field differences, and inlet/outlet mass flux balance.

Verification requires gas/solid electric field and deposition-time errors below 0.2%, trajectory errors below 0.2% of the 10 mm travel distance, carrier velocity L2 error below 2%, inlet/outlet mass imbalance below 1e-6, and unchanged solved carrier fields during particle tracking. Density reconstructed by heRhoThermo at restart may differ by less than 1e-8 relative; this is checked separately. The verification calculates drift from the actual copied gas field and the local N2 Sutherland coefficients. These checks use native output rather than generating a synthetic trajectory.

## Physical scope

This is an electrostatic deposition example with a solved carrier flow. It does not demonstrate space-charge feedback or a self-consistent EHD flow driven by ions: one parcel is tracked in frozen fields, particle momentum coupling is off, and no space charge is inserted into the Poisson equation. Surface charge accumulation, dielectric charging, image force, breakdown, Brownian diffusion, and realistic aerosol statistics are outside this example. The comparison isolates the voltage division introduced by the dielectric layer.

The electric-field snapshot is generated automatically by Allrun/check_results.py. After a completed run, python3 plotElectricField.py regenerates it independently from native C/E/U fields and cloud CSV files. Each image pixel is the magnitude of a solved cell field; no analytical flat-fill field is substituted.
