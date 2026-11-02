# multiRegionChargedNanoParcelEFFoam

Multi-region compressible reacting flow and solid CHT with charged reacting nano parcels and a dielectric potential solve in fluid/solid regions. Choose this solver when electrode potentials and dielectric interfaces must determine E. Region meshes are static.

## Choosing between the solvers

| Solver | Mesh/thermal scope | Electric field |
|---|---|---|
| reactingChargedNanoParcelFoam | One fluid region; dynamic-mesh framework | Prescribed E |
| chtReactingNanoParcelFoam | Static fluid/solid regions, CHT | Prescribed E for CoulombSep |
| multiRegionChargedNanoParcelEFFoam | Static fluid/solid regions, CHT | Dielectric potential solve, or frozen E |

These retain reacting thermophysical/cloud infrastructure. Skipping equations does not remove startup thermo fields or create a minimal isothermal kinematic solver.

The release adds independent dilute-species gas diffusivities and explicit dilute/legacy selection, plus nanoStokesDrag and charge IO/injection where missing. Existing thermal/reaction/film functionality is retained. The sibling reactingNanoParcelFoam solver and its shared libraries are not modified or required to be rebuilt.

## Quick start with the supplied cases

After building, copy a case outside `solver/` and run it:

~~~sh
mkdir -p "$FOAM_RUN"
cp -a tutorials/diluteSpecies "$FOAM_RUN/multiRegionChargedNanoParcelEFFoam-dilute"
multiRegionChargedNanoParcelEFFoam -case "$FOAM_RUN/multiRegionChargedNanoParcelEFFoam-dilute"
python3 tests/check_cases.py "$FOAM_RUN/multiRegionChargedNanoParcelEFFoam-dilute"
# Or copy tutorials/chargedDrift10nm and run the same executable.
# The tutorials include a small ready-to-run mesh.
bash tests/Allrun
~~~

Validation scripts require Python 3 and NumPy. Tutorial gas properties deliberately
use constant density/viscosity to isolate transport verification; replace them with
validated thermophysical data for a physical DEG application.

## Repository layout and build

~~~text
README.md, LICENSE
Allwmake, Allwclean       repository-level dispatchers
solver/                  full build sources and native build scripts
  include/diluteSpeciesDiffusion.H
  lagrangian/            bundled custom cloud sources
tutorials/               self-contained native cases, outside compiled sources
tests/                   case generation, numerical checks and release results
~~~

Fluid/solid equations are in solver/fluid and solver/solid.

The electrostatic solve/interface condition is in solver/electrostatics.

Use the **OpenCFD/OpenFOAM.com** family; Foundation numbered releases are a different code line. Intended versions are listed below with actual validation status.

A development installation needs source headers, wmake, runtime libraries, a compatible C++ compiler, make and the installation's MPI/runtime dependencies. Python 3 supports validation scripts, not the solver runtime.

~~~sh
source /path/to/OpenFOAM-v2406/etc/bashrc
cd /path/to/multiRegionChargedNanoParcelEFFoam
./Allwmake
command -v multiRegionChargedNanoParcelEFFoam
multiRegionChargedNanoParcelEFFoam -help
~~~

The root dispatcher builds solver/Allwmake. Native build scripts compile the bundled private Lagrangian libraries before the executable and build mrChargedEFElectrostatics first. Outputs use FOAM_USER_APPBIN/LIBBIN; private library names start with mrChargedEF. Build scripts use the installed native lagrangian core; do not replace its library with a bundled custom basic library or rename private outputs to common shared-library names.

Build and run under the same sourced environment. Use separate source copies and bin/lib outputs per release: identical WM_OPTIONS does not mean compatible ABI. Do not carry old objects, generated lnInclude or binaries into another release.

~~~sh
./Allwclean
./Allwmake
multiRegionChargedNanoParcelEFFoam -case /path/to/working-case
# After configuring and running the appropriate mesh decomposition:
mpirun -np 2 multiRegionChargedNanoParcelEFFoam -parallel -case /path/to/working-case
~~~

The parallel command is an invocation example, not evidence for every case/model. Copy tutorials into a separate working directory before meshing/running; generated fields and logs do not belong in solver/.

## Native case structure

~~~text
case/
  0/
    fluid/       U, p, p_rgh, T, carrier species Yi
    solid/       solid thermo/temperature fields
  constant/
    regionProperties, g
    fluid/
      polyMesh/
      thermophysicalProperties
      chemistryProperties, combustionProperties
      turbulenceProperties, radiationProperties, surfaceFilmProperties
      speciesDiffusionProperties
      reactingCloud1Properties, particlePositions  [when parcels enabled]
    solid/
      polyMesh/
      thermophysicalProperties
  system/
    controlDict, fvSolution
    fluid/fvSchemes, fluid/fvSolution
    solid/fvSchemes, solid/fvSolution
~~~

Region names are illustrative. Set constant/regionProperties, for example:

~~~foam
regions
(
    fluid (fluid)
    solid (solid)
);
~~~

Use valid reacting fluid thermo, solid thermo and thermal interface conditions. Multiple fluids/solids, or no solids, can be declared. Each region owns a registry: its cloud/positions/diffusion dictionary belongs under `constant/<fluid>`, and its field under `<time>/<fluid>`. Do not put region inputs at the case root.

Electric runs additionally need phiE and electricProperties for every enabled electric region; frozen-field clouds need E. See below.

For inert carrier tests, disable chemistry/combustion with their native model selections; flow freezing does not disable species reaction sources. When no film is needed, the Eulerian surface-film framework still constructs independently of parcel activity. Supply:

~~~foam
FoamFile { format ascii; class dictionary; object surfaceFilmProperties; }
surfaceFilmModel none;
~~~

## Split physics controls

In `system/<fluid>/fvSolution/PIMPLE`:

~~~foam
PIMPLE
{
    solveFlow false;
    solveSpecies true;
    solveEnergy false;
    solveParcels false;
    nCorrectors 1;
    nNonOrthogonalCorrectors 0;
}
~~~

Defaults: flow/species/energy/parcels enabled. solveFlow overrides the legacy frozenFlow flag (default false). Flow/species/energy flags are reread; solveParcels is startup-only, and changing it requires restarting. Frozen flow may still solve species and/or energy; fluid energy is independent of the species flag.

solveParcels false avoids SLGThermo/cloud construction and its cloud/liquid/solid parcel inputs. Carrier thermo, radiation and Eulerian film remain. Enabled clouds evolve once per physical step before outer thermal/fluid iterations; source terms are reused within that step.

A steadyState ddt scheme requires solveParcels false. Obtain a steady carrier/CHT solution, then restart with transient schemes and solveParcels true for physical flight/residence times. Steady Time is an iteration label. Active regions need matching steady/transient formulations.

Global system/fvSolution supplies PIMPLE/nOuterCorrectors. Implicit assembled thermal coupling requires solveEnergy true for participating regions. Freezing energy requires compatible explicit thermal interface treatment; it cannot leave an active assembled energy interface in a frozen region.

The solid region also accepts PIMPLE/solveEnergy (default true), independent of fluid energy.

Equation skipping does not remove startup T/thermo requirements, and is not equivalent to disabling parcel phase-change/heat-transfer submodels.

## Dielectric electric solve

The solver computes div(epsilonR grad(phiE))=0 and E=-grad(phiE), with no volume free-charge term. epsilon0 cancels in the charge-free equation. In each enabled electric region, provide `constant/<region>/electricProperties`:

~~~foam
FoamFile { format ascii; class dictionary; object electricProperties; }
epsilonR 4;
~~~

Use positive epsilonR (for example 1 for gas). phiE has dimensions [1 2 -3 0 0 -1 0]; E has [1 1 -3 0 0 -1 0]. Electrode phiE patches use fixedValue. At fluid/solid interfaces use electricPotentialCoupled on both sides, with value uniform 0 as the initial boundary value. The mesh patches must be mapped to the correct neighbor region/patch. The field condition alone does not create mapping.

The interface uses epsilonR and each side's face/cell distance to couple potential and normal displacement. Normal E is generally discontinuous for unequal permittivities. Supply a phiE linear solver in each region fvSolution.

Per-region PIMPLE/solveElectrostatics defaults true. Global system/fvSolution/PIMPLE entries are maxElectricIterations (default 20000) and electricTolerance (default 1e-10). The tolerance monitors normalized internal-potential change; it is not an electric-field/flux accuracy criterion.

The electric solve runs before tracking and at physical time levels. Coupling iterations do not advance particle time. With solveElectrostatics false, supplied E is kept; enabled frozen-field fluid clouds require E at the initial/restart time. If all electric solves are off, phiE/electricProperties may be absent. A two-stage run can converge/write phiE/E, then restart with the field frozen and chosen fluid flags; measure flight time relative to injection/restart, not an iteration label.

No ion/charge transport, corona, particle space-charge or deposited surface-charge feedback is solved.


## Independent dilute-species transport

Species fields Yi are **mass fractions**. Place speciesDiffusionProperties in the carrier's constant registry: constant/speciesDiffusionProperties for the single-region solver; `constant/<fluid>/speciesDiffusionProperties` for a multi-region solver. A complete explicit example is:

~~~foam
FoamFile { format ascii; class dictionary; object speciesDiffusionProperties; }
model dilute;
turbulentSchmidt 0.9;
species
{
    H2O
    {
        D [0 2 -1 0 0 0 0] 2.5e-5;
        Tref 298.15;
        pref 101325;
        temperatureExponent 1.75;
        pressureExponent 1;
    }
    DEG
    {
        D [0 2 -1 0 0 0 0] 8e-6;
        Tref 298.15;
        pref 101325;
        temperatureExponent 1.75;
        pressureExponent 1;
    }
}
~~~

The numbers are illustrative, not calibrated material data. H2O and DEG must be actual carrier gas species. Specify the dominant background in thermophysicalProperties, for example inertSpecie N2;. Every other **active** species requires its own D entry; no N2 entry is required.

~~~text
Di(T,p) = Dref,i (T/Tref,i)^temperatureExponent,i
                  (pref,i/p)^pressureExponent,i
Gammai = rho Di + mut / turbulentSchmidt       [mut = muEff - mu]

d(rho Yi)/dt + div(phi Yi) - div(Gammai grad(Yi)) = Si
Ybackground = max(0, 1 - sum(Yi for solved active non-background species))
~~~

phi is mass flux. Si includes enabled parcel, reaction, surface-film and fvOptions sources. D is a molecular diffusivity in m2/s; rho D and mut have dynamic-viscosity dimensions. Molecular species diffusivity is independent of momentum diffusivity.

| Entry | Units / type | Default or requirement |
|---|---|---|
| model | dilute or legacy | Absent/empty dictionary: legacy; nonempty dictionary without model: dilute |
| turbulentSchmidt | scalar, dimensionless | 0.9; positive, finite |
| species.<name>.D | dimensioned scalar, m2/s | Mandatory for each active non-background specie; positive, finite |
| species.<name>.Tref | plain scalar, K | 298.15; positive, finite |
| species.<name>.pref | plain scalar, Pa absolute | 101325; positive, finite |
| temperatureExponent | plain scalar | 0; finite |
| pressureExponent | plain scalar | 0; finite |

D may also be a bare SI scalar where the dimensionedScalar reader permits it. Tref/pref are plain-scalar entries: do not attach dimension sets. Zero exponents give constant D. Nonzero exponents require positive local absolute T/p. Restart after editing coefficients; they are read at startup, even when the species solve is disabled.

An absent dictionary or explicit model legacy; preserves the old Gamma=muEff behavior, effectively molecular/turbulent Schmidt numbers of one, and prints a legacy notice. It does not give independent gas diffusivities. Invalid explicit dilute input (missing species/D, invalid D or model) fails instead of silently reverting.

### Approximation and composition limits

This is independent dilute mass-fraction Fick transport against a dominant background, not Maxwell-Stefan multicomponent diffusion. Background closure implies the balancing background flux; no separate background PDE is solved. Cross-diffusion, Soret, barodiffusion and Eulerian species electromigration are absent. The pressure exponent changes D only.

Use trace concentrations with total dilute mass fraction much smaller than one. All non-background species contributing to the mixture should be active: closure sums the solved active species. Nonzero inactive species and sums above one are outside this dilute setup. Clipping is a numerical bound, not a valid replacement for conservative composition.

For ideal-gas partial pressure/saturation calculations, convert mass to mole fractions:

~~~text
xi = (Yi/Wi) / sum_j(Yj/Wj)
pi = xi p
Si = pi / psat,i(T)
~~~

Do not directly substitute mass fraction Yi for mole fraction xi.

## Charged nano parcels

The cloud is **reactingCloud1**, not reactCloud1. The native dictionary is reactingCloud1Properties. Each parcel carries a signed chargeNumber:

~~~text
q = chargeNumber * 1.602176634e-19 C
FCoulomb = q E
~~~

This is charge per physical particle, not per computational parcel; nParticle sets the number of represented particles. Positive charge drifts along E; negative charge against E. The runtime force is case-sensitive CoulombSep. Its E entry chooses a named electric field, default E.

chargeNumber is written with Lagrangian fields and read on restart. An old field-based cloud without it defaults to zero charge. Keep this file for charged restarts; binary parcels are not guaranteed interchangeable with stock types or another label/scalar build.

### Drag selection

~~~foam
particleForces
{
    nanoStokesDrag { meanFreePath 6.5e-8; }
    CoulombSep { E E; }
}
~~~

nanoStokesDrag accepts meanFreePath **or** lambda, one plain SI scalar in metres. Supply exactly one; negative values fail. The configured mean free path is constant, not automatically scaled with T/p.

~~~text
Kn = 2 lambda/dp
Cc = 1 + Kn [1.257 + 0.4 exp(-1.1/Kn)]       (Cc=1 for lambda=0)
Fdrag = (3 pi mu dp/Cc) (U-Up)
tauP = rhoP dp^2 Cc/(18 mu)
terminal electric drift = q Cc E/(3 pi mu dp)
~~~

This is a low-Re spherical-particle model. sphereDrag remains a continuum sphere correlation with **no Cunningham correction** in this source. Choose one drag model: enabling both adds their drags. Brownian motion, gravity, turbulent dispersion, collisions, image force and dielectrophoresis are not automatically included by either model.

Use analytical velocity integration for stiff nanoparticle drag, then verify position/contact-time convergence. Stability alone is not trajectory accuracy. Choose native patchInteractionModel behavior (stick/escape/rebound) explicitly; the tracked particle is a point unless a separate wall model says otherwise.

### Native cloud dictionary: inert 10 nm solid particle

This example uses no reactions/phase change and one physical particle per parcel. It still uses the reacting-multiphase composition infrastructure. Provide matching NaCl solid properties in SLGThermo input, or replace NaCl with a material for which you provide consistent data. Patch names below must match your mesh.

~~~foam
FoamFile { format ascii; class dictionary; object reactingCloud1Properties; }
    solution
    {
        active          true;
        coupled         no;
        transient       yes;
        cellValueSourceCorrection no;
        maxCo           0.3;

        sourceTerms
        {
            schemes
            {
                rho             explicit 1;
                U               explicit 1;
                Yi              explicit 1;
                h               explicit 1;
                radiation       explicit 1;
            }
        }

        interpolationSchemes
        {
            rho             cell;
            U               cellPoint;
            E               cell;
            thermo:mu       cell;
            mu              cell;
            T               cell;
            Cp              cell;
            kappa           cell;
            p               cell;
        }

        integrationSchemes
        {
            U               analytical;
            T               analytical;
        }
    }

    constantProperties
    {
        parcelTypeId        -1;
        rhoMin              1e-15;
        rho0                1000;
        T0                  298.15;
        TMin                200;
        TMax                5000;
        Cp0                 850;
        epsilon0            0;
        f0                  0;
        pMin                1000;
        minParcelMass       1e-30;
        volumeUpdateMethod  constantRho;
        TDevol              1e30;
        LDevol              0;
        hRetentionCoeff     0;
        lockParcelDiameterSet false;
        lockParcelDiameter  0;
    }

    subModels
    {
        particleForces
        {
            nanoStokesDrag { meanFreePath 6.5e-8; }
            CoulombSep
            {
                E E;
            }
        }

        injectionModels
        {
            centreParticle
            {
                type manualInjection;
                positionsFile particlePositions;
                SOI 0;
                parcelBasisType fixed;
                nParticle 1;
                massTotal 0;
                U0 (0 0 0);
                parcelCharge 1;
                sizeDistribution { type fixedValue; fixedValueDistribution { value 1e-8; } }
            }
        }

        dispersionModel none;
        patchInteractionModel localInteraction;
        chemicalReactionModel none;
        collisionModel none;
        heatTransferModel none;
        compositionModel singleMixtureFraction;
        phaseChangeModel none;
        devolatilisationModel none;
        surfaceReactionModel none;
        stochasticCollisionModel none;
        radiationModel none;
        surfaceFilmModel none;
        radiation off;

        localInteractionCoeffs
        {
            patches
            (
                inlet { type escape; }
                outlet { type escape; }
                wall { type stick; }
            );
        }

        singleMixtureFractionCoeffs
        {
            phases
            (
                gas {}
                liquid {}
                solid
                {
                    NaCl 1;
                }
            );
            YGasTot0        0;
            YLiquidTot0     0;
            YSolidTot0      1;
        }
    }


cloudFunctions {}

~~~

The manual positions file belongs to the same constant registry as the cloud dictionary (constant/particlePositions or `constant/<fluid>/particlePositions`):

~~~foam
FoamFile { format ascii; class vectorField; object particlePositions; }
1
(
    (0.005 0.005 0.005)
)
~~~

Place the point inside your mesh. manualInjection injects the listed positions at SOI. parcelCharge sets the signed charge **number**, not coulombs; its default is zero. sizeDistribution/value is diameter in metres. Use physical integral charge states when appropriate; scalar storage also permits an explicit effective-charge model.

For a patch source replace the manual entry with:

~~~foam
patchSource
{
    type patchInjection;
    patch inlet;
    SOI 0;
    duration 0.1;
    parcelsPerSecond 100;
    parcelBasisType fixed;
    nParticle 1;
    massTotal 0;
    velocityType fixedValue;
    U0 (0.1 0 0);
    flowRateProfile constant 1;
    parcelCharge -1;
    sizeDistribution
    {
        type fixedValue;
        fixedValueDistribution { value 1e-8; }
    }
}
~~~

Both manualInjection and patchInjection set chargeNumber from parcelCharge. Do not assume every other injector implements it.

solution/coupled false disables parcel source feedback. true enables the implemented mass/species/momentum/enthalpy source coupling. CoulombSep is a non-coupled external force: direct qE is not added again as a gas body force. This switch does not enable electric space-charge feedback.

## Phase change and DEG applications

The existing phase-change runtime model is liquidEvapFuchsKnudsen. This release does not supply calibrated DEG gas/liquid material functions.

A DEG application requires user-supplied consistent gas thermochemistry, molecular weight, liquid-property functions (density, Cp, vapour pressure, surface tension, latent heat and diffusivity), and a matching cloud composition. The phase-change source evaluates:

~~~text
Dab = liquids.properties()[liquidId].D(pc, Ts)
~~~

That liquidProperties diffusivity is **separate** from the Eulerian gas D in speciesDiffusionProperties. Changing one does not override the other. Keep their definitions consistent for the physical problem.

The native coefficient dictionary is liquidEvapFuchsKnudsenCoeffs. Existing entries include gamma, alpham, solution and activityCoefficient; Hoff additionally reads ic. This illustrates the dictionary shape, not calibrated DEG/NaCl data:

~~~foam
phaseChangeModel liquidEvapFuchsKnudsen;
liquidEvapFuchsKnudsenCoeffs
{
    gamma 6.5e-8; // mean gas free path, m
    alpham 1.0;
    solution (DEG NaCl);
    activityCoefficient Hoff;
    ic 1.0;
}
~~~

Every species/phase name must map to the carrier and cloud composition. UNIFAC is not implemented. All three solvers support one liquid, for example solution (H2O);, as well as a two-component (liquid solid) list. Pure-liquid numerical verification is recorded with the release results. No pure-liquid option supplies missing DEG material data. A cloud with phaseChangeModel none still requires valid properties for phases declared by its composition model.

Freezing carrier energy does not switch off parcel heat transfer or phase change: choose those cloud submodels separately and assess the resulting thermal assumption.

## Validation status and limits

Built and numerically verified with **OpenCFD OpenFOAM v2406, v2412 and v2506** on Ubuntu 24.04 (WSL2), GCC 13.3, double precision and 32-bit labels. The v2406 build is 630d60de3b-20240620; v2412/v2506 use official patch 260127 packages. This statement covers the tested builds and cases, not every platform or third-party model combination:

| OpenCFD release | Native build | Two-species diffusion / T,p law / invalid input | Charged drift / tracking | Status |
|---|---|---|---|---|
| v2406 | PASS | PASS | PASS | Numerically verified |
| v2412 | PASS | PASS | PASS | Numerically verified |
| v2506 | PASS | PASS | PASS | Numerically verified |

The diffusion test compares two trace species in N2 with independent analytical decays A(t)=A(0) exp(-D k^2 t). Other cases exercise temperature/pressure exponents, absent/explicit legacy settings and invalid-input rejection. They do not validate DEG material functions or a condensation experiment.

Machine-readable results are in [validation/](validation/). The release checks include independent-D transient decay, the temperature/pressure law, turbulent Schmidt mixing, legacy compatibility, invalid-input rejection, steady frozen fields, signed 10 nm drift, restart, and benign pure-water phase change. CHT/EF versions also pass a fluid-solid conduction test. CHT passes a 10000-substep advancing-trajectory regression. The EF version additionally passes dielectric/flow/deposition cases and a two-rank particle ownership/charge transfer test on each release. Dynamic-mesh operation, reactive chemistry accuracy, DEG material data, and every two-way source balance are not covered by these verification cases.

Before application:
- Check species mass-fraction sums, trace loading, species activation, absolute T/p, property validity and boundary conditions.
- Check selected forces and phase/chemistry submodels in startup output.
- Verify mesh/time-step convergence of trajectories, contact times and coupled source balances.
- Decide whether neglected Brownian motion, gravity, turbulence, collisions and additional electric forces are physically acceptable.
- The electric model does not add charging kinetics, ion/corona generation, deposited surface-charge accumulation, space-charge, image forces or dielectric polarization forces.

## License

GPL-3.0-or-later, consistent with the bundled OpenFOAM-derived sources. Retain upstream notices and license terms when redistributing modified code.

## Reproducing the transport checks

Run the supplied acceptance cases:

~~~sh
python3 tests/make_cases.py --solver multiRegionChargedNanoParcelEFFoam
bash tests/Allrun --solver multiRegionChargedNanoParcelEFFoam
~~~

Keep run cases outside solver/. The generator's independent expected values and acceptance summaries should be retained with the release report.

## Changes relevant when migrating an existing case

- Keep native `reactingCloud1Properties`, `singleMixtureFractionCoeffs`, and
  `inertSpecie`. A new `speciesDiffusionProperties` dictionary selects independent
  molecular diffusivities; absent dictionaries retain the legacy coefficient.
- The species Laplacian keeps the legacy `laplacian(muEff,<species>)` scheme
  lookup name, including in dilute mode. A `default Gauss linear orthogonal;`
  entry is sufficient for the supplied mesh.
- Transient injection defers movement to the normal cloud tracking step and
  preserves the initial partial timestep. This avoids advancing a new parcel twice.
- Pure-liquid Fuchs-Knudsen input no longer reads an absent solid component.
- Private distribution/parcel libraries keep their own output names and distribution
  namespaces. Do not load an older similarly modified cloud library in the same
  process, or add it to `controlDict/libs`.

## Included tutorial cases

| Case | Purpose |
|---|---|
| `tutorials/diluteSpecies` | Independent DEG/H2O tracer diffusion in N2; constant-property analytic verification |
| `tutorials/chargedDrift10nm` | One +1e, 10 nm particle; prescribed field and Stokes-Cunningham drift |
| `tutorials/pureLiquidH2O` | Pure-liquid Fuchs-Knudsen runtime and native parcel composition check |
| `tutorials/solidHeatConduction` | Fluid/solid temperature solution and interface heat flux |
| `tutorials/parallelPlates10nm` | Steady field/flow, then transient deposition; with/without 1 mm dielectric |

Run a copied case with the executable; run `bash tests/Allrun` for automatic
generation and numerical assertions. The reported 10 nm terminal drift is
0.1038719232 m/s at 50000 V/m. Finite-step displacement differs from the exact
relaxation solution by about 0.059%, and the terminal velocity matches to floating
point precision. Do not interpret the illustrative DEG diffusivity as measured data.

For the two-stage electrode examples, run `bash tutorials/parallelPlates10nm/Allrun`.
The 20 mm gas gap gives 50000 V/m without dielectric and 49382.7160494 V/m with
the 1 mm, epsilonR=4 dielectric. Measured deposition times are 0.0962723668 s
and 0.0974757710 s, with absolute theory differences below 3.4e-8 s.
Use `RUN_MPI=1 bash tests/Allrun` to include the two-rank ownership-transfer test.
