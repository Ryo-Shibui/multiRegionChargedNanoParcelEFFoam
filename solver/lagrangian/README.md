# Reacting charged multiphase cloud libraries

These are isolated copies of the reactingChargedNanoParcelFoam cloud sources.
Only basicReactingMultiphaseCloud is instantiated; its native reacting,
multiphase, heat transfer, phase change, collision, injection, patch interaction
and cloud function model selections are retained. Standard liblagrangian supplies
basic particle infrastructure; the custom charged KinematicParcel remains local.

Build with the v2406 Opt environment sourced, then run ./Allwmake.
Outputs: libmrChargedEFDistribution.so, libmrChargedEFIntermediate.so,
libmrChargedEFTurbulence.so, all in FOAM_USER_LIBBIN. No standard/custom existing
library is overwritten.

## Charge and drag

CoulombSep uses q = chargeNumber * 1.602176634e-19 C. patchInjection and
manualInjection accept parcelCharge in elementary-charge units (default zero).
chargeNumber uses the reference parcel field IO and ASCII/binary/processor
serialization. The original sphereDrag is the native finite-Reynolds drag with
no slip correction. Select nanoStokesDrag INSTEAD of sphereDrag for low-Re
nanoparticles; selecting both adds both forces.

```
particleForces
{
    nanoStokesDrag { meanFreePath 6.6e-8; }
    CoulombSep { E E; }
}
```

nanoStokesDrag requires an explicit scalar meanFreePath or lambda in metres;
zero gives continuum Stokes drag. Exactly one key is allowed. Its coefficient is
18*mu*mass/(rho*d*d*Cc), where Kn=2*lambda/d and
Cc=1+Kn*(1.257+0.4*exp(-1.1/Kn)). This is low-Re Stokes drag, not a finite-Re
correction. The specified lambda is constant, not inferred from temperature or
pressure. Use solution/integrationSchemes/U analytical for the stiff velocity
relaxation of nanometre particles.

## Diagnostics

Add to the top-level cloudFunctions dictionary:

```
chargedDiagnostics { type chargedParticleDiagnostics; }
```

CSV files are appended under
postProcessing/lagrangian/<cloud>/<model>/<region>/processor<N>/.
Header for trajectory.csv and deposition.csv:

```
time,origProc,origId,x,y,z,d,chargeNumber,patch
```

trajectory.csv records each native parcel move substep, including inactive
parcels. deposition.csv records the first retained inactive wall hit for each
particle identity in the run; the position, diameter and charge are captured
at the hit and the time is time.value()-deltaT+stepFraction*deltaT. The native
postPatch hook runs before patch interaction, so deposition is confirmed only
when the retained parcel is subsequently inactive. postEvolve also confirms
hits exactly at stepFraction=1. Outlet/processor patches, escaped parcels and
rebound events do not count as deposition. Use native localInteraction stick
on depositing walls. Other cloud functions that freeze a parcel after a wall
hit can also satisfy the inactive check; this diagnostic targets the native
localInteraction workflow.

Ranks use distinct CSV streams and particle origProc/origId identities are
preserved through transfers. Files flush at each cloud evolution, independent
of field output intervals. On restart, rows after the checkpoint time are removed and the earlier
deposited identities are restored before appending. A fresh t=0 run clears
previous rows. Existing stuck
particles loaded on restart do not create new deposition events. Diagnostics
assume global transient time stepping; they do not add a Lagrangian subcycle or
change the native position-before-velocity update. Resolve early particle
acceleration and wall-hit timing through a time-step convergence study.

The private distribution class and implementation namespace are prefixed
mrChargedEF to avoid duplicate runtime registrations when native surface-film
libraries also load standard libdistributionModels. Dictionary type choices
and public cloud types are unchanged.
