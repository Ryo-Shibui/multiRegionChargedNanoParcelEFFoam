/*---------------------------------------------------------------------------*\
  Transient continuity equation with optional reacting-parcel source.
  This file must not be solved for steadyState ddt; solveFluid/pEqn guard it.
\*---------------------------------------------------------------------------*/

{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phi)
     ==
        surfaceFilm.Srho()
      + fvOptions(rho)
    );

    if (solveParcels)
    {
        rhoEqn -= parcelsPtr->Srho(rho);
    }

    rhoEqn.solve();
    fvOptions.correct(rho);
}
