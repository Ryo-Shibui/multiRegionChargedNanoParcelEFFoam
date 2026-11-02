    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    UEqn =
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );

    // SU(U) was on the RHS in reactingNanoParcelFoam.  fvMatrix stores the
    // equation in residual form, so subtracting the returned matrix adds the
    // same RHS contribution without ever dereferencing a disabled cloud.
    if (solveParcels)
    {
        UEqn -= parcelsPtr->SU(U);
    }

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (momentumPredictor)
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            ),
            U.select(finalIter)
        );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
