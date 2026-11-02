{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence.alphaEff(), he)
     ==
        rho*(U&g)
      + rad.Sh(thermo, he)
      + Qdot
      + fvOptions(rho, he)
    );

    // surfaceFilm.Sh() returns a DimensionedField in OpenFOAM-v2406.
    // Add it to the assembled fvMatrix instead of directly to rho*(U&g),
    // otherwise the expression becomes volScalarField + DimensionedField.
    EEqn -= surfaceFilm.Sh();

    if (solveParcels)
    {
        EEqn -= parcelsPtr->Sh(he);
    }

    EEqn.relax();

    fvOptions.constrain(EEqn);

    if (coupled)
    {
        fvMatrixAssemblyPtr->addFvMatrix(EEqn);
    }
    else
    {
        EEqn.solve(he.select(finalIter));
        fvOptions.correct(he);

        thermo.correct();
        rad.correct();

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;
    }
}
