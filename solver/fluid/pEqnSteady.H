/*---------------------------------------------------------------------------*\
  Steady pressure correction adapted from OpenFOAM-v2406
  chtMultiRegionSimpleFoam/fluid/pEqn.H.

  solveParcels=true is rejected at startup for steadyState operation, so this
  equation deliberately contains no Lagrangian mass source.
\*---------------------------------------------------------------------------*/

{
    // Synchronize density with the current thermo state
    rho = thermo.rho();
    rho = max(rho, rhoMin);
    rho = min(rho, rhoMax);
    rho.relax();
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phig
    (
        -rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(rho*HbyA)
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    const bool closedVolume = adjustPhi(phiHbyA, U, p_rgh);

    phiHbyA += phig;

    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

    const dimensionedScalar compressibility = fvc::domainIntegrate(psi);
    const bool compressible = (compressibility.value() > SMALL);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rhorAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference
        (
            pRefCell,
            getRefCellValue(p_rgh, pRefCell)
        );

        p_rghEqn.solve(p_rgh.select(finalIter));

        if (nonOrth == nNonOrthCorr)
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA
              + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
            K = 0.5*magSqr(U);
        }
    }

    p = p_rgh + rho*gh;
    pressureControl.limit(p);

    if (closedVolume)
    {
        if (!compressible)
        {
            p += dimensionedScalar
            (
                "p",
                p.dimensions(),
                pRefValue - getRefCellValue(p, pRefCell)
            );
        }
        else
        {
            p += (targetMass - fvc::domainIntegrate(psi*p))
                /compressibility;
        }

        p_rgh = p - rho*gh;
        p_rgh.correctBoundaryConditions();
    }

    rho = thermo.rho();
    rho = max(rho, rhoMin);
    rho = min(rho, rhoMax);
    rho.relax();

    Info<< "Min/max rho:" << min(rho).value() << ' '
        << max(rho).value() << endl;
}
