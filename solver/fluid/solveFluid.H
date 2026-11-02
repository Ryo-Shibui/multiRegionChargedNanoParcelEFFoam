if (finalIter)
{
    mesh.data().setFinalIteration(true);
}

if (frozenFlow)
{
    // Fixed carrier velocity/pressure mode. Retain the original
    // reactingNanoParcelFoam capability to advance species and energy.
    if (mesh.transient() && oCorr == 0 && solveSpecies)
    {
        #include "rhoEqn.H"
    }

    if (solveSpecies)
    {
        #include "YEqn.H"
    }

    if (solveEnergy)
    {
        if (!solveSpecies)
        {
            reaction.correct();
            Qdot = reaction.Qdot();
        }
        #include "EEqn.H"
    }

    if (!coupled)
    {
        rho = thermo.rho();
    }
}
else
{
    // rhoEqn is a transient density equation. With steadyState ddt it has no
    // time diagonal and must not be solved. The steady pressure equation
    // handles mass conservation instead.
    if (mesh.transient() && oCorr == 0)
    {
        #include "rhoEqn.H"
    }

    #include "UEqn.H"

    if (solveSpecies)
    {
        #include "YEqn.H"
    }

    if (solveEnergy)
    {
        if (!solveSpecies)
        {
            reaction.correct();
            Qdot = reaction.Qdot();
        }
        #include "EEqn.H"
    }

    if (!coupled)
    {
        Info<< "\nSolving for fluid region " << fluidRegions[i].name() << endl;

        // --- Pressure-velocity correction
        for (int corr=0; corr<nCorr; corr++)
        {
            if (mesh.steady())
            {
                #include "pEqnSteady.H"
            }
            else
            {
                #include "pEqn.H"
            }
        }

        turbulence.correct();
        rho = thermo.rho();
    }
}

if (coupled && !solveEnergy)
{
    FatalErrorInFunction
        << "Energy-coupled CHT requires solveEnergy true in fluid region "
        << fluidRegions[i].name() << exit(FatalError);
}

if (finalIter)
{
    mesh.data().setFinalIteration(false);
}

if (!solveEnergy && solveSpecies)
{
    volScalarField frozenT
    (
        IOobject("frozenT",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false),
        thermo.T()
    );
    thermo.he() = thermo.he(thermo.p(), frozenT);
    thermo.correct();
    thermo.T() == frozenT;
}
