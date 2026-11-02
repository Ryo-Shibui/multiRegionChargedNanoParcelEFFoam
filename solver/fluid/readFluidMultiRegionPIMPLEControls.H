    const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

    const int nCorr =
        pimple.getOrDefault<int>("nCorrectors", 1);

    const int nNonOrthCorr =
        pimple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const bool momentumPredictor =
        pimple.getOrDefault("momentumPredictor", true);

    // Re-read the Eulerian solve controls at run time.
    // solveFlow takes precedence over the legacy frozenFlow entry.
    bool solveFlow = !frozenFlowFluid[i];
    if (pimple.readIfPresent("solveFlow", solveFlow))
    {
        frozenFlowFluid[i] = !solveFlow;
    }
    else
    {
        pimple.readIfPresent("frozenFlow", frozenFlowFluid[i]);
    }

    pimple.readIfPresent("solveSpecies", solveSpeciesFluid[i]);
    pimple.readIfPresent("solveEnergy", solveEnergyFluid[i]);

    // solveParcels is intentionally startup-only because false->true would
    // require constructing SLGThermo and the cloud after field creation.
    const bool requestedSolveParcels =
        pimple.getOrDefault("solveParcels", solveParcelsFluid[i]);

    if (requestedSolveParcels != solveParcelsFluid[i])
    {
        WarningInFunction
            << "Ignoring run-time change of solveParcels for region "
            << mesh.name() << ". solveParcels is read at solver startup; "
            << "restart the solver to change it." << endl;
    }
