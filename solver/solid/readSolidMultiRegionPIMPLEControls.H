const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

int nNonOrthCorr =
    pimple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);

solveEnergySolid[i] = pimple.getOrDefault<bool>("solveEnergy", true);
