forAll(fluidRegions, i)
{
    if (!solveEnergyFluid[i]) continue;
    rhoThermo& thermo = thermoFluid[i];
    radiation::radiationModel& rad = radiation[i];
    fv::options& fvOptions = fluidFvOptions[i];
    volScalarField& he = thermo.he();
    fvOptions.correct(he);
    thermo.correct();
    rhoFluid[i] = thermo.rho();
    rad.correct();
}

forAll(solidRegions, i)
{
    if (!solveEnergySolid[i]) continue;
    solidThermo& thermo = thermos[i];
    fv::options& fvOptions = solidHeatSources[i];
    volScalarField& h = thermo.he();
    fvOptions.correct(h);
    thermo.correct();
}

