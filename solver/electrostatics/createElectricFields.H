
List<fvMesh*> electricMeshes(fluidRegions.size()+solidRegions.size());
forAll(fluidRegions,i) electricMeshes[i]=&fluidRegions[i];
forAll(solidRegions,i) electricMeshes[i+fluidRegions.size()]=&solidRegions[i];
PtrList<volScalarField> epsilonR(electricMeshes.size()), phiE(electricMeshes.size());
PtrList<volVectorField> electricField(electricMeshes.size());
List<bool> solveElectric(electricMeshes.size(),true);
forAll(electricMeshes,i)
{
    fvMesh& em=*electricMeshes[i];
    solveElectric[i]=em.solutionDict().subDict("PIMPLE").getOrDefault<bool>("solveElectrostatics",true);
    IOobject propIO("electricProperties",runTime.constant(),em,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE);
    const bool haveProps=propIO.typeHeaderOk<IOdictionary>(true);
    if(solveElectric[i] && !haveProps)
      FatalErrorInFunction<<"Missing electricProperties in "<<em.name()<<exit(FatalError);
    scalar er=1;
    if(haveProps){IOdictionary props(propIO);er=props.get<scalar>("epsilonR");}
    if(er<=0)FatalErrorInFunction<<"epsilonR must be positive in "<<em.name()<<exit(FatalError);
    epsilonR.set(i,new volScalarField(IOobject("epsilonR",runTime.timeName(),em,IOobject::NO_READ,IOobject::AUTO_WRITE),em,dimensionedScalar(dimless,er)));
    IOobject potentialIO("phiE",runTime.timeName(),em,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE);
    if(potentialIO.typeHeaderOk<volScalarField>(true))
    {
      potentialIO.readOpt(IOobject::MUST_READ);
      phiE.set(i,new volScalarField(potentialIO,em));
    }
    else
    {
      if(solveElectric[i])FatalErrorInFunction<<"Missing phiE in "<<em.name()<<exit(FatalError);
      phiE.set(i,new volScalarField(potentialIO,em,dimensionedScalar(dimensionSet(1,2,-3,0,0,-1,0),0)));
    }
    if(phiE[i].dimensions()!=dimensionSet(1,2,-3,0,0,-1,0))
      FatalErrorInFunction<<"phiE must have dimensions of voltage"<<exit(FatalError);
    IOobject eIO("E",runTime.timeName(),em,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE);
    if(eIO.typeHeaderOk<volVectorField>(true))
    {
      eIO.readOpt(IOobject::MUST_READ);
      electricField.set(i,new volVectorField(eIO,em));
    }
    else
    {
      if(!solveElectric[i] && i<fluidRegions.size() && solveParcelsFluid[i])
        FatalErrorInFunction<<"Frozen-field parcel tracking requires an E field in "<<em.name()<<exit(FatalError);
      electricField.set(i,new volVectorField(eIO,em,dimensionedVector(dimensionSet(1,1,-3,0,0,-1,0),Zero)));
    }
    if(electricField[i].dimensions()!=dimensionSet(1,1,-3,0,0,-1,0))
      FatalErrorInFunction<<"E must have dimensions of voltage/length"<<exit(FatalError);
}
