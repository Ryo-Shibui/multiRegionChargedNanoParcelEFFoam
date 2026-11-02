
{
    bool anyElectric=false;
    forAll(electricMeshes,i)
    {
      solveElectric[i]=electricMeshes[i]->solutionDict().subDict("PIMPLE").getOrDefault<bool>("solveElectrostatics",true);
      anyElectric=anyElectric||solveElectric[i];
    }
    if(anyElectric)
    {
      fvSolution globalSolution(runTime);
      const dictionary& ec=globalSolution.subDict("PIMPLE");
      const label maxIt=ec.getOrDefault<label>("maxElectricIterations",20000);
      const scalar tol=ec.getOrDefault<scalar>("electricTolerance",1e-10);
      bool electricConverged=false;
      for(label iteration=0;iteration<maxIt;++iteration)
      {
        scalar change=0;
        forAll(electricMeshes,i)if(solveElectric[i])
        {
          auto& potential=phiE[i];
          const scalarField previous(potential.primitiveField());
          fvScalarMatrix electricEqn(fvm::laplacian(epsilonR[i],potential));
          electricEqn.solve();
          change=max(change,gMax(mag(potential.primitiveField()-previous))/max(gMax(mag(potential.primitiveField())),scalar(1)));
        }
        forAll(electricMeshes,i)if(solveElectric[i])phiE[i].correctBoundaryConditions();
        if(change<tol){electricConverged=true;break;}
      }
      if(!electricConverged)FatalErrorInFunction<<"Electrostatic region coupling did not converge"<<exit(FatalError);
      forAll(electricMeshes,i)if(solveElectric[i])
      {
        electricField[i]=-fvc::grad(phiE[i]);
        electricField[i].correctBoundaryConditions();
      }
    }
}
