/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019,2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    multiRegionChargedNanoParcelEFFoam

Group
    grpHeatTransferSolvers

Description
    Transient multi-region solver for compressible reacting fluid flow,
    solid heat conduction and conjugate heat transfer, augmented with the
    reactingNanoParcelFoam Lagrangian reacting-multiphase cloud models.

    Each fluid region owns an independent reacting cloud and optional surface
    film. Parcel-to-carrier mass, species, momentum and enthalpy source terms
    are coupled to the Eulerian equations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "diluteSpeciesDiffusion.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "pressureControl.H"
#include "surfaceFilmModel.H"
#include "SLGThermo.H"
#include "cloudMacros.H"

#ifndef CLOUD_BASE_TYPE
    #define CLOUD_BASE_TYPE ReactingMultiphase
    #define CLOUD_BASE_TYPE_NAME "reacting"
#endif

#include CLOUD_INCLUDE_FILE(CLOUD_BASE_TYPE)
#define basicReactingTypeCloud CLOUD_TYPE(CLOUD_BASE_TYPE)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc,char* argv[])
{
    argList::addNote("Multi-region electrostatic CHT and charged reacting nano parcels; fvSolution split physics controls.");
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    PtrList<diluteSpeciesDiffusion> speciesDiffusionFluid(fluidRegions.size());
    forAll(fluidRegions,i)
        speciesDiffusionFluid.set(i,new diluteSpeciesDiffusion
        (fluidRegions[i],thermoFluid[i].composition(),
         thermoFluid[i].get<word>("inertSpecie")));
    #include "createElectricFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"

    bool steadyRun=true, transientRun=true;
    forAll(fluidRegions,i)
    {
      steadyRun=steadyRun&&fluidRegions[i].steady();
      transientRun=transientRun&&fluidRegions[i].transient();
      Info<<"Controls "<<fluidRegions[i].name()
          <<" solveFlow="<<!frozenFlowFluid[i]<<" solveSpecies="<<solveSpeciesFluid[i]
          <<" solveEnergy="<<solveEnergyFluid[i]<<" solveParcels="<<solveParcelsFluid[i]<<endl;
    }
    forAll(solidRegions,i)if(solveEnergySolid[i])
    {
      steadyRun=steadyRun&&solidRegions[i].steady();
      transientRun=transientRun&&solidRegions[i].transient();
    }
    if(!steadyRun&&!transientRun)
      FatalErrorInFunction<<"Active CHT regions must use matching steady/transient ddt schemes"<<exit(FatalError);
    if(transientRun)
    {
      #include "compressibleMultiRegionCourantNo.H"
      #include "solidRegionDiffusionNo.H"
      #include "surfaceFilmMultiRegionCourantNo.H"
      #include "setInitialMultiRegionDeltaT.H"
    }
    #include "createCoupledRegions.H"
    #include "solveElectricPotential.H"

    while(runTime.run())
    {
      #include "readTimeControls.H"
      #include "readSolidTimeControls.H"
      #include "readPIMPLEControls.H"
      forAll(fluidRegions,i)
      {
        fvMesh& mesh=fluidRegions[i];
        #include "readFluidMultiRegionPIMPLEControls.H"
      }
      forAll(solidRegions,i)
        solveEnergySolid[i]=solidRegions[i].solutionDict().subDict("PIMPLE").getOrDefault<bool>("solveEnergy",true);
      if(coupled)
      {
        forAll(fluidRegions,i)if(!solveEnergyFluid[i])
          FatalErrorInFunction<<"Use explicit thermal interface coupling when freezing a temperature region"<<exit(FatalError);
        forAll(solidRegions,i)if(!solveEnergySolid[i])
          FatalErrorInFunction<<"Use explicit thermal interface coupling when freezing a temperature region"<<exit(FatalError);
      }
      if(transientRun)
      {
        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "surfaceFilmMultiRegionCourantNo.H"
        #include "setMultiRegionDeltaT.H"
      }
      ++runTime;
      Info<<"Time = "<<runTime.timeName()<<nl<<endl;
      #include "solveElectricPotential.H"
      if(transientRun)forAll(fluidRegions,i)
      {
        if(solveParcelsFluid[i])parcelsFluid[i].evolve();
        surfaceFilmFluid[i].evolve();
        scalar massRate=fvc::domainIntegrate(surfaceFilmFluid[i].Srho()).value();
        if(solveParcelsFluid[i])massRate+=fvc::domainIntegrate(parcelsFluid[i].Srho()).value();
        targetMassFluid[i]+=runTime.deltaTValue()*massRate;
      }
      // The custom loop owns iteration storage, including one-corrector SIMPLE.
      forAll(fluidRegions,i)
      {
        #include "storeOldFluidFields.H"
      }
      for(int oCorr=0;oCorr<nOuterCorr;++oCorr)
      {
        const bool finalIter=transientRun&&(oCorr==nOuterCorr-1);
        forAll(fluidRegions,i)
        {
          fvMesh& mesh=fluidRegions[i];
          #include "readFluidMultiRegionPIMPLEControls.H"
          #include "setRegionFluidFields.H"
          #include "solveFluid.H"
        }
        forAll(solidRegions,i)
        {
          fvMesh& mesh=solidRegions[i];
          #include "readSolidMultiRegionPIMPLEControls.H"
          #include "setRegionSolidFields.H"
          #include "solveSolid.H"
        }
        if(coupled)
        {
          fvMatrixAssemblyPtr->solve();
          #include "correctThermos.H"
          forAll(fluidRegions,i)
          {
            fvMesh& mesh=fluidRegions[i];
            #include "readFluidMultiRegionPIMPLEControls.H"
            #include "setRegionFluidFields.H"
            if(!frozenFlow)
            {
              for(int corr=0;corr<nCorr;++corr)
              {
                if(mesh.steady()){
                  #include "pEqnSteady.H"
                }else{
                  #include "pEqn.H"
                }
              }
              turbulence.correct();
            }
          }
          fvMatrixAssemblyPtr->clear();
        }
      }
      runTime.write();
      runTime.printExecutionTime(Info);
    }
    runTime.writeNow();
    Info<<"End"<<endl;
    return 0;
}
