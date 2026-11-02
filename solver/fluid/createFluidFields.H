// Initialise fluid field pointer lists
PtrList<rhoReactionThermo> thermoFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<compressible::turbulenceModel> turbulenceFluid(fluidRegions.size());
PtrList<CombustionModel<rhoReactionThermo>> reactionFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<radiation::radiationModel> radiation(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<volScalarField> dpdtFluid(fluidRegions.size());
PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>
    fieldsFluid(fluidRegions.size());
PtrList<volScalarField> QdotFluid(fluidRegions.size());

PtrList<fvVectorMatrix> UEqFluid(fluidRegions.size());

List<scalar> targetMassFluid(fluidRegions.size());
List<bool> frozenFlowFluid(fluidRegions.size(), false);
List<bool> solveSpeciesFluid(fluidRegions.size(), true);
List<bool> solveEnergyFluid(fluidRegions.size(), true);
List<bool> solveParcelsFluid(fluidRegions.size(), true);

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());

PtrList<dimensionedScalar> rhoMinFluid(fluidRegions.size());
PtrList<dimensionedScalar> rhoMaxFluid(fluidRegions.size());

PtrList<pressureControl> pressureControls(fluidRegions.size());

// Declare dependent objects after their carrier fields so that C++ destroys
// them first at solver shutdown (reverse declaration order).
PtrList<SLGThermo> slgThermoFluid(fluidRegions.size());
PtrList<basicReactingTypeCloud> parcelsFluid(fluidRegions.size());
PtrList<regionModels::surfaceFilmModel> surfaceFilmFluid(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoFluid\n" << endl;
    thermoFluid.set(i, rhoReactionThermo::New(fluidRegions[i]).ptr());
    thermoFluid[i].validate(args.executable(), "h", "e");

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].rho()
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(word::null, dimLength, Zero)
        )
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        compressible::turbulenceModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to reactionFluid\n" << endl;
    reactionFluid.set
    (
        i,
        CombustionModel<rhoReactionThermo>::New
        (
            thermoFluid[i],
            turbulenceFluid[i]
        )
    );

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // Force p_rgh to be consistent with p
    p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    Info<< "    Adding to radiationFluid\n" << endl;
    radiation.set
    (
        i,
        radiation::radiationModel::New(thermoFluid[i].T())
    );

    targetMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    Info<< "    Adding to dpdtFluid\n" << endl;
    dpdtFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                fluidRegions[i]
            ),
            fluidRegions[i],
            dimensionedScalar(thermoFluid[i].p().dimensions()/dimTime, Zero)
        )
    );

    Info<< "    Adding to fieldsFluid\n" << endl;
    fieldsFluid.set
    (
        i,
        new multivariateSurfaceInterpolationScheme<scalar>::fieldTable
    );
    forAll(thermoFluid[i].composition().Y(), j)
    {
        fieldsFluid[i].add(thermoFluid[i].composition().Y()[j]);
    }
    fieldsFluid[i].add(thermoFluid[i].he());

    Info<< "    Adding to QdotFluid\n" << endl;
    QdotFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
        )
    );

    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");

    // Preserve the useful split-solve controls from reactingNanoParcelFoam.
    // solveFlow takes precedence when present; otherwise retain the native
    // chtMultiRegionFoam frozenFlow control.
    bool solveFlow = true;
    if (pimpleDict.readIfPresent("solveFlow", solveFlow))
    {
        frozenFlowFluid[i] = !solveFlow;
    }
    else
    {
        pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);
    }

    pimpleDict.readIfPresent("solveSpecies", solveSpeciesFluid[i]);
    pimpleDict.readIfPresent("solveEnergy", solveEnergyFluid[i]);
    pimpleDict.readIfPresent("solveParcels", solveParcelsFluid[i]);

    if (fluidRegions[i].steady() && solveParcelsFluid[i])
    {
        FatalErrorInFunction
            << "solveParcels=true is not supported with steadyState ddt in "
            << "fluid region " << fluidRegions[i].name() << nl
            << "Use solveParcels false for the steady CHT/carrier solution. "
            << "Then restart from the converged fields with a transient ddt "
            << "scheme and solveParcels true for physical-time parcel tracking."
            << exit(FatalError);
    }

    rhoMaxFluid.set
    (
        i,
        new dimensionedScalar("rhoMax", dimDensity, GREAT, pimpleDict)
    );

    rhoMinFluid.set
    (
        i,
        new dimensionedScalar("rhoMin", dimDensity, Zero, pimpleDict)
    );

    pressureControls.set
    (
        i,
        new pressureControl(thermoFluid[i].p(), rhoFluid[i], pimpleDict, false)
    );

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );

    UEqFluid.set
    (
        i,
        new fvVectorMatrix(UFluid[i], dimForce)
    );

    turbulenceFluid[i].validate();

    // Construct the Lagrangian thermo/cloud only when explicitly enabled.
    // This is important for steady carrier/CHT calculations: when
    // solveParcels=false, no cloud dictionary, SLGThermo liquid/solid data,
    // injector, cloud function or parcel submodel is constructed.
    if (solveParcelsFluid[i])
    {
        Info<< "    Adding to SLGThermo\n" << endl;
        slgThermoFluid.set
        (
            i,
            new SLGThermo(fluidRegions[i], thermoFluid[i])
        );

        Info<< "    Constructing reacting nano cloud\n" << endl;
        parcelsFluid.set
        (
            i,
            new basicReactingTypeCloud
            (
                "reactingCloud1",
                rhoFluid[i],
                UFluid[i],
                g,
                slgThermoFluid[i]
            )
        );
    }
    else
    {
        Info<< "    Reacting nano cloud disabled (solveParcels false)\n"
            << endl;
    }

    // Keep the Eulerian surface-film model independent from solveParcels.
    // With the usual 'none' film this is a no-op; retaining it preserves the
    // original chtMultiRegionFoam/reactingParcelFoam film functionality.
    Info<< "    Constructing surface film model\n" << endl;
    surfaceFilmFluid.set
    (
        i,
        regionModels::surfaceFilmModel::New(fluidRegions[i], g).ptr()
    );

    pRefCellFluid[i] = -1;
    pRefValueFluid[i] = 0.0;

    if (p_rghFluid[i].needReference())
    {
        setRefCell
        (
            thermoFluid[i].p(),
            p_rghFluid[i],
            pimpleDict,
            pRefCellFluid[i],
            pRefValueFluid[i]
        );
    }

}
