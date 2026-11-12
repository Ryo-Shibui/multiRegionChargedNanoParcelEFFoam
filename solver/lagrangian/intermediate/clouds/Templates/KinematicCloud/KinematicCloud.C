/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "KinematicCloud.H"
#include "integrationScheme.H"
#include "interpolation.H"
#include "subCycleTime.H"

#include "InjectionModelList.H"
#include "DispersionModel.H"
#include "PatchInteractionModel.H"
#include "StochasticCollisionModel.H"
#include "SurfaceFilmModel.H"
#include "profiling.H"

#include "PackingModel.H"
#include "ParticleStressModel.H"
#include "DampingModel.H"
#include "IsotropyModel.H"
#include "TimeScaleModel.H"
#include "depositionPair.H"
#include "Pstream.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::KinematicCloud<CloudType>::setModels()
{
    dispersionModel_.reset
    (
        DispersionModel<KinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    patchInteractionModel_.reset
    (
        PatchInteractionModel<KinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    stochasticCollisionModel_.reset
    (
        StochasticCollisionModel<KinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    surfaceFilmModel_.reset
    (
        SurfaceFilmModel<KinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    packingModel_.reset
    (
        PackingModel<KinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    dampingModel_.reset
    (
        DampingModel<KinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    isotropyModel_.reset
    (
        IsotropyModel<KinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    UIntegrator_.reset
    (
        integrationScheme::New
        (
            "U",
            solution_.integrationSchemes()
        ).ptr()
    );
}


template<class CloudType>
template<class TrackCloudType>
void Foam::KinematicCloud<CloudType>::solve
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
    addProfiling(prof, "cloud::solve");

    log = solution_.log();

    if (solution_.steadyState())
    {
        cloud.storeState();

        cloud.preEvolve(td);

        evolveCloud(cloud, td);

        if (solution_.coupled())
        {
            cloud.relaxSources(cloud.cloudCopy());
        }
    }
    else
    {
        cloud.preEvolve(td);

        evolveCloud(cloud, td);

        if (solution_.coupled())
        {
            cloud.scaleSources();
        }
    }

    cloud.info();

    cloud.postEvolve(td);

    if (solution_.steadyState())
    {
        cloud.restoreState();
    }
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::buildCellOccupancy()
{
    if (!cellOccupancyPtr_)
    {
        cellOccupancyPtr_.reset
        (
            new List<DynamicList<parcelType*>>(mesh_.nCells())
        );
    }
    else if (cellOccupancyPtr_().size() != mesh_.nCells())
    {
        // If the size of the mesh has changed, reset the
        // cellOccupancy size

        cellOccupancyPtr_().setSize(mesh_.nCells());
    }

    List<DynamicList<parcelType*>>& cellOccupancy = cellOccupancyPtr_();

    for (auto& list : cellOccupancy)
    {
        list.clear();
    }

    for (parcelType& p : *this)
    {
        cellOccupancy[p.cell()].append(&p);
    }
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::updateCellOccupancy()
{
    // Only build the cellOccupancy if the pointer is set, i.e. it has
    // been requested before.

    if (cellOccupancyPtr_)
    {
        buildCellOccupancy();
    }
}


template<class CloudType>
template<class TrackCloudType>
void Foam::KinematicCloud<CloudType>::evolveCloud
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
    if (solution_.coupled())
    {
        cloud.resetSourceTerms();
    }

    if (solution_.transient())
    {
        label preInjectionSize = this->size();

        this->surfaceFilm().inject(cloud);

        // Update the cellOccupancy if the size of the cloud has changed
        // during the injection.
        if (preInjectionSize != this->size())
        {
            updateCellOccupancy();
            preInjectionSize = this->size();
        }

        injectors_.inject(cloud, td);

        // Assume that motion will update the cellOccupancy as necessary
        // before it is required.
        cloud.motion(cloud, td);

        stochasticCollision().update(td, solution_.trackTime());
    }
    else
    {
//        this->surfaceFilm().injectSteadyState(cloud);

        injectors_.injectSteadyState(cloud, td, solution_.trackTime());

        td.part() = parcelType::trackingData::tpLinearTrack;
        CloudType::move(cloud, td, solution_.trackTime());
    }
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    Log_<< endl;

    if (debug)
    {
        this->writePositions();
    }

    this->dispersion().cacheFields(false);

    this->patchInteraction().postEvolve();

    forces_.cacheFields(false);

    functions_.postEvolve(td);

    solution_.nextIter();

    if (this->db().time().writeTime())
    {
        outputProperties_.writeObject
        (
            IOstreamOption
            (
                IOstreamOption::ASCII,
                this->db().time().writeCompression()
            ),
            true
        );
    }

    if (this->dampingModel().active())
    {
        this->dampingModel().cacheFields(false);
    }
    if (this->packingModel().active())
    {
        this->packingModel().cacheFields(false);
    }
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::cloudReset(KinematicCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    rndGen_ = c.rndGen_;

    forces_.transfer(c.forces_);

    functions_.transfer(c.functions_);

    injectors_.transfer(c.injectors_);

    dispersionModel_.reset(c.dispersionModel_.ptr());
    patchInteractionModel_.reset(c.patchInteractionModel_.ptr());
    stochasticCollisionModel_.reset(c.stochasticCollisionModel_.ptr());
    surfaceFilmModel_.reset(c.surfaceFilmModel_.ptr());

    packingModel_.reset(c.packingModel_.ptr());
    dampingModel_.reset(c.dampingModel_.ptr());
    isotropyModel_.reset(c.isotropyModel_.ptr());

    UIntegrator_.reset(c.UIntegrator_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KinematicCloud<CloudType>::KinematicCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    bool readFields
)
:
    CloudType(rho.mesh(), cloudName, false),
    kinematicCloud(),
    cloudCopyPtr_(nullptr),
    mesh_(rho.mesh()),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_MODIFIED,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),
    outputProperties_
    (
        IOobject
        (
            cloudName + "OutputProperties",
            mesh_.time().timeName(),
            "uniform"/cloud::prefix/cloudName,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),
    solution_(mesh_, particleProperties_.subDict("solution")),
    constProps_(particleProperties_),
    subModelProperties_
    (
        particleProperties_.subOrEmptyDict
        (
            "subModels",
            keyType::REGEX,
            solution_.active()
        )
    ),
    rndGen_(Pstream::myProcNo()),
    cellOccupancyPtr_(),
    cellLengthScale_(mag(cbrt(mesh_.V()))),
    rho_(rho),
    U_(U),
    mu_(mu),
    g_(g),
    pAmbient_(0.0),
    forces_
    (
        *this,
        mesh_,
        subModelProperties_.subOrEmptyDict
        (
            "particleForces",
            keyType::REGEX,
            solution_.active()
        ),
        solution_.active()
    ),
    functions_
    (
        *this,
        particleProperties_.subOrEmptyDict("cloudFunctions"),
        solution_.active()
    ),
    injectors_
    (
        subModelProperties_.subOrEmptyDict("injectionModels"),
        *this
    ),
    dispersionModel_(nullptr),
    patchInteractionModel_(nullptr),
    stochasticCollisionModel_(nullptr),
    surfaceFilmModel_(nullptr),

    packingModel_(nullptr),
    dampingModel_(nullptr),
    isotropyModel_(nullptr),

    UIntegrator_(nullptr),
    rhokTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::scopedName(this->name(), "rhokTrans"),
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedScalar(dimMass, Zero)
        )
    ),
    UTrans_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                IOobject::scopedName(this->name(), "UTrans"),
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedVector(dimMass*dimVelocity, Zero)
        )
    ),
    UCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::scopedName(this->name(), "UCoeff"),
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedScalar(dimMass, Zero)
        )
    ),
    log(true)
{
    const dictionary deposition
    (
        particleProperties_.subOrEmptyDict("depositedParticleInteractions")
    );
    depositedParticleInteractions_ = deposition.getOrDefault("enabled", false);
    depositionMaxRelativeTravel_ =
        deposition.getOrDefault<scalar>("maxRelativeTravel", 0.1);
    if
    (
        depositedParticleInteractions_
     && (solution_.steadyState() || !(depositionMaxRelativeTravel_ > 0
         && depositionMaxRelativeTravel_ <= 1))
    )
    {
        FatalErrorInFunction
            << "depositedParticleInteractions requires transient tracking "
            << "and 0 < maxRelativeTravel <= 1"
            << exit(FatalError);
    }

    if (solution_.active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
            this->deleteLostParticles();
        }
    }

    if (solution_.resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::KinematicCloud<CloudType>::KinematicCloud
(
    KinematicCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c.mesh_, name, c),
    kinematicCloud(),
    cloudCopyPtr_(nullptr),
    mesh_(c.mesh_),
    particleProperties_(c.particleProperties_),
    outputProperties_(c.outputProperties_),
    solution_(c.solution_),
    constProps_(c.constProps_),
    subModelProperties_(c.subModelProperties_),
    rndGen_(c.rndGen_, true),
    cellOccupancyPtr_(nullptr),
    cellLengthScale_(c.cellLengthScale_),
    depositionSources_(c.depositionSources_),
    depositedParticleInteractions_(c.depositedParticleInteractions_),
    depositionMaxRelativeTravel_(c.depositionMaxRelativeTravel_),
    rho_(c.rho_),
    U_(c.U_),
    mu_(c.mu_),
    g_(c.g_),
    pAmbient_(c.pAmbient_),
    forces_(c.forces_),
    functions_(c.functions_),
    injectors_(c.injectors_),
    dispersionModel_(c.dispersionModel_->clone()),
    patchInteractionModel_(c.patchInteractionModel_->clone()),
    stochasticCollisionModel_(c.stochasticCollisionModel_->clone()),
    surfaceFilmModel_(c.surfaceFilmModel_->clone()),

    packingModel_(c.packingModel_->clone()),
    dampingModel_(c.dampingModel_->clone()),
    isotropyModel_(c.isotropyModel_->clone()),

    UIntegrator_(c.UIntegrator_->clone()),
    rhokTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::scopedName(this->name(), "rhokTrans"),
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            c.rhokTrans_()
        )
    ),
    UTrans_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                IOobject::scopedName(this->name(), "UTrans"),
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            c.UTrans_()
        )
    ),
    UCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::scopedName(this->name(), "UCoeff"),
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            c.UCoeff_()
        )
    ),
    log(c.log)
{}


template<class CloudType>
Foam::KinematicCloud<CloudType>::KinematicCloud
(
    const fvMesh& mesh,
    const word& name,
    const KinematicCloud<CloudType>& c
)
:
    CloudType(mesh, name, IDLList<parcelType>()),
    kinematicCloud(),
    cloudCopyPtr_(nullptr),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            name + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ),
    outputProperties_
    (
        IOobject
        (
            name + "OutputProperties",
            mesh_.time().timeName(),
            "uniform"/cloud::prefix/name,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ),
    solution_(mesh),
    constProps_(),
    subModelProperties_(),
    rndGen_(),
    cellOccupancyPtr_(nullptr),
    cellLengthScale_(c.cellLengthScale_),
    rho_(c.rho_),
    U_(c.U_),
    mu_(c.mu_),
    g_(c.g_),
    pAmbient_(c.pAmbient_),
    forces_(*this, mesh),
    functions_(*this),
    injectors_(*this),
    dispersionModel_(nullptr),
    patchInteractionModel_(nullptr),
    stochasticCollisionModel_(nullptr),
    surfaceFilmModel_(nullptr),

    packingModel_(nullptr),
    dampingModel_(nullptr),
    isotropyModel_(nullptr),

    UIntegrator_(nullptr),
    rhokTrans_(nullptr),
    UTrans_(nullptr),
    UCoeff_(nullptr),
    log(c.log)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::KinematicCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    // If rho0 is given in the const properties
    if (constProps_.rho0() != -1)
    {
        parcel.rho() = constProps_.rho0();
    }
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    const scalar carrierDt = mesh_.time().deltaTValue();
    parcel.stepFraction() = (carrierDt - lagrangianDt)/carrierDt;
    if (solution_.transient())
        parcel.setInitialStepFraction(parcel.stepFraction());

    if (parcel.typeId() == -1)
    {
        parcel.typeId() = constProps_.parcelTypeId();
    }

    if (parcel.rho() == -1)
    {
        FatalErrorInFunction
            << "The kinematic cloud needs rho0 in the constantProperties "
            << " dictionary. " << nl
            << abort(FatalError);
    }
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<KinematicCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::resetSourceTerms()
{
    rhokTrans().field() = Zero;
    UTrans().field() = Zero;
    UCoeff().field() = Zero;
}


template<class CloudType>
template<class Type>
void Foam::KinematicCloud<CloudType>::relax
(
    DimensionedField<Type, volMesh>& field,
    const DimensionedField<Type, volMesh>& field0,
    const word& name
) const
{
    const scalar coeff = solution_.relaxCoeff(name);
    field = field0 + coeff*(field - field0);
}


template<class CloudType>
template<class Type>
void Foam::KinematicCloud<CloudType>::scale
(
    DimensionedField<Type, volMesh>& field,
    const word& name
) const
{
    const scalar coeff = solution_.relaxCoeff(name);
    field *= coeff;
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::relaxSources
(
    const KinematicCloud<CloudType>& cloudOldTime
)
{
    this->relax(rhokTrans_(), cloudOldTime.rhokTrans(), "rhok");
    this->relax(UTrans_(), cloudOldTime.UTrans(), "U");
    this->relax(UCoeff_(), cloudOldTime.UCoeff(), "U");
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::scaleSources()
{
    this->scale(rhokTrans_(), "rhok");
    this->scale(UTrans_(), "U");
    this->scale(UCoeff_(), "U");
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::checkDepositionParcel
(
    const parcelType& p
) const
{
    if (p.nParticle() != 1 || !(p.d() > 0) || !std::isfinite(p.d()))
    {
        FatalErrorInFunction
            << "depositedParticleInteractions requires nParticle=1 and a "
            << "finite positive diameter; parcel "
            << p.origProc() << ':' << p.origId()
            << ", nParticle=" << p.nParticle() << ", d=" << p.d()
            << exit(FatalError);
    }
}


template<class CloudType>
bool Foam::KinematicCloud<CloudType>::isDepositedParticle(const parcelType& p) const
{
    if (p.active())
    {
        return false;
    }
    if (p.depositionSite() == parcelType::dpLegacy)
    {
        // Older field-based checkpoints have no deposition marker. The
        // boundary-face tetrahedron and coordinates survive positions I/O.
        const label facei = p.currentTetIndices().face();
        return facei >= mesh_.nInternalFaces()
            && mag(p.coordinates()[0]) <= SMALL
            && isA<wallPolyPatch>
               (mesh_.boundaryMesh()[mesh_.boundaryMesh().whichPatch(facei)]);
    }
    return p.depositionSite() == parcelType::dpWall
        || p.depositionSite() == parcelType::dpParticle;
}


template<class CloudType>
Foam::vector Foam::KinematicCloud<CloudType>::depositedParticleForce
(
    const parcelType& p
) const
{
    vector force(Zero);
    const scalar q = p.chargeNumber()*depositionPair::elementaryCharge;
    if (q == 0)
    {
        return force;
    }

    const point position = p.position();
    const scalar radius = 0.5*p.d();
    const scalar epsilonR = mesh_.lookupObject<volScalarField>("epsilonR")[p.cell()];
    for (const auto& rankSources : depositionSources_)
    {
        for (const auto& source : rankSources)
        {
            const scalar sourceCharge = source.second().second();
            if (sourceCharge != 0)
            {
                const vector separation = position - source.first();
                force += depositionPair::forceCoefficient
                (
                    q, sourceCharge, epsilonR, magSqr(separation),
                    sqr(radius + source.second().first())
                )*separation;
            }
        }
    }
    return force;
}


template<class CloudType>
Foam::scalar Foam::KinematicCloud<CloudType>::limitDepositionStep
(
    const parcelType& p,
    const point& start,
    const vector& displacement,
    bool& contact
) const
{
    const scalar radius = 0.5*p.d();
    const scalar displacementSquared = magSqr(displacement);
    scalar nearestSquared = GREAT;
    scalar firstHit = 2;
    for (const auto& rankSources : depositionSources_)
    {
        for (const auto& source : rankSources)
        {
            const vector separation = start - source.first();
            const scalar distanceSquared = magSqr(separation);
            nearestSquared = min(nearestSquared, distanceSquared);
            firstHit = min
            (
                firstHit,
                depositionPair::firstContact
                (
                    distanceSquared, separation & displacement,
                    displacementSquared, sqr(radius + source.second().first())
                )
            );
        }
    }

    scalar step = 1;
    if (displacementSquared > 0 && nearestSquared < GREAT)
    {
        step = min(step, depositionMaxRelativeTravel_
            *sqrt(nearestSquared)/sqrt(displacementSquared));
    }
    contact = firstHit <= step;
    return min(firstHit, step);
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    // force calculation of mesh dimensions - needed for parallel runs
    // with topology change due to lazy evaluation of valid mesh dimensions
    label nGeometricD = mesh_.nGeometricD();

    Log_<< "\nSolving" << nGeometricD << "-D cloud " << this->name() << endl;

    if (depositedParticleInteractions_)
    {
        DynamicList<depositSource> local;
        for (const auto& p : *this)
        {
            if (isDepositedParticle(p))
            {
                checkDepositionParcel(p);
                const scalar radius = 0.5*p.d();
                point centre = p.position();
                if (p.depositionSite() != parcelType::dpParticle)
                {
                    // Native stick stores the tracking point on the wall.
                    const vector area =
                        mesh_.faceAreas()[p.currentTetIndices().face()];
                    centre -= radius*area/mag(area);
                }
                local.append
                (
                    depositSource
                    (
                        centre,
                        Pair<scalar>(radius, p.chargeNumber()*depositionPair::elementaryCharge)
                    )
                );
            }
        }
        depositionSources_.setSize(Pstream::nProcs());
        depositionSources_[Pstream::myProcNo()].transfer(local);
        Pstream::allGatherList(depositionSources_);
    }

    this->dispersion().cacheFields(true);
    forces_.cacheFields(true);

    pAmbient_ = constProps_.dict().template
        getOrDefault<scalar>("pAmbient", pAmbient_);

    if (this->dampingModel().active() || this->packingModel().active())
    {
         const_cast<typename parcelType::trackingData&>(td).updateAverages(*this);
    }

    if (this->dampingModel().active())
    {
        this->dampingModel().cacheFields(true);
    }
    if (this->packingModel().active())
    {
        this->packingModel().cacheFields(true);
    }

    updateCellOccupancy();

    functions_.preEvolve(td);
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::evolve()
{
    if (solution_.canEvolve())
    {
        typename parcelType::trackingData td(*this);
        solve(*this, td);
    }
}


template<class CloudType>
template<class TrackCloudType>
void Foam::KinematicCloud<CloudType>::motion
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
    td.part() = parcelType::trackingData::tpLinearTrack;
    CloudType::move(cloud, td, solution_.trackTime());

    if (isotropyModel_->active())
    {
        td.updateAverages(cloud);
        isotropyModel_->calculate();
    }

    updateCellOccupancy();
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::patchData
(
    const parcelType& p,
    const polyPatch& pp,
    vector& nw,
    vector& Up
) const
{
    p.patchData(nw, Up);

    // If this is a wall patch, then there may be a non-zero tangential velocity
    // component; the lid velocity in a lid-driven cavity case, for example. We
    // want the particle to interact with this velocity, so we look it up in the
    // velocity field and use it to set the wall-tangential component.
    if (isA<wallPolyPatch>(pp))
    {
        const label patchi = pp.index();
        const label patchFacei = pp.whichFace(p.face());

        // We only want to use the boundary condition value  only if it is set
        // by the boundary condition. If the boundary values are extrapolated
        // (e.g., slip conditions) then they represent the motion of the fluid
        // just inside the domain rather than that of the wall itself.
        if (U_.boundaryField()[patchi].fixesValue())
        {
            const vector Uw1(U_.boundaryField()[patchi][patchFacei]);
            const vector& Uw0 =
                U_.oldTime().boundaryField()[patchi][patchFacei];

            const scalar f = p.currentTimeFraction();

            const vector Uw(Uw0 + f*(Uw1 - Uw0));

            const tensor nnw(nw*nw);

            Up = (nnw & Up) + Uw - (nnw & Uw);
        }
    }
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::updateMesh()
{
    updateCellOccupancy();
    injectors_.updateMesh();
    cellLengthScale_ = mag(cbrt(mesh_.V()));
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    updateMesh();
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::info()
{
    const vector linearMomentum =
        returnReduce(linearMomentumOfSystem(), sumOp<vector>());

    const scalar linearKineticEnergy =
        returnReduce(linearKineticEnergyOfSystem(), sumOp<scalar>());

    const label nTotParcel = returnReduce(this->size(), sumOp<label>());

    const scalar particlePerParcel =
    (
        nTotParcel
      ? (returnReduce(totalParticlePerParcel(), sumOp<scalar>()) / nTotParcel)
      : 0
    );

    Log_<< "Cloud: " << this->name() << nl
        << "    Current number of parcels       = " << nTotParcel << nl
        << "    Current mass in system          = "
        << returnReduce(massInSystem(), sumOp<scalar>()) << nl
        << "    Linear momentum                 = " << linearMomentum << nl
        << "   |Linear momentum|                = " << mag(linearMomentum) << nl
        << "    Linear kinetic energy           = " << linearKineticEnergy << nl
        << "    Average particle per parcel     = " << particlePerParcel << nl;


    injectors_.info();
    this->surfaceFilm().info();
    this->patchInteraction().info();

    if (this->packingModel().active())
    {
        tmp<volScalarField> alpha = this->theta();

        if (this->db().time().writeTime())
        {
            alpha().write();
        }

        const scalar alphaMin = gMin(alpha().primitiveField());
        const scalar alphaMax = gMax(alpha().primitiveField());

        Log_<< "    Min cell volume fraction        = " << alphaMin << nl
            << "    Max cell volume fraction        = " << alphaMax << endl;

        if (alphaMax < SMALL)
        {
            return;
        }

        scalar nMin = GREAT;

        forAll(this->mesh().cells(), celli)
        {
            const label n = this->cellOccupancy()[celli].size();

            if (n > 0)
            {
                const scalar nPack = n*alphaMax/alpha()[celli];

                if (nPack < nMin)
                {
                    nMin = nPack;
                }
            }
        }

        reduce(nMin, minOp<scalar>());

        Log_<< "    Min dense number of parcels     = " << nMin << endl;
    }
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::readObjects(const objectRegistry& obr)
{
    parcelType::readObjects(*this, obr);
}


template<class CloudType>
void Foam::KinematicCloud<CloudType>::writeObjects(objectRegistry& obr) const
{
    parcelType::writeObjects(*this, obr);
}


// ************************************************************************* //
