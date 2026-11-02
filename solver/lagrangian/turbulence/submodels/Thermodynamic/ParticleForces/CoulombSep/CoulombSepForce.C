/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "CoulombSepForce.H"
#include "demandDrivenData.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CoulombSepForce<CloudType>::CoulombSepForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    EName_(this->coeffs().template getOrDefault<word>("E", "E")),
    EInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::CoulombSepForce<CloudType>::CoulombSepForce
(
    const CoulombSepForce& pf
)
:
    ParticleForce<CloudType>(pf),
    EName_(pf.EName_),
    EInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CoulombSepForce<CloudType>::~CoulombSepForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CoulombSepForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const volVectorField& E =
            this->mesh().template lookupObject<volVectorField>(EName_);

        EInterpPtr_ = interpolation<vector>::New
        (
            this->owner().solution().interpolationSchemes(),
            E
        ).ptr();
    }
    else
    {
        deleteDemandDrivenData(EInterpPtr_);
    }
}


template<class CloudType>
Foam::forceSuSp Foam::CoulombSepForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const interpolation<vector>& EInterp = *EInterpPtr_;

    const vector E =
        EInterp.interpolate(p.coordinates(), p.currentTetIndices());
    const scalar elementaryCharge = 1.602176634e-19;

    return forceSuSp
    (
        p.chargeNumber()*elementaryCharge*E,
        Zero
    );
}


// ************************************************************************* //
