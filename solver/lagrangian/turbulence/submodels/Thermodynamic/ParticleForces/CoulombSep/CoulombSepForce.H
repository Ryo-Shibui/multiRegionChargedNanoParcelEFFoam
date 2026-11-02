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

#ifndef Foam_CoulombSepForce_H
#define Foam_CoulombSepForce_H

#include "ParticleForce.H"
#include "interpolation.H"

namespace Foam
{

template<class CloudType>
class CoulombSepForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Name of electric field
        const word EName_;

        //- Electric-field interpolator
        const interpolation<vector>* EInterpPtr_;


public:

    //- Runtime type information
    TypeName("CoulombSep");


    // Constructors

        //- Construct from mesh
        CoulombSepForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Copy construct
        CoulombSepForce(const CoulombSepForce& pf);

        //- No copy assignment
        void operator=(const CoulombSepForce<CloudType>&) = delete;

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new CoulombSepForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CoulombSepForce();


    // Member Functions

        //- Cache fields
        virtual void cacheFields(const bool store);

        //- Calculate the non-coupled force
        virtual forceSuSp calcNonCoupled
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td,
            const scalar dt,
            const scalar mass,
            const scalar Re,
            const scalar muc
        ) const;
};

} // End namespace Foam

#ifdef NoRepository
    #include "CoulombSepForce.C"
#endif

#endif

// ************************************************************************* //
