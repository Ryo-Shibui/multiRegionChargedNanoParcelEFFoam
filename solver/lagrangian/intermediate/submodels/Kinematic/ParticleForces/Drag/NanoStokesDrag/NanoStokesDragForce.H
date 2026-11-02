/* GPL-3.0-or-later. Stokes drag with Cunningham correction, SI units. */
#ifndef NanoStokesDragForce_H
#define NanoStokesDragForce_H
#include "ParticleForce.H"
namespace Foam
{
template<class CloudType>
class NanoStokesDragForce : public ParticleForce<CloudType>
{
    const scalar lambda_;
public:
    TypeName("nanoStokesDrag");
    NanoStokesDragForce(CloudType& owner, const fvMesh& mesh, const dictionary& dict)
    :
        ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
        lambda_
        (
            this->coeffs().found("meanFreePath")
          ? this->coeffs().template get<scalar>("meanFreePath")
          : this->coeffs().template get<scalar>("lambda")
        )
    {
        if (lambda_ < 0)
        {
            FatalErrorInFunction << "meanFreePath/lambda must be non-negative"
                << exit(FatalError);
        }
        if (this->coeffs().found("meanFreePath") && this->coeffs().found("lambda"))
        {
            FatalErrorInFunction << "Specify only meanFreePath or lambda"
                << exit(FatalError);
        }
    }
    NanoStokesDragForce(const NanoStokesDragForce& f)
    : ParticleForce<CloudType>(f), lambda_(f.lambda_) {}
    virtual autoPtr<ParticleForce<CloudType>> clone() const
    {
        return autoPtr<ParticleForce<CloudType>>(new NanoStokesDragForce(*this));
    }
    virtual forceSuSp calcCoupled
    (
        const typename CloudType::parcelType& p,
        const typename CloudType::parcelType::trackingData&,
        const scalar, const scalar mass, const scalar, const scalar muc
    ) const
    {
        const scalar Kn = 2*lambda_/p.d();
        const scalar Cc = Kn > VSMALL
          ? 1 + Kn*(1.257 + 0.4*exp(-1.1/Kn)) : 1;
        return forceSuSp(Zero, mass*18*muc/(p.rho()*sqr(p.d())*Cc));
    }
};
}
#endif
