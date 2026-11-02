/*---------------------------------------------------------------------------*\
    Dilute non-background species transport. GPL-3.0-or-later.
    Each executable owns this header; no standard OpenFOAM source is modified.
\*---------------------------------------------------------------------------*/
#ifndef diluteSpeciesDiffusion_H
#define diluteSpeciesDiffusion_H
#include "fvCFD.H"
#include "basicSpecieMixture.H"
#include <cmath>

namespace Foam
{
class diluteSpeciesDiffusion
{
    const fvMesh& mesh_;
    IOdictionary dict_;
    word model_;
    scalar turbulentSchmidt_;
    scalarField D_, Tref_, pref_, n_, m_;

    scalar value(const label i, const scalar T, const scalar p) const
    {
        scalar result=D_[i];
        if (n_[i] != 0)
        {
            if (!(T > 0)) FatalErrorInFunction
                << "Positive absolute temperature required by speciesDiffusionProperties"
                << exit(FatalError);
            result *= pow(T/Tref_[i], n_[i]);
        }
        if (m_[i] != 0)
        {
            if (!(p > 0)) FatalErrorInFunction
                << "Positive absolute pressure required by speciesDiffusionProperties"
                << exit(FatalError);
            result *= pow(pref_[i]/p, m_[i]);
        }
        if (!(result > 0) || !std::isfinite(result)) FatalErrorInFunction
            << "Non-positive/non-finite species diffusivity" << exit(FatalError);
        return result;
    }

public:
    diluteSpeciesDiffusion
    (
        const fvMesh& mesh, const basicSpecieMixture& composition,
        const word& background
    )
    : mesh_(mesh),
      dict_(IOobject("speciesDiffusionProperties",mesh.time().constant(),mesh,
          IOobject::READ_IF_PRESENT,IOobject::NO_WRITE,false)),
      model_(dict_.getOrDefault<word>("model",dict_.empty()?"legacy":"dilute")),
      turbulentSchmidt_(dict_.getOrDefault<scalar>("turbulentSchmidt",0.9)),
      D_(composition.species().size(),0), Tref_(D_.size(),298.15),
      pref_(D_.size(),101325), n_(D_.size(),0), m_(D_.size(),0)
    {
        if (model_=="legacy")
        {
            Info<< "Species diffusion in "<<mesh.name()
                <<": legacy muEff (Sc=ScT=1); configure speciesDiffusionProperties"
                <<" for independent dilute-species diffusivities."<<endl;
            return;
        }
        if (model_!="dilute") FatalIOErrorInFunction(dict_)
            <<"model must be dilute or legacy"<<exit(FatalIOError);
        if (!(turbulentSchmidt_>0) || !std::isfinite(turbulentSchmidt_))
            FatalIOErrorInFunction(dict_)<<"turbulentSchmidt must be positive and finite"
                <<exit(FatalIOError);
        if (!composition.species().found(background)) FatalIOErrorInFunction(dict_)
            <<"Background inertSpecie "<<background<<" is not a gas species"
                <<exit(FatalIOError);
        const dictionary& species=dict_.subDict("species");
        Info<<"Species diffusion in "<<mesh.name()<<": dilute, background "
            <<background<<", turbulentSchmidt "<<turbulentSchmidt_<<endl;
        forAll(composition.species(),i)
        {
            const word& name=composition.species()[i];
            if (name==background || !composition.active(i)) continue;
            const dictionary& coefficients=species.subDict(name);
            D_[i]=dimensionedScalar("D",dimArea/dimTime,coefficients).value();
            Tref_[i]=dimensionedScalar("Tref",dimTemperature,
                coefficients.getOrDefault<scalar>("Tref",298.15)).value();
            pref_[i]=coefficients.getOrDefault<scalar>("pref",101325);
            n_[i]=coefficients.getOrDefault<scalar>("temperatureExponent",0);
            m_[i]=coefficients.getOrDefault<scalar>("pressureExponent",0);
            if (!(D_[i]>0) || !(Tref_[i]>0) || !(pref_[i]>0)
                || !std::isfinite(D_[i]) || !std::isfinite(Tref_[i])
                || !std::isfinite(pref_[i]) || !std::isfinite(n_[i])
                || !std::isfinite(m_[i]))
                FatalIOErrorInFunction(coefficients)
                    <<"D, Tref and pref must be positive and finite; exponents finite"
                    <<exit(FatalIOError);
            Info<<"  "<<name<<": Dref="<<D_[i]<<" m2/s, Tref="<<Tref_[i]
                <<" K, pref="<<pref_[i]<<" Pa, temperatureExponent="<<n_[i]
                <<", pressureExponent="<<m_[i]<<endl;
        }
    }

    tmp<volScalarField> coefficient
    (
        const label i, const volScalarField& rho, const volScalarField& T,
        const volScalarField& p, const volScalarField& muEff,
        const volScalarField& mu
    ) const
    {
        if (model_=="legacy") return 1.0*muEff;
        auto tD=volScalarField::New
        (
            "diluteD_"+Foam::name(i), IOobject::NO_REGISTER, mesh_,
            dimensionedScalar("D",dimArea/dimTime,D_[i])
        );
        volScalarField& D=tD.ref();
        forAll(D,celli) D[celli]=value(i,T[celli],p[celli]);
        forAll(D.boundaryField(),patchi)
            forAll(D.boundaryField()[patchi],facei)
                D.boundaryFieldRef()[patchi][facei]=value
                (i,T.boundaryField()[patchi][facei],p.boundaryField()[patchi][facei]);
        // muEff-mu is the dynamic turbulent viscosity returned by the selected
        // turbulence model. D is an independent molecular property, not mu/rho.
        return rho*D+(muEff-mu)/turbulentSchmidt_;
    }
};
}
#endif
