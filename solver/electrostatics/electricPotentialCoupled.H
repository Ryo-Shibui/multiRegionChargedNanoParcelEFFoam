
// GPL-3.0-or-later. Mixed interface follows OpenFOAM coupled baffle formulation.
#ifndef electricPotentialCoupled_H
#define electricPotentialCoupled_H
#include "mixedFvPatchFields.H"
namespace Foam {
class electricPotentialCoupledFvPatchScalarField : public mixedFvPatchScalarField
{
public:
 TypeName("electricPotentialCoupled");
 electricPotentialCoupledFvPatchScalarField(const fvPatch& p,const DimensionedField<scalar,volMesh>& f)
 :mixedFvPatchScalarField(p,f) {refValue()=0;refGrad()=0;valueFraction()=1;}
 electricPotentialCoupledFvPatchScalarField(const fvPatch& p,const DimensionedField<scalar,volMesh>& f,const dictionary& d)
 :mixedFvPatchScalarField(p,f) {readValueEntry(d,IOobjectOption::MUST_READ);refValue()=*this;refGrad()=0;valueFraction()=1;}
 electricPotentialCoupledFvPatchScalarField(const electricPotentialCoupledFvPatchScalarField& b,const fvPatch& p,const DimensionedField<scalar,volMesh>& f,const fvPatchFieldMapper& m)
 :mixedFvPatchScalarField(b,p,f,m){}
 electricPotentialCoupledFvPatchScalarField(const electricPotentialCoupledFvPatchScalarField& b)
 :mixedFvPatchScalarField(b){}
 electricPotentialCoupledFvPatchScalarField(const electricPotentialCoupledFvPatchScalarField& b,const DimensionedField<scalar,volMesh>& f)
 :mixedFvPatchScalarField(b,f){}
 virtual tmp<fvPatchScalarField> clone() const {return tmp<fvPatchScalarField>(new electricPotentialCoupledFvPatchScalarField(*this));}
 virtual tmp<fvPatchScalarField> clone(const DimensionedField<scalar,volMesh>& f) const {return tmp<fvPatchScalarField>(new electricPotentialCoupledFvPatchScalarField(*this,f));}
 virtual void updateCoeffs();
 virtual void write(Ostream& os) const {mixedFvPatchScalarField::write(os);}
};
}
#endif
