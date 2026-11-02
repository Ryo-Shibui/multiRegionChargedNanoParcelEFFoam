
#include "electricPotentialCoupled.H"
#include "mappedPatchBase.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
namespace Foam {
makePatchTypeField(fvPatchScalarField,electricPotentialCoupledFvPatchScalarField);
void electricPotentialCoupledFvPatchScalarField::updateCoeffs()
{
 if(updated()) return;
 const int oldTag=UPstream::incrMsgType();
 const auto& mp=refCast<const mappedPatchBase>(patch().patch());
 const auto& nm=refCast<const fvMesh>(mp.sampleMesh());
 const label ni=mp.samplePolyPatch().index();
 const auto& np=nm.boundary()[ni];
 const auto& nf=np.lookupPatchField<volScalarField>(internalField().name());
 scalarField nbr(nf.patchInternalField());
 scalarField nk(np.lookupPatchField<volScalarField>("epsilonR")*np.deltaCoeffs());
 mp.distribute(nbr);mp.distribute(nk);
 scalarField mine(patch().lookupPatchField<volScalarField>("epsilonR")*patch().deltaCoeffs());
 refValue()=nbr;refGrad()=0;valueFraction()=nk/(nk+mine);
 mixedFvPatchScalarField::updateCoeffs();
 UPstream::msgType(oldTag);
}
}
