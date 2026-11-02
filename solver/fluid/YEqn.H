tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);

if (Y.size())
{
    mvConvection = tmp<fv::convectionScheme<scalar>>
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );
}

{
    reaction.correct();
    Qdot = reaction.Qdot();
    volScalarField Yt
    (
        IOobject("Yt", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("Yt", dimless, 0)
    );

    forAll(Y, specieI)
    {
        if (specieI != inertIndex && composition.active(specieI))
        {
            volScalarField& Yi = Y[specieI];

            tmp<volScalarField> gammaSpecies=speciesDiffusionFluid[i].coefficient
            (specieI,rho,thermo.T(),thermo.p(),turbulence.muEff(),thermo.mu());

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian
                (
                    gammaSpecies(), Yi,
                    "laplacian(muEff," + Yi.name() + ")"
                )
             ==
                surfaceFilm.Srho(specieI)
              + reaction.R(Yi)
              + fvOptions(rho, Yi)
            );

            if (solveParcels)
            {
                YiEqn -= parcelsPtr->SYi(specieI, Yi);
            }

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.clamp_min(0);
            Yt += Yi;
        }
    }

    if (Y.size())
    {
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].clamp_min(0);
    }
}
