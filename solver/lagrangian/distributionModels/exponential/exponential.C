/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "exponential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mrChargedEFDistributionModels
{
    defineTypeNameAndDebug(exponential, 0);
    addToRunTimeSelectionTable(mrChargedEFDistributionModel, exponential, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mrChargedEFDistributionModels::exponential::exponential
(
    const dictionary& dict,
    Random& rndGen
)
:
    mrChargedEFDistributionModel(typeName, dict, rndGen),
    lambda_(distributionModelDict_.get<scalar>("lambda"))
{
    if (lambda_ < VSMALL)
    {
        FatalErrorInFunction
            << "Rate parameter cannot be equal to or less than zero:" << nl
            << "    lambda = " << lambda_
            << exit(FatalError);
    }

    check();
}


Foam::mrChargedEFDistributionModels::exponential::exponential(const exponential& p)
:
    mrChargedEFDistributionModel(p),
    lambda_(p.lambda_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::mrChargedEFDistributionModels::exponential::sample() const
{
    const scalar u = rndGen_.sample01<scalar>();
    const scalar qMin = exp(-lambda_*minValue_);
    const scalar qMax = exp(-lambda_*maxValue_);
    return -(scalar(1)/lambda_)*log(qMin + u*(qMax - qMin));
}


Foam::scalar Foam::mrChargedEFDistributionModels::exponential::meanValue() const
{
    return scalar(1)/lambda_;
}


// ************************************************************************* //
