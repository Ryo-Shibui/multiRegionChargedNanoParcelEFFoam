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

#include "uniform.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mrChargedEFDistributionModels
{
    defineTypeNameAndDebug(uniform, 0);
    addToRunTimeSelectionTable(mrChargedEFDistributionModel, uniform, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mrChargedEFDistributionModels::uniform::uniform
(
    const dictionary& dict,
    Random& rndGen
)
:
    mrChargedEFDistributionModel(typeName, dict, rndGen)
{
    check();
}


Foam::mrChargedEFDistributionModels::uniform::uniform(const uniform& p)
:
    mrChargedEFDistributionModel(p)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::mrChargedEFDistributionModels::uniform::sample() const
{
    return rndGen_.position<scalar>(minValue_, maxValue_);
}


Foam::scalar Foam::mrChargedEFDistributionModels::uniform::meanValue() const
{
    return 0.5*(minValue_ + maxValue_);
}


// ************************************************************************* //
