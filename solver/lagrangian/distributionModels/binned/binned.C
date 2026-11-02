/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "binned.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mrChargedEFDistributionModels
{
    defineTypeNameAndDebug(binned, 0);
    addToRunTimeSelectionTable(mrChargedEFDistributionModel, binned, dictionary);
}
}


const char* Foam::mrChargedEFDistributionModels::binned::header =
    "(bin probability)";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mrChargedEFDistributionModels::binned::initialise()
{
    const label nSample(xy_.size());

    // Convert values to integral values
    for (label bini = 1; bini < nSample; ++bini)
    {
        xy_[bini][1] += xy_[bini - 1][1];
    }

    // Normalise
    scalar sumProb = xy_.last()[1];

    if (sumProb < VSMALL)
    {
        FatalErrorInFunction
            << type() << "distribution: "
            << "The sum of elements in the second column cannot be zero." << nl
            << "sum = " << sumProb
            << exit(FatalError);
    }

    forAll(xy_, bini)
    {
        xy_[bini][1] /= sumProb;
    }

    // Calculate the mean value
    label bini = 0;
    forAll(xy_, i)
    {
        if (xy_[i][1] > 0.5)
        {
            bini = i;
            break;
        }
    }

    meanValue_ = xy_[bini][1];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mrChargedEFDistributionModels::binned::binned
(
    const dictionary& dict,
    Random& rndGen
)
:
    mrChargedEFDistributionModel(typeName, dict, rndGen),
    xy_(distributionModelDict_.lookup("distribution")),
    meanValue_(0)
{
    minValue_ = xy_[0][0];
    maxValue_ = xy_[xy_.size()-1][0];

    check();

    initialise();
}


Foam::mrChargedEFDistributionModels::binned::binned
(
    const UList<scalar>& sampleData,
    const scalar binWidth,
    Random& rndGen
)
:
    mrChargedEFDistributionModel(typeName, dictionary::null, rndGen),
    xy_(),
    meanValue_(0)
{
    minValue_ = GREAT;
    maxValue_ = -GREAT;
    forAll(sampleData, i)
    {
        minValue_ = min(minValue_, sampleData[i]);
        maxValue_ = max(maxValue_, sampleData[i]);
    }

    const label bin0 = floor(minValue_/binWidth);
    const label bin1 = ceil(maxValue_/binWidth);
    const label nBin = bin1 - bin0;

    if (nBin == 0)
    {
        WarningInFunction
            << "Data cannot be binned - zero bins generated" << nl
            << "   Bin width   : " << binWidth << nl
            << "   Sample data : " << sampleData
            << endl;

        return;
    }

    // Populate bin boundaries and initialise occurrences
    xy_.setSize(nBin);
    forAll(xy_, bini)
    {
        xy_[bini][0] = (bin0 + bini)*binWidth;
        xy_[bini][1] = 0;
    }

    // Bin the data
    forAll(sampleData, i)
    {
        // Choose the nearest bin
        label bini = floor(sampleData[i]/binWidth) - bin0;
        label binii = min(bini + 1, nBin - 1);

        scalar d1 = mag(sampleData[i] - xy_[bini][0]);
        scalar d2 = mag(xy_[binii][0] - sampleData[i]);

        if (d1 < d2)
        {
            xy_[bini][1]++;
        }
        else
        {
            xy_[binii][1]++;
        }
    }

    initialise();
}


Foam::mrChargedEFDistributionModels::binned::binned(const binned& p)
:
    mrChargedEFDistributionModel(p),
    xy_(p.xy_),
    meanValue_(p.meanValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::mrChargedEFDistributionModels::binned::sample() const
{
    const scalar u = rndGen_.sample01<scalar>();

    for (label i = 0; i < xy_.size() - 1; ++i)
    {
        if (xy_[i][1] > u)
        {
            return xy_[i][0];
        }
    }

    return maxValue_;
}


Foam::scalar Foam::mrChargedEFDistributionModels::binned::meanValue() const
{
    return meanValue_;
}


void Foam::mrChargedEFDistributionModels::binned::readData(Istream& is)
{
//    mrChargedEFDistributionModel::readData(is);
    is  >> xy_;
}


void Foam::mrChargedEFDistributionModels::binned::writeData(Ostream& os) const
{
//    mrChargedEFDistributionModel::writeData(os);
    os  << xy_ ;
}


Foam::dictionary Foam::mrChargedEFDistributionModels::binned::writeDict
(
    const word& dictName
) const
{
//    dictionary dict = mrChargedEFDistributionModel::writeDict(dictName);
    dictionary dict(dictName);
    dict.add("distribution", xy_);

    return dict;
}


void Foam::mrChargedEFDistributionModels::binned::readDict(const dictionary& dict)
{
//    mrChargedEFDistributionModel::readDict(dict);
    dict.readEntry("distribution", xy_);
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const mrChargedEFDistributionModels::binned& b
)
{
    os.check(FUNCTION_NAME);

    b.writeData(os);
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, mrChargedEFDistributionModels::binned& b)
{
    is.check(FUNCTION_NAME);

    b.readData(is);
    return is;
}


// ************************************************************************* //
