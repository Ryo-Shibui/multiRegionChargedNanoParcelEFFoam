// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef depositionPair_H
#define depositionPair_H

#include <algorithm>
#include <cmath>

namespace Foam
{
namespace depositionPair
{

constexpr double elementaryCharge = 1.602176634e-19;
constexpr double coulombConstant = 8.987551792261171e9;

// Multiply by the centre-to-centre displacement. The lower bound
// protects against roundoff at a contact endpoint.
template<class Scalar>
Scalar forceCoefficient
(
    Scalar q,
    Scalar sourceCharge,
    Scalar epsilonR,
    Scalar distanceSquared,
    Scalar contactRadiusSquared
)
{
    const Scalar r2 = std::max(distanceSquared, contactRadiusSquared);
    return Scalar(coulombConstant)*q*sourceCharge/(epsilonR*r2*std::sqrt(r2));
}


// First intersection along a straight segment, parameterized from 0 to 1.
// Return 2 for no contact within the segment.
template<class Scalar>
Scalar firstContact
(
    Scalar distanceSquared,
    Scalar displacementDotSeparation,
    Scalar displacementSquared,
    Scalar contactRadiusSquared
)
{
    const Scalar c = distanceSquared - contactRadiusSquared;
    if (c <= 0)
    {
        return 0;
    }
    const Scalar b = displacementDotSeparation;
    if (b >= 0 || displacementSquared <= 0)
    {
        return 2;
    }
    const Scalar discriminant = std::fma(-displacementSquared, c, b*b);
    if (discriminant < 0)
    {
        return 2;
    }
    const Scalar contact = c/(-b + std::sqrt(discriminant));
    return contact <= 1 ? contact : Scalar(2);
}

}
}
#endif
