// GPL-3.0-or-later. Standalone checks for the solver's pair kernels.
#include "../../solver/lagrangian/intermediate/clouds/Templates/KinematicCloud/depositionPair.H"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace pair = Foam::depositionPair;

void require(bool condition, const char* message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

struct Landing
{
    double x, y, time;
    bool particle;
};

// Same deterministic split as the native tracker, in a uniform field and
// stationary gas. The single deposit is fixed at (0, radius).
Landing trace
(
    double x,
    double seedCharge,
    bool enabled,
    double dt,
    std::ostream* trajectory = nullptr
)
{
    constexpr double radius = 24e-9;
    constexpr double charge = 300*pair::elementaryCharge;
    constexpr double field = -1.25e6;
    constexpr double viscosity = 1.81e-5;
    constexpr double diameter = 2*radius;
    constexpr double pi = 3.141592653589793;
    constexpr double lambda = 6.5e-8;
    const double cc = 1 + 2*lambda/diameter
        *(1.257 + 0.4*std::exp(-1.1*diameter/(2*lambda)));
    const double drag = 3*pi*viscosity*diameter/cc;
    const double tau = 1000*pi*diameter*diameter*diameter/(6*drag);
    double y = 1.2e-6;
    double ux = 0;
    double uy = charge*field/drag;
    double time = 0;
    const double offset = x;
    while (time < 4e-7)
    {
        if (trajectory)
        {
            *trajectory << offset*1e9 << ',' << enabled << ',' << time*1e9
                << ',' << x*1e9 << ',' << y*1e9 << '\n';
        }
        const double dx = dt*ux;
        const double dy = dt*uy;
        const double wall = dy < 0 && y+dy <= 0 ? -y/dy : 2.0;
        const double sepY = y-radius;
        const double distance2 = x*x + sepY*sepY;
        const double displacement2 = dx*dx + dy*dy;
        const double pairHit = enabled ? pair::firstContact
        (
            distance2, x*dx + sepY*dy, displacement2, diameter*diameter
        ) : 2.0;
        const double limit = enabled && displacement2 > 0
            ? 0.1*std::sqrt(distance2/displacement2) : 2.0;
        const double fraction = std::min({1.0, wall, pairHit, limit});
        x += fraction*dx;
        y += fraction*dy;
        time += fraction*dt;
        if (pairHit <= fraction || wall <= fraction)
        {
            if (trajectory)
            {
                *trajectory << offset*1e9 << ',' << enabled << ',' << time*1e9
                    << ',' << x*1e9 << ',' << y*1e9 << '\n';
            }
            return {x, y, time, pairHit <= fraction};
        }
        const double r2 = x*x + (y-radius)*(y-radius);
        const double k = enabled ? pair::forceCoefficient
        (
            charge, seedCharge, 1.0, r2, diameter*diameter
        ) : 0.0;
        const double relaxation = -std::expm1(-fraction*dt/tau);
        ux += relaxation*(k*x/drag - ux);
        uy += relaxation*((charge*field + k*(y-radius))/drag - uy);
    }
    throw std::runtime_error("trajectory did not land");
}

int main(int argc, char** argv)
{
    require(std::abs(pair::firstContact(25., -50., 100., 1.) - .4) < 1e-14,
        "a segment crossing a sphere must stop at its first surface");
    require(std::abs(pair::firstContact(26., -50., 100., 1.) - .5) < 1e-14,
        "a tangent must count as contact");
    require(pair::firstContact(29., -50., 100., 1.) > 1.,
        "a segment missing a sphere must continue");
    require(pair::firstContact(1., -1., 1., 4.) == 0.,
        "an initial overlap must be caught");
    require(pair::firstContact(25., -5., 1., 1.) > 1.,
        "an intersection beyond the step must be ignored");
    require(pair::forceCoefficient(2., 3., 1., 4., 1.) > 0.,
        "like charges must repel");
    const double force = pair::forceCoefficient(2., -3., 2., 4., 1.)*2;
    require(std::abs(force/(-6.740663844195878e9) - 1.) < 1e-12,
        "opposite charges and permittivity must match Coulomb's law");

    std::ofstream trajectory;
    if (argc > 1)
    {
        trajectory.open(argv[1]);
        require(bool(trajectory), "could not open trajectory output");
        trajectory << std::setprecision(12) << "offset_nm,enabled,t_ns,x_nm,y_nm\n";
    }
    auto* output = trajectory.is_open() ? &trajectory : nullptr;
    std::cout << std::setprecision(12)
        << "offset_nm,off_x_nm,on_x_nm,on_minus_off_nm,off_t_ns,on_t_ns,refined_shift_nm\n";
    for (double offset : {250e-9, 500e-9, 1000e-9})
    {
        const Landing off = trace(offset, 300*pair::elementaryCharge, false, 5e-11, output);
        const Landing on = trace(offset, 300*pair::elementaryCharge, true, 5e-11, output);
        const Landing refined = trace(offset, 300*pair::elementaryCharge, true, 2.5e-11);
        require(!off.particle && !on.particle, "repulsive cases must reach the wall");
        require(std::abs(off.x-offset) < 1e-15, "unforced lateral motion");
        require(on.x > off.x + 1e-9, "same-sign deposit must deflect the trajectory");
        require(std::abs(on.x-refined.x) < 2e-9, "landing must converge with smaller dt");
        std::cout << offset*1e9 << ',' << off.x*1e9 << ',' << on.x*1e9
            << ',' << (on.x-off.x)*1e9 << ',' << off.time*1e9 << ','
            << on.time*1e9 << ',' << (refined.x-on.x)*1e9 << '\n';
    }
    const Landing hit = trace(0, 0, true, 5e-11);
    require(hit.particle && std::abs(hit.y-72e-9) < 1e-15,
        "neutral seed must still be contacted at the sum of radii");
    std::cout << "neutral-seed contact y_nm=" << hit.y*1e9 << '\n';
}
