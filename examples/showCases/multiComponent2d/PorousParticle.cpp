/** \file
 * Simulation of a 2D Rayleigh-Taylor instability, which describes a symmetry breakdown,
 * as a heavy fluid, initially located on top of a light one, starts penetrating the
 * latter. This is an illutration of the Shan/Chen multi-component model. Note that the
 * single-component multi-phase model is also implemented in Palabos, and a corresponding
 * sample program is provided. Also, note that the multi-component model can be used to
 * couple more than two phases, using the same approach as the one shown here.
 **/

#include "palabos2D.h"
#include "palabos2D.hh"
#include <cstdlib>
#include <iostream>

using namespace plb;
using namespace std;

// Use double-precision arithmetics
typedef double T;
// Use a grid which additionally to the f's stores two variables for
//   the external force term.
#define DESCRIPTOR descriptors::ForcedShanChenD2Q9Descriptor

plint Nx = 504;
plint Ny = 504;
T rhoH = 1.0;
T rhoL = 0.8;
plint hini = 50;
Box2D inlet(0, Nx-1, 0, 0);
Box2D outlet(0,Nx-1, Ny-1, Ny-1);

void rayleighTaylorSetup(MultiBlockLattice2D<T, DESCRIPTOR> &heavyFluid,
                         MultiBlockLattice2D<T, DESCRIPTOR> &lightFluid,
                         MultiScalarField2D<bool> &boolMask,
                         OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition,
                         T rho0, T rho1,
                         T vIn,
                         plint iT)
{
    Array<T,2> zeroVelocity((T)0.,(T)0.);
    T zeroDensity = (T)1.e-6;
    if (iT==0) ///initial
    {
    // Bounce-back
    defineDynamics(heavyFluid, boolMask, new BounceBack<T, DESCRIPTOR>(rho0), true);
    defineDynamics(lightFluid, boolMask, new BounceBack<T, DESCRIPTOR>(rho1), true);
    ///heavy:
    initializeAtEquilibrium(heavyFluid, Box2D(0,Nx-1, hini, Ny-1), rhoH, zeroVelocity);
    initializeAtEquilibrium(heavyFluid, Box2D(0,Nx-1, 0, hini-1), zeroDensity, zeroVelocity);
    //initializeAtEquilibrium(heavyFluid, Box2D(0,Nx-1, 0, 0), zeroDensity, zeroVelocity);

    ///light:
    initializeAtEquilibrium(lightFluid, Box2D(0,Nx-1, 0, hini-1), rhoL, zeroVelocity);
    initializeAtEquilibrium(lightFluid, Box2D(0,Nx-1, hini, Ny-1), zeroDensity, zeroVelocity);
    //initializeAtEquilibrium(lightFluid, Box2D(0,Nx-1, 0, 0), rhoL, zeroVelocity);

    heavyFluid.initialize();
    lightFluid.initialize();

    }
    if(iT>=1000)
    {
    ///heavy:
    //boundaryCondition->setPressureConditionOnBlockBoundaries(heavyFluid, inlet);
    //setBoundaryDensity(heavyFluid, inlet, zeroDensity);
    //initializeAtEquilibrium(heavyFluid, inlet, zeroDensity, zeroVelocity);

    boundaryCondition->addVelocityBoundary1P(outlet, heavyFluid, boundary::outflow);
    //initializeAtEquilibrium(heavyFluid, outlet, rhoH, zeroVelocity);

    ///light:
    boundaryCondition->addVelocityBoundary1N(inlet, lightFluid);
    setBoundaryVelocity(lightFluid, inlet, Array<T, 2>((T)0., vIn));
    //initializeAtEquilibrium(lightFluid, inlet, rhoL, Array<T, 2>((T)0., vIn));

    boundaryCondition->addVelocityBoundary1P(outlet, lightFluid, boundary::outflow);
    //initializeAtEquilibrium(lightFluid, outlet, zeroDensity, zeroVelocity);

    heavyFluid.initialize();
    lightFluid.initialize();

    }
}

void writeGifs(MultiBlockLattice2D<T, DESCRIPTOR> &heavyFluid,
               MultiBlockLattice2D<T, DESCRIPTOR> &lightFluid, plint iT)
{
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("v_heavy_", iT, 6),
                               *computeVelocityNorm(heavyFluid));
    imageWriter.writeScaledGif(createFileName("rho_heavy_", iT, 6),
                               *computeDensity(heavyFluid));
    imageWriter.writeScaledGif(createFileName("v_light_", iT, 6),
                               *computeVelocityNorm(lightFluid));
    imageWriter.writeScaledGif(createFileName("rho_light_", iT, 6),
                               *computeDensity(lightFluid));
}
void writeVTKs(MultiBlockLattice2D<T, DESCRIPTOR> &heavyFluid,
               MultiBlockLattice2D<T, DESCRIPTOR> &lightFluid, plint iT)
{
    T dx = 1.0/50.0;
    T dt = dx/10.0;
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iT, 6), dx);
    vtkOut.writeData<2,float>(*computeVelocity(heavyFluid), "velocity_H", dx/dt);
    vtkOut.writeData<float>(*computeDensity(heavyFluid), "density_H", 1.);
    vtkOut.writeData<2,float>(*computeVelocity(lightFluid), "velocity_L", dx/dt);
    vtkOut.writeData<float>(*computeDensity(lightFluid), "density_L", 1.);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    srand(global::mpi().getRank());
    util::delete_path("./tmp/");
    const T vIn = 1.e-2;
    const T G = 2.2;
    T rho0 = 0.6;
    T rho1 = 1-rho0;
    const plint maxIter = 50000;
    const plint saveIter = 100;
    const plint statIter = 100;

    const T omega1 = 1.;
    const T omega2 = 1.;

    MultiScalarField2D<bool> boolMask(Nx, Ny);

    plb_ifstream ifile("PPGeo.dat");
    ifile >> boolMask;

    // Use regularized BGK dynamics to improve numerical stability (but note that
    //   BGK dynamics works well too).
    MultiBlockLattice2D<T, DESCRIPTOR> heavyFluid(
        Nx, Ny, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega1));
    MultiBlockLattice2D<T, DESCRIPTOR> lightFluid(
        Nx, Ny, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega2));

    // Boundary conditions.
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *
        boundaryCondition = createLocalBoundaryCondition2D<T, DESCRIPTOR>();

    heavyFluid.periodicity().toggle(0,true);
    lightFluid.periodicity().toggle(0,true);
    //heavyFluid.periodicity().toggle(1,true);
    //lightFluid.periodicity().toggle(1,true);

    // Store a pointer to all lattices (two in the present application) in a vector to
    //   create the Shan/Chen coupling therm. The heavy fluid being at the first place
    //   in the vector, the coupling term is going to be executed at the end of the call
    //   to collideAndStream() or stream() for the heavy fluid.*/
    vector<MultiBlockLattice2D<T, DESCRIPTOR> *> blockLattices;
    blockLattices.push_back(&heavyFluid);
    blockLattices.push_back(&lightFluid);
    std::vector<T> constOmegaValues;
    constOmegaValues.push_back(omega1);
    constOmegaValues.push_back(omega2);
    plint processorLevel = 1;
    integrateProcessingFunctional(
        new ShanChenMultiComponentProcessor2D<T,DESCRIPTOR> (G,constOmegaValues),
            Box2D(0, Nx-1, 0, Ny-1), blockLattices, processorLevel);

    pcout << "Starting simulation" << endl;
    // Main loop over time iterations.
    for (plint iT = 0; iT < maxIter; ++iT)
    {
        rayleighTaylorSetup(heavyFluid, lightFluid, boolMask, boundaryCondition, rho0, rho1, vIn, iT);
        heavyFluid.collideAndStream();
        lightFluid.collideAndStream();
        if (iT % saveIter == 0)
        {
            writeGifs(heavyFluid, lightFluid, iT);
            writeVTKs(heavyFluid, lightFluid, iT);
        }
        // Time iteration for the heavy fluid must come after the light fluid,
        //   because the coupling is executed here. You should understand this as follows.
        //   The effect of the coupling is to compute the interaction force between
        //   species, and to precompute density and momentum for each species. This must
        //   be executed *before* collide-and-streaming the fluids, because the collision
        //   step needs to access all these values. In the present case, it is done after
        //   both collide-and-stream step, which means, before the collide-and-stream of
        //   the next iteration (it's the same if you are before or after; the important
        //   point is not to be between the two collide-and-streams of the light and heavy
        //   fluid. As for the initial condition, the coupling is initially performed once
        //   during the function call to heavyFluid.initialize().

        if (iT % statIter == 0)
        {
            pcout << "Average density fluid one = "
                  << getStoredAverageDensity<T>(heavyFluid);
            pcout << ", average density fluid two = "
                  << getStoredAverageDensity<T>(lightFluid) << endl;
        }
    }
    delete boundaryCondition;
}
