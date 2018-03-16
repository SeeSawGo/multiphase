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
#include <dirent.h>

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
Box2D inlet(0, Nx-1, 0, 0);
Box2D outlet(0,Nx-1, Ny-1, Ny-1);

void rayleighTaylorSetup(MultiBlockLattice2D<T, DESCRIPTOR> &heavyFluid,
                         MultiScalarField2D<bool> &boolMask,
                         OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition,
                         T vIn,
                         plint iT)
{
    Array<T,2> zeroVelocity((T)0.,(T)0.);
    T zeroDensity = (T)1.e-6;
    if (iT==0) ///initial
    {
    // Bounce-back
    defineDynamics(heavyFluid, boolMask, new BounceBack<T, DESCRIPTOR>, true);
    ///heavy:
    initializeAtEquilibrium(heavyFluid, heavyFluid.getBoundingBox(), rhoH, zeroVelocity);
    //initializeAtEquilibrium(heavyFluid, Box2D(0,Nx-1, 0, 199), zeroDensity, zeroVelocity);
    //initializeAtEquilibrium(heavyFluid, Box2D(0,Nx-1, 0, 0), zeroDensity, zeroVelocity);

    heavyFluid.initialize();

    }
    if(iT==0)
    {
    ///heavy:
    boundaryCondition->addVelocityBoundary1N(inlet, heavyFluid);
    initializeAtEquilibrium(heavyFluid, inlet, rhoH, Array<T,2>((T)0.0, vIn));

    boundaryCondition->addVelocityBoundary1P(outlet, heavyFluid, boundary::outflow);
    initializeAtEquilibrium(heavyFluid, outlet, rhoH, zeroVelocity);

    //defineDynamics(heavyFluid, boolMask, new BounceBack<T, DESCRIPTOR>, true);
    heavyFluid.initialize();

    }
}

void writeGifs(MultiBlockLattice2D<T, DESCRIPTOR> &heavyFluid, plint iT)
{
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("v_heavy_", iT, 6),
                               *computeVelocityNorm(heavyFluid));
    imageWriter.writeScaledGif(createFileName("rho_heavy_", iT, 6),
                               *computeDensity(heavyFluid));
}
void writeVTKs(MultiBlockLattice2D<T, DESCRIPTOR> &heavyFluid, plint iT)
{
    T dx = 1.0/50.0;
    T dt = dx/10.0;
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iT, 6), dx);
    vtkOut.writeData<2,float>(*computeVelocity(heavyFluid), "velocity_H", dx/dt);
    vtkOut.writeData<float>(*computeDensity(heavyFluid), "density_H", 1.);
}

void delete_path(const char *path)
{

    DIR *pDir = NULL;
    struct dirent *dmsg;
    char szFileName[128];
    char szFolderName[128];

    strcpy(szFolderName, path);
    strcat(szFolderName, "/%s");
    if ((pDir = opendir(path)) != NULL)
    {
        // 遍历目录并删除文件
        while ((dmsg = readdir(pDir)) != NULL)
        {
            if (strcmp(dmsg->d_name, ".") != 0 && strcmp(dmsg->d_name, "..") != 0)
            {
                sprintf(szFileName, szFolderName, dmsg->d_name);
                string tmp = szFileName;
                //如果是文件夹，名称中不包含"."
                if (tmp.find(".") == -1)
                {
                    delete_path(szFileName);
                }
                remove(szFileName);
            }
        }
    }
    if (pDir != NULL)
    {
        closedir(pDir);
    }
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("F://tmp//");
    srand(global::mpi().getRank());
    delete_path("F://tmp//");
    const T vIn = 0.005;
    const plint maxIter = 50000;
    const plint saveIter = 10;
    const plint statIter = 10;

    const T omega1 = 1.;

    MultiScalarField2D<bool> boolMask(Nx, Ny);

    plb_ifstream ifile("PPGeo.dat");
    ifile >> boolMask;

    // Use regularized BGK dynamics to improve numerical stability (but note that
    //   BGK dynamics works well too).
    MultiBlockLattice2D<T, DESCRIPTOR> heavyFluid(
        Nx, Ny, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega1));

    // Boundary conditions.
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *
        boundaryCondition = createLocalBoundaryCondition2D<T, DESCRIPTOR>();

    heavyFluid.periodicity().toggle(0,true);

    // Store a pointer to all lattices (two in the present application) in a vector to
    //   create the Shan/Chen coupling therm. The heavy fluid being at the first place
    //   in the vector, the coupling term is going to be executed at the end of the call
    //   to collideAndStream() or stream() for the heavy fluid.*/
    pcout << "Starting simulation" << endl;
    // Main loop over time iterations.
    for (plint iT = 0; iT < maxIter; ++iT)
    {
        rayleighTaylorSetup(heavyFluid, boolMask, boundaryCondition, vIn, iT);
        heavyFluid.collideAndStream();
        if (iT % saveIter == 0)
        {
            writeGifs(heavyFluid, iT);
            writeVTKs(heavyFluid, iT);
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
                  << getStoredAverageDensity<T>(heavyFluid)<<endl;
        }
    }
    delete boundaryCondition;
}
