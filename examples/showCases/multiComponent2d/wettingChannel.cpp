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
#include "poiseuille.h"
#include "poiseuille.hh"
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

plint Nx = 400;
plint Ny = 100;
plint Nw = 100;
plint w1 = 19; // inlet1
plint w2 = 48; // intlet2
T rhoH = 1.0;
T rhoL = 0.8;
plint initialH = 60;
Box2D mainPipe(1, Nx-2, 1, w2);
Box2D verticalPipe1(Nw, Nw+w1-1, w2+1, initialH); //heavy
Box2D verticalPipe2(Nw, Nw+w1-1, initialH+1, Ny-2);//light
Box2D inlet1(Nw, Nw+w1-1, Ny-1, Ny-1); //light
Box2D inlet2(0,0,1,w2);//heavy
Box2D outlet(Nx-1,Nx-1, 1, w2);

void rayleighTaylorSetup(MultiBlockLattice2D<T, DESCRIPTOR> &heavyFluid,
                         MultiBlockLattice2D<T, DESCRIPTOR> &lightFluid,
                         MultiScalarField2D<bool> &boolMask,
                         OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition,
                         IncomprFlowParam<T>& parameter_mainpipe, IncomprFlowParam<T>& parameter_verticalpipe,
                         T rho0, T rho1, T v1, T v2,
                         plint iT)
{
    Array<T,2> zeroVelocity((T)0.,(T)0.);
    T zeroDensity = (T)1.e-10;
    // Bounce-back
    defineDynamics(heavyFluid, boolMask, new BounceBack<T, DESCRIPTOR>(rho0), true);
    defineDynamics(lightFluid, boolMask, new BounceBack<T, DESCRIPTOR>(rho1), true);
    if (iT==0) ///initial
    {

    ///heavy:
    initializeAtEquilibrium(heavyFluid, mainPipe, rhoH, zeroVelocity);
    initializeAtEquilibrium(heavyFluid, verticalPipe1, rhoH, zeroVelocity);
    initializeAtEquilibrium(heavyFluid, verticalPipe2, zeroDensity, zeroVelocity);
    //initializeAtEquilibrium(heavyFluid, inlet1, zeroDensity, zeroVelocity);
    initializeAtEquilibrium(heavyFluid, inlet2, rhoH, zeroVelocity);
    initializeAtEquilibrium(heavyFluid, outlet, zeroDensity, zeroVelocity);

    ///light:
    initializeAtEquilibrium(lightFluid, mainPipe, zeroDensity, zeroVelocity);
    initializeAtEquilibrium(lightFluid, verticalPipe1, zeroDensity, zeroVelocity);
    initializeAtEquilibrium(lightFluid, verticalPipe2, rhoL, zeroVelocity);
    //initializeAtEquilibrium(lightFluid, inlet1, rhoL, zeroVelocity);
    initializeAtEquilibrium(lightFluid, inlet2, zeroDensity, zeroVelocity);
    initializeAtEquilibrium(lightFluid, outlet, zeroDensity, zeroVelocity);

    boundaryCondition->addVelocityBoundary1P(inlet1, lightFluid);
    setBoundaryVelocity(lightFluid, inlet1, Array<T, 2>((T)0., -v1));
    initializeAtEquilibrium(lightFluid, inlet1, rhoL, Array<T, 2>((T)0., -v1));

    heavyFluid.initialize();
    lightFluid.initialize();

    }
    if(iT==0)
    {
    ///heavy:
    //boundaryCondition->setPressureConditionOnBlockBoundaries(heavyFluid, inlet1);
    //setBoundaryDensity(heavyFluid, inlet1, zeroDensity);
    //boundaryCondition->addVelocityBoundary1P(inlet1, heavyFluid, boundary::outflow);
    //initializeAtEquilibrium(heavyFluid, inlet1, zeroDensity, zeroVelocity);
    //initializeAtEquilibrium(heavyFluid, verticalPipe2, zeroDensity, zeroVelocity);

    //boundaryCondition->setPressureConditionOnBlockBoundaries(heavyFluid, inlet2);
    //setBoundaryDensity(heavyFluid, inlet2, 2*rhoH);
    //initializeAtEquilibrium(heavyFluid, inlet2, 2*rhoH, zeroVelocity);
    boundaryCondition->addVelocityBoundary0N(Box2D(0,0,0,w2+2), heavyFluid);
    setBoundaryVelocity(heavyFluid, Box2D(0,0,0,w2+2), PoiseuilleVelocity<T>(parameter_mainpipe));
    initializeAtEquilibrium(heavyFluid, Box2D(0,0,0,w2+2), PoiseuilleVelocityAndDensity<T,DESCRIPTOR>(parameter_mainpipe));

    //boundaryCondition->addVelocityBoundary0P(outlet, heavyFluid);
    //setBoundaryVelocity(heavyFluid, outlet, Array<T,2>((T)0., (T)0.));
    boundaryCondition->addVelocityBoundary0P(outlet, heavyFluid, boundary::outflow);
    //boundaryCondition->setPressureConditionOnBlockBoundaries(heavyFluid, outlet);
    //setBoundaryDensity(heavyFluid, outlet, zeroDensity);
    //initializeAtEquilibrium(heavyFluid, outlet, zeroDensity, zeroVelocity);

    ///light:
    //boundaryCondition->addVelocityBoundary1P(inlet1, lightFluid);
    //setBoundaryVelocity(lightFluid, inlet1, Array<T, 2>((T)0., -v1));
    
    //initializeAtEquilibrium(lightFluid, inlet1, rhoL, Array<T, 2>((T)0., -v1));

    //boundaryCondition->setPressureConditionOnBlockBoundaries(lightFluid, inlet2);
    //setBoundaryDensity(lightFluid, inlet2, zeroDensity);
    boundaryCondition->addVelocityBoundary0N(inlet2, lightFluid, boundary::outflow);
    //initializeAtEquilibrium(lightFluid, inlet2, zeroDensity, zeroVelocity);

    //boundaryCondition->setPressureConditionOnBlockBoundaries(lightFluid, outlet);
    //setBoundaryDensity(lightFluid, outlet, zeroDensity);
    boundaryCondition->addVelocityBoundary0P(outlet, lightFluid, boundary::outflow);
    //initializeAtEquilibrium(lightFluid, outlet, zeroDensity, zeroVelocity);
    }
    if(iT>100000)
    {
        boundaryCondition->addVelocityBoundary1P(Box2D(Nw-1, Nw+w1, Ny-1, Ny-1), lightFluid);
        initializeAtEquilibrium(heavyFluid, verticalPipe2, zeroDensity, zeroVelocity);
        setBoundaryVelocity(lightFluid, Box2D(Nw-1, Nw+w1, Ny-1, Ny-1), PoiseuilleVelocity<T>(parameter_verticalpipe));
        initializeAtEquilibrium(lightFluid, Box2D(Nw-1, Nw+w1, Ny-1, Ny-1), PoiseuilleVelocityAndDensity<T,DESCRIPTOR>(parameter_verticalpipe));

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
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iT, 6), 1.);
    vtkOut.writeData<2,float>(*computeVelocity(heavyFluid), "velocity_H", 1.);
    vtkOut.writeData<float>(*computeDensity(heavyFluid), "density_H", 1.);
    vtkOut.writeData<2,float>(*computeVelocity(lightFluid), "velocity_L", 1.);
    vtkOut.writeData<float>(*computeDensity(lightFluid), "density_L", 1.);
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
    global::directories().setOutputDir("./tmp/");
    srand(global::mpi().getRank());
    delete_path("./tmp/");
    const T v1 = 0.001;
    const T v2 = 0.005;
    const T G = 3.;
    T rho0 = 0.7;
    T rho1 = (T)1.0-rho0;
    IncomprFlowParam<T> parameter_mainpipe(
            (T) v2,  // uMax
            (T) 300.,  // Re
            w2+2,       // N
            5.,        // lx
            1.         // ly 
    );
    IncomprFlowParam<T> parameter_verticalpipe(
            (T) -v1,  // uMax
            (T) 300.,  // Re
            w1+2,       // N
            5.,        // lx
            1.         // ly 
    );
    const plint maxIter = 50000;
    const plint saveIter = 100;
    const plint statIter = 100;

    const T omega1 = 1.;
    const T omega2 = 1.;

    MultiScalarField2D<bool> boolMask(Nx, Ny);

    plb_ifstream ifile("Geo.dat");
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

    //heavyFluid.periodicity().toggle(0,true);
    //lightFluid.periodicity().toggle(0,true);

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
            Box2D(0, Nx-2, 0, Ny-1), blockLattices, processorLevel);//contain boundary or not? why?
    /*    
    integrateProcessingFunctional(
        new ShanChenMultiComponentProcessor2D<T,DESCRIPTOR> (G,constOmegaValues),
            inlet1, blockLattices, processorLevel);
    integrateProcessingFunctional(
        new ShanChenMultiComponentProcessor2D<T,DESCRIPTOR> (G,constOmegaValues),
            inlet2, blockLattices, processorLevel);
        /*
    integrateProcessingFunctional(
        new ShanChenMultiComponentProcessor2D<T,DESCRIPTOR> (G,constOmegaValues),
            outlet, blockLattices, processorLevel);         
    integrateProcessingFunctional(
        new ShanChenMultiComponentProcessor2D<T,DESCRIPTOR> (G,constOmegaValues),
            mainPipe, blockLattices, processorLevel);
    integrateProcessingFunctional(
        new ShanChenMultiComponentProcessor2D<T,DESCRIPTOR> (G,constOmegaValues),
            verticalPipe1, blockLattices, processorLevel);
    integrateProcessingFunctional(
        new ShanChenMultiComponentProcessor2D<T,DESCRIPTOR> (G,constOmegaValues),
            verticalPipe2, blockLattices, processorLevel);
        */    

    

    //SetBCs(heavyFluid,lightFluid,boundaryCondition_heavy,boundaryCondition_light, v1, v2);
    pcout << "Starting simulation" << endl;
    // Main loop over time iterations.
    for (plint iT = 0; iT < maxIter; ++iT)
    {
        rayleighTaylorSetup(heavyFluid, lightFluid, boolMask, boundaryCondition,parameter_mainpipe, parameter_verticalpipe,rho0, rho1, v1, v2, iT);
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
