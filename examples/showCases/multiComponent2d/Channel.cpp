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
//#include "poiseuille.h"
//#include "poiseuille.hh"
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
plint w1 = 50; // inlet1
plint w2 = 50; // inlet2

void rayleighTaylorSetup( MultiBlockLattice2D<T, DESCRIPTOR>& heavyFluid,
                          T rho0,
                          MultiScalarField2D<bool>& boolMask)
{
    // Bounce-back
    defineDynamics(heavyFluid, boolMask, new BounceBack<T, DESCRIPTOR>, true );

    heavyFluid.initialize();
}

void writeGifs(MultiBlockLattice2D<T, DESCRIPTOR>& heavyFluid,
               plint iT)
{
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("rho_heavy_", iT, 6),
                               *computeVelocityNorm(heavyFluid));
}

void delete_path(const char* path){  
    
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
                if (tmp.find(".") == -1){  
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
    const T rho0 = 0.0;
    const T v1     =0.005;
    const T v2     =0.01;
    const plint maxIter  = 16000;
    const plint saveIter = 100;
    const plint statIter = 100;

    // Liquid one parameters
    IncomprFlowParam<T> parameters1 (
        v1,  // uMax
        (T) 10.,  // Re
        5,        // N
        80,        // lx
        20.         // ly
    );
    const T omega1 = 1;//parameters1.getOmega();
     pcout << "omega1 = "
                  << omega1<<endl;

    MultiScalarField2D<bool> boolMask(Nx, Ny);

    plb_ifstream ifile("Geo.dat");
    ifile >> boolMask;

    // Use regularized BGK dynamics to improve numerical stability (but note that
    //   BGK dynamics works well too).
    MultiBlockLattice2D<T, DESCRIPTOR> heavyFluid (
        Nx,Ny, new CompleteBGKdynamics<T, DESCRIPTOR>(omega1) );

    // Boundary conditions.
    
    // Heavy Fluid
    
    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition1 = createLocalBoundaryCondition2D<T,DESCRIPTOR>();
    //vertical Pipe inlet:
    Box2D inlet1H(Nw, Nw+w1-1, Ny-1, Ny-1);
    boundaryCondition1->addPressureBoundary1P(inlet1H, heavyFluid);
    setBoundaryDensity(heavyFluid, inlet1H, 1.5);
    //boundaryCondition1->addVelocityBoundary1P(inlet1H, heavyFluid);
    //boundaryCondition1->addExternalVelocityCornerNP(Nw-1,Ny-1, heavyFluid);
    //boundaryCondition1->addExternalVelocityCornerPP(Nw+w1,Ny-1, heavyFluid);
    //Box2D topH(Nw-1, Nw+w1, Ny-1, Ny-1);
    //setBoundaryVelocity(heavyFluid, topH, Array<T,2>(-v1,0.));
    //horizontal pipe inlet:
    Box2D inlet2H(0, 0, 1, w2-2);
    boundaryCondition1->addVelocityBoundary0N(inlet2H, heavyFluid);
    boundaryCondition1->addExternalVelocityCornerNP(0,0, heavyFluid);
    boundaryCondition1->addExternalVelocityCornerNN(0,w2-1, heavyFluid);
    Box2D leftH(0,0, 0,w2-1);
    setBoundaryVelocity(heavyFluid, leftH, Array<T,2>(v2,0.));
    //horizontal pipe outletH:
    Box2D outletH(Nx-1,Nx-1, 1,w2-2);
    boundaryCondition1->addVelocityBoundary0P(outletH, heavyFluid, boundary::outflow);
    boundaryCondition1->addExternalVelocityCornerPP(Nx-1, 0, heavyFluid, boundary::outflow);
    boundaryCondition1->addExternalVelocityCornerPN(Nx-1, w2-1, heavyFluid, boundary::outflow);
    

    /*
    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition1 = createLocalBoundaryCondition2D<T,DESCRIPTOR>();
    //vertical Pipe inlet:
    Box2D inlet1(Nw, Nw+49, Ny-1, Ny-1);
    boundaryCondition1->addVelocityBoundary1P(inlet1, heavyFluid);
    //boundaryCondition1->addExternalVelocityCornerNP(Nw-1,Ny-1, heavyFluid);
    //boundaryCondition1->addExternalVelocityCornerPP(Nw+50,Ny-1, heavyFluid);
    //Box2D top(Nw-1, Nw+50, Ny-1, Ny-1);
    setBoundaryVelocity(heavyFluid, inlet1, Array<T,2>(0., -v1));
    //horizontal pipe inlet:
    
    Box2D inlet2(0, 0, 1, 48);
    boundaryCondition1->addVelocityBoundary0N(inlet2, heavyFluid);
    boundaryCondition1->addExternalVelocityCornerNP(0,0, heavyFluid);
    boundaryCondition1->addExternalVelocityCornerNN(0,49, heavyFluid);
    Box2D left(0,0, 0,49);
    setBoundaryVelocity(heavyFluid, left, Array<T,2>(v2, 0.));
    
    //horizontal pipe outlet:
    Box2D outlet(Nx-1,Nx-1, 1,48);
    boundaryCondition1->addVelocityBoundary0P(outlet, heavyFluid, boundary::outflow);
    boundaryCondition1->addExternalVelocityCornerPP(Nx-1, 0, heavyFluid, boundary::outflow);
    boundaryCondition1->addExternalVelocityCornerPN(Nx-1, 49, heavyFluid, boundary::outflow);
    */

    rayleighTaylorSetup(heavyFluid, rho0, boolMask);

    pcout << "Starting simulation" << endl;
    // Main loop over time iterations.
    for (plint iT=0; iT<maxIter; ++iT)
    {
        if (iT%saveIter==0)
        {
            writeGifs(heavyFluid,iT);
        }

        heavyFluid.collideAndStream();

        if (iT%statIter==0)
        {
            pcout << "Average density fluid one = "
                  << getStoredAverageEnergy<T>(heavyFluid)<<endl;
        }
    }
}
