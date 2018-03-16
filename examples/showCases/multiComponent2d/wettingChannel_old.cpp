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
plint w1 = 15; // inlet1
plint w2 = 50; // intlet2
/// Initial condition: heavy fluid on top, light fluid on bottom.
/** This functional is going to be used as an argument to the function "applyIndexed",
 *  to setup the initial condition. For efficiency reasons, this approach should
 *  always be preferred over explicit space loops in end-user codes.
 */
template <typename T, template <typename U> class Descriptor>
class TwoLayerInitializer : public OneCellIndexedWithRandFunctional2D<T, Descriptor>
{
  public:
    TwoLayerInitializer(bool isLiquid1_, plint ny_)
        : isLiquid1(isLiquid1_), ny(ny_)
    {
    }
    TwoLayerInitializer<T, Descriptor> *clone() const
    {
        return new TwoLayerInitializer<T, Descriptor>(*this);
    }
    virtual void execute(plint iX, plint iY, T rand_val, Cell<T, Descriptor> &cell) const
    {
        T densityFluctuations = 0;
        T almostNoFluid = 1e-4;
        T rho = 1.;
        Array<T, 2> zeroVelocity(0., 0.);
        if (isLiquid1 && (iY < ny))
        {
            rho = 1 + rand_val * densityFluctuations;
        }
        else if ((!isLiquid1) && iY >=ny)
        {
            rho = 0.8 + rand_val * densityFluctuations;
        }
        else
        {
            rho = almostNoFluid;
        }
        iniCellAtEquilibrium(cell, rho, zeroVelocity);
    }

  private:
    bool isLiquid1;
    plint ny;
};

void rayleighTaylorSetup(MultiBlockLattice2D<T, DESCRIPTOR> &heavyFluid,
                         MultiBlockLattice2D<T, DESCRIPTOR> &lightFluid,
                         MultiScalarField2D<bool> &boolMask,
                         T rho0, T rho1,
                         T force)
{
    // The setup is: periodicity along horizontal direction, bounce-back on top
    // and bottom. The upper half is initially filled with fluid 1 + random noise,
    // and the lower half with fluid 2. Only fluid 1 experiences a forces,
    // directed downwards.
    plint nx = heavyFluid.getNx();
    plint ny = heavyFluid.getNy();

    // Bounce-back
    defineDynamics(heavyFluid, boolMask, new BounceBack<T, DESCRIPTOR>(rho0), true);
    defineDynamics(lightFluid, boolMask, new BounceBack<T, DESCRIPTOR>(rho1), true);

    // Initialize Liquid one.
    applyIndexed(heavyFluid, Box2D(0, nx-1, 1, ny-1),
                 new TwoLayerInitializer<T, DESCRIPTOR>(true, 200));
    // Initialize Liquid two.
    applyIndexed(lightFluid, Box2D(1, nx-2, 1, ny-1),
                 new TwoLayerInitializer<T, DESCRIPTOR>(false, 200));

    // Let's have gravity acting on the heavy fluid only. This represents a situation
    //   where the molecular mass of the light species is very small, and thus the
    //   action of gravity on this species is negligible.
    
    setExternalVector(heavyFluid, heavyFluid.getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,2>(0.,0.));
    setExternalVector(lightFluid, lightFluid.getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,2>(0.,-force));
}
void SetBCs(MultiBlockLattice2D<T, DESCRIPTOR> &heavyFluid,
            MultiBlockLattice2D<T, DESCRIPTOR> &lightFluid,
            OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition_heavy,
            OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition_light,
            T v1, T v2)
{
    
    //vertical Pipe inlet:
    /*
    Box2D inlet1H(Nw, Nw+w1-2, Ny-1, Ny-1);
    boundaryCondition_heavy->setPressureConditionOnBlockBoundaries(heavyFluid, inlet1H);
    setBoundaryDensity(heavyFluid, inlet1H, 2.);
    */
    //horizontal pipe inlet:
    Box2D inlet2H(1, 1, 1, w2-2);
    //boundaryCondition_heavy->setPressureConditionOnBlockBoundaries(heavyFluid, inlet2H);
    //setBoundaryDensity(heavyFluid, inlet2H, 10.);
    boundaryCondition_heavy->addVelocityBoundary0N(inlet2H, heavyFluid);
    setBoundaryVelocity(heavyFluid, inlet2H, Array<T,2>(v2, 0.));
    //horizontal pipe outletH:
    Box2D outletH(Nx, Nx, 1, w2-2);
    boundaryCondition_heavy->addVelocityBoundary0P(outletH, heavyFluid, boundary::outflow);
    //setBoundaryVelocity(heavyFluid, outletH, Array<T,2>(v2, 0.));
    //boundaryCondition_heavy->setPressureConditionOnBlockBoundaries(heavyFluid, outletH);
    //setBoundaryDensity(heavyFluid, outletH, 0.0);

    
    //vertical Pipe inlet:
    
    Box2D inlet1L(Nw, Nw+w1-2, Ny-1, Ny-1);
    boundaryCondition_light->addVelocityBoundary1P(inlet1L, lightFluid);
    setBoundaryVelocity(lightFluid, inlet1L, Array<T, 2>(0., -v1));
    //boundaryCondition_light->setPressureConditionOnBlockBoundaries(lightFluid, inlet1L);
    //setBoundaryDensity(lightFluid, inlet1L, 1.6);
    //horizontal pipe inlet:
    
    //Box2D inlet2L(0, 0, 1, w2-2);
    //boundaryCondition_light->addVelocityBoundary0N(inlet2L, lightFluid);
    //setBoundaryVelocity(lightFluid, inlet2L, Array<T, 2>(v1, 0.));
    //boundaryCondition_light->setPressureConditionOnBlockBoundaries(lightFluid, inlet2L);
    //setBoundaryDensity(lightFluid, inlet2L, 1.);
    //horizontal pipe outletL:
    //Box2D outletL(Nx-1, Nx-1, 1, w2-2);
    //boundaryCondition_light->addVelocityBoundary0P(outletL, lightFluid, boundary::outflow);
    //setBoundaryDensity(lightFluid, outletL, 0.);
    
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
    const T v1 = 0.05;
    const T v2 = 0.05;
    const T G = 3.2;
    T rho0 = 0.6;
    T rho1 = 1-rho0; 
    T force = 0.0;//0.5/(T)Ny;
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

    //heavyFluid.periodicity().toggle(0,true);
    lightFluid.periodicity().toggle(0,true);
    // Boundary conditions.
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *
        boundaryCondition_heavy = createLocalBoundaryCondition2D<T, DESCRIPTOR>();
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *
        boundaryCondition_light = createLocalBoundaryCondition2D<T, DESCRIPTOR>();

   
   
    // Store a pointer to all lattices (two in the present application) in a vector to
    //   create the Shan/Chen coupling therm. The heavy fluid being at the first place
    //   in the vector, the coupling term is going to be executed at the end of the call
    //   to collideAndStream() or stream() for the heavy fluid.*/
    vector<MultiBlockLattice2D<T, DESCRIPTOR> *> blockLattices;
    blockLattices.push_back(&heavyFluid);
    blockLattices.push_back(&lightFluid);

    // The argument "constOmegaValues" to the Shan/Chen processor is optional,
    //   and is used for efficiency reasons only. It tells the data processor
    //   that the relaxation times are constant, and that their inverse must be
    //   computed only once.
    std::vector<T> constOmegaValues;
    constOmegaValues.push_back(omega1);
    constOmegaValues.push_back(omega2);
    plint processorLevel = 1;
    integrateProcessingFunctional(
        new ShanChenMultiComponentProcessor2D<T,DESCRIPTOR> (G,constOmegaValues),
            Box2D(1, Nx+10, 1, Ny-2), blockLattices, processorLevel);

    rayleighTaylorSetup(heavyFluid, lightFluid, boolMask, rho0, rho1, force);
    SetBCs(heavyFluid,lightFluid,boundaryCondition_heavy,boundaryCondition_light, v1, v2);

    lightFluid.initialize();
    heavyFluid.initialize();
    
    pcout << "Starting simulation" << endl;
    // Main loop over time iterations.
    for (plint iT = 0; iT < maxIter; ++iT)
    {
        if (iT % saveIter == 0)
        {
            writeGifs(heavyFluid, lightFluid, iT);
            writeVTKs(heavyFluid, lightFluid, iT);
        }
        //SetBCs(heavyFluid,lightFluid,boundaryCondition_heavy,boundaryCondition_light, v1, v2);
        // Time iteration for the light fluid.
        lightFluid.collideAndStream();
        heavyFluid.collideAndStream();
        
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
    delete boundaryCondition_heavy;
    delete boundaryCondition_light;
}
