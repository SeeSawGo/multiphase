/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "palabos2D.h"
#include "palabos2D.hh"

#include "poiseuille.h"
#include "poiseuille.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D2Q9Descriptor
void writeVTKs(MultiBlockLattice2D<T, DESCRIPTOR> &lattice,plint iT)
{
    T dx = 1.0/50.0;
    T dt = dx/10.0;
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iT, 6), dx);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity_H", dx/dt);
    vtkOut.writeData<float>(*computeDensity(lattice), "density_H", 1.);
}


int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    util::delete_path("./tmp/");
    // Define numeric parameters.
    IncomprFlowParam<T> parameters (
        (T) 1e-2,  // uMax
        (T) 300.,  // Re
        40,        // N
        12.575,        // lx
        12.575         // ly
    );

    plint Nx = parameters.getNx();
    plint Ny = parameters.getNy();
    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
              Nx, Ny,
              new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    MultiScalarField2D<bool> boolMask(Nx, Ny);

    plb_ifstream ifile("PPGeo.dat");
    ifile >> boolMask;

    defineDynamics(lattice, boolMask, new BounceBack<T,DESCRIPTOR>, true);

    Box2D inlet(0, 0, 0, Ny-1);
    Box2D outlet(Nx-1, Nx-1, 0, Ny-1);
    boundaryCondition->addVelocityBoundary0N(inlet, lattice);
    setBoundaryVelocity(lattice, inlet, Array<T,2>(parameters.getLatticeU(),0.));

    boundaryCondition->addVelocityBoundary0P(outlet, lattice, boundary::outflow);
    //setBoundaryVelocity(lattice, outlet, Array<T,2>(1e-2,0.));
    //boundaryCondition->setPressureConditionOnBlockBoundaries(lattice, outlet);
    //setBoundaryDensity(lattice, outlet, (T)0.5);
    //createPoiseuilleBoundaries(lattice, parameters, *boundaryCondition);
    lattice.initialize();

    // Main loop over time iterations.
    for (plint iT=0; iT<800000; ++iT) {
        if (iT%1000==0) {
            pcout << "Writing image at dimensionless time " << iT*parameters.getDeltaT() << endl;
            ImageWriter<T> imageWriter("leeloo");
            imageWriter.writeScaledGif (
                    createFileName("velocity", iT, 6),
                    *computeVelocityNorm(lattice) );
            pcout << computeAverageEnergy(lattice) << endl;
            writeVTKs(lattice, iT);
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
