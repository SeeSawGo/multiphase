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

/** \file
 * Helper functions for domain initialization -- header file.
 */
#include "plbWrapper/block/multiBlockInfo2D.h"
#include "plbWrapper/block/multiBlockInfo2D.hh"

namespace plb {

template class MultiNTensorInfo2D<PRECOMP_T>;


MultiBlockInfo2D::MultiBlockInfo2D(MultiBlock2D const& multiBlock)
    : nx(),
      ny(),
      numBlocks(),
      numAllocatedCells()
{
    getMultiBlockInfo(multiBlock, nx, ny, numBlocks,
                      smallest, largest, numAllocatedCells);
}

plint MultiBlockInfo2D::getNx() const {
    return nx;
}

plint MultiBlockInfo2D::getNy() const {
    return ny;
}

plint MultiBlockInfo2D::getNumBlocks() const {
    return numBlocks;
}

plint MultiBlockInfo2D::getNumAllocatedCells() const {
    return numAllocatedCells;
}

Box2D MultiBlockInfo2D::getSmallestBlock() const {
    return smallest;
}

Box2D MultiBlockInfo2D::getLargestBlock() const {
    return largest;
}

}  // namespace plb
