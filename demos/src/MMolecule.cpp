/*
  Copyright (c) 2015 Lester Hedges <lester.hedges+vmmc@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "MMolecule.hpp"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal.h"
#include "up_down_pass.h"
#include "verify.h"
#include "kernel.h"
#include "thread.h"

using exafmm::Bodies;
using exafmm::BoundBox;
using exafmm::Bounds;
using exafmm::BuildTree;
using exafmm::Cells;
using exafmm::Traversal;
using exafmm::UpDownPass;

MMolecule::MMolecule(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_,
    NodeInfo * nodeinfo_) :
    Model(box_, particles_, cells_, maxInteractions_, interactionEnergy_, interactionRange_)
{
    const exafmm::real_t cycle = 2 * M_PI;
    if(nodeinfo_ != NULL) {
        args.threads = nodeinfo_->threads;
    }
    pbodies = new Bodies(/*put something in here */);
    pbodies2 = new Bodies();
    pjbodies = new Bodies();
    pbuffer = new Bodies();
    pboundBox = new BoundBox(args.nspawn);
    pbounds = new Bounds;
    pbuildTree = new BuildTree(args.ncrit, args.nspawn);
    pcells = new Cells();
    pjcells = new Cells();
    ptraversal = new Traversal(args.nspawn, args.images);
    pupDownPass = new UpDownPass(args.theta, args.useRmax, args.useRopt);
    num_threads(args.threads);

    exafmm::kernel::eps2 = 0.0;
    exafmm::kernel::setup();
    exafmm::logger::verbose = args.verbose;
    exafmm::logger::printTitle("FMM Parameters");
    //args.print(exafmm::logger::stringLength, P);

    pbuffer->reserve(pbodies->size());    

    // Work out the potential shift.
    potentialShift = std::pow(1.0/interactionRange, 12) - std::pow(1/interactionRange, 6);
}

double MMolecule::computePairEnergy(unsigned int particle1,
				    double position1[],
				    double orientation1[],
				    unsigned int particle2,
				    double position2[],
				    double orientation2[])
{
    // Separation vector.
    std::vector<double> sep(box.dimension);

    // Calculate separation.
    for (unsigned int i=0;i<box.dimension;i++)
        sep[i] = position1[i] - position2[i];

    // Enforce minimum image.
    box.minimumImage(sep);

    double normSqd = 0;

    // Calculate squared norm of vector.
    for (unsigned int i=0;i<box.dimension;i++)
        normSqd += sep[i]*sep[i];

    // Particles interact.
    if (normSqd < squaredCutOffDistance)
    {
        double r2Inv = 1.0 / normSqd;
        double r6Inv = r2Inv*r2Inv*r2Inv;
        return 4.0*interactionEnergy*((r6Inv*r6Inv) - r6Inv - potentialShift);
    }
    else return 0;
}

double MMolecule::nonPairwiseCallback(unsigned int particle,
				      double position[],
				      double orientation[])
{
   return 0.0;
}

