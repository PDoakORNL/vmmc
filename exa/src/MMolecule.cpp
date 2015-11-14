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
#include "nodeinfo.h"
#include "FloatingPoint.h"

using exafmm::Bodies;
using exafmm::B_iter;
using exafmm::BoundBox;
using exafmm::Bounds;
using exafmm::BuildTree;
using exafmm::Cells;
using exafmm::Traversal;
using exafmm::UpDownPass;

MMolecule::MMolecule(
    Box& box_,
    std::vector<Particle>& atoms_,
    CellList& cells_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_,
    NodeInfo * nodeinfo_) :
  Model(box_, atoms_,
	       cells_,
	       maxInteractions_,
	       interactionEnergy_,
	       interactionRange_)
{
    if(nodeinfo_ != NULL) {
        args.threads = nodeinfo_->threads;
    }
    
    pbodies = new Bodies(particles.size());
    pbodies2 = new Bodies(particles.size());
    pjbodies = new Bodies(particles.size());
    pbuffer = new Bodies(particles.size());
    pboundBox = new BoundBox(args.nspawn);
    pbounds = new Bounds;
    pbuildTree = new BuildTree(args.ncrit, args.nspawn);
    pcells = new Cells;
    pjcells = new Cells;
    ptraversal = new Traversal(args.nspawn, args.images);
    pupDownPass = new UpDownPass(args.theta, args.useRmax, args.useRopt);
    num_threads(args.threads);

    exafmm::kernel::eps2 = 0.0;
    exafmm::kernel::setup();
    exafmm::logger::verbose = args.verbose;
    exafmm::logger::printTitle("FMM Parameters");
    //args.print(exafmm::logger::stringLength, P);

    // Work out the potential shift.
    potentialShift = std::pow(1.0/interactionRange, 12) - std::pow(1/interactionRange, 6);

}

void MMolecule::initBodies()
{
    for (unsigned int i = 0; i < particles.size(); i++) {
	(*pbodies)[i].SRC = particles[i].charge;
	for (int j = 0; j < 3; j++) {
	    (*pbodies)[i].X[j] = particles[i].position[j];
	    //std::cout << (*pbodies)[i].X[j] << " " << particles[i].position[j] << std::endl;
	}
    }
}

void MMolecule::updateBodies()
{
    for (unsigned int i = 0; i < particles.size(); i++) {
	for (int j = 0; j < 3; j++) {
	    (*pbodies)[i].X[j] = particles[i].position[j];
	}
    }
    // std::cout << "updated bodies" << std::endl;
}

double MMolecule::computeEnergy(unsigned int particle, double position[], double orientation[])
{
      double energy = 0;      // energy counter
      energy = Model::computeEnergy(particle,position,orientation);
      //std::cout << "energy pair:" << energy;
      double nonPair;
      nonPair = this->nonPairwiseCallback(particle,position,orientation);  
      energy += nonPair;
      //std::cout << " + es " << std::setprecision(16) << nonPair << std::endl;
      return energy;
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
  // std::cout << std::setprecision(16) <<position[0] << " " << position[1] << " " << position[2] << std::endl;
  // std::cout << std::setprecision(16) << (*pbodies)[particle].X[0] << " " << (*pbodies)[particle].X[1]
  // 	    << " " << (*pbodies)[particle].X[2] << std::endl;
  // std::cout << "accepted" << accepted << std::endl;
  FloatingPoint<double> lhs0(position[0]);
  FloatingPoint<double> rhs0((*pbodies)[particle].X[0]);
  FloatingPoint<double> lhs1(position[1]);
  FloatingPoint<double> rhs1((*pbodies)[particle].X[1]);
  FloatingPoint<double> lhs2(position[2]);
  FloatingPoint<double> rhs2((*pbodies)[particle].X[2]);
  
  if (0) {// ! lhs0.AlmostEquals(rhs0) ||
      // ! lhs1.AlmostEquals(rhs1) ||
      // ! lhs2.AlmostEquals(rhs2)) {
    //std::cout << "clean" << accepted << std::endl;
	pjbodies->push_back((*pbodies)[particle]);
	//pjbodies->back().TRG[0] = pjbodies->back().SRC;
	//pjbodies->back().TRG[1] = pjbodies->back().X[0];
	//pjbodies->back().TRG[2] = pjbodies->back().X[1];
	//pjbodies->back().TRG[3] = pjbodies->back().X[2];
	//std::cout << pjbodies->back().SRC << " " << pjbodies->back().TRG[0] << std::endl;
	initTargets(pjbodies);
	*pbounds = pboundBox->getBounds(*pbodies);
	*pbounds = pboundBox->getBounds(*pjbodies,*pbounds);
	*pjcells = pbuildTree->buildTree(*pjbodies, *pbuffer, *pbounds);
	pupDownPass->upwardPass(*pjcells);
	ptraversal->traverse(*pcells, *pjcells, cycle, args.dual, false);
	pupDownPass->downwardPass(*pcells);
	double result = 0;
	result = (pjbodies->back().TRG[0]) * (pjbodies->back().SRC);
	//std::cout << result << std::endl;
	if (result!=result) {
	  std::cout << particle << " produced: " << result << std::endl;
	}
	return result * c_ftov;
  } else if (1) { //(accepted) {
    //std::cout << "rebuild main tree" << std::endl;
        initTargets(pbodies);
        *pbounds = pboundBox->getBounds(*pbodies);
	*pcells = pbuildTree->buildTree(*pbodies, *pbuffer, *pbounds);
	//pjbodies->clear();
	pupDownPass->upwardPass(*pcells);
	ptraversal->initListCount(*pcells);
	ptraversal->initWeight(*pcells);
	ptraversal->traverse(*pcells, *pcells, cycle, args.dual, args.mutual);
	pupDownPass->downwardPass(*pcells);
	ptraversal->writeList(*pcells,0);
	//ptraversal->normalize(*pbodies);
	accepted=0;
	double result = 0;
	result = ((*pbodies)[particle].TRG[0] * (*pbodies)[particle].SRC);
	if (result!=result) {
	  std::cout << particle << " produced: " << result << std::endl;
	}
	return result * c_ftov;
    } else {
    return ((*pbodies)[particle].TRG[0] * (*pbodies)[particle].SRC) * c_ftov;
	std::cout << "unexpected non-pairwise result" << std::endl;
    }
    
}

void MMolecule::initTargets(Bodies * ptbodies) {
  for (B_iter B=ptbodies->begin(); B!=ptbodies->end(); B++) {     // Loop over bodies
	  B->TRG = 0;                                             //  Clear target values
	  B->IBODY = B-ptbodies->begin();                            //  Initial body numbering
	  B->WEIGHT = 1;                                          //  Initial weight
	}                                                         // End loop over bodies
}

void MMolecule::applyPostMoveUpdates(unsigned int particle, double position[], double orientation[])
{
  Model::applyPostMoveUpdates(particle, position, orientation);
  for (unsigned int i=0;i<box.dimension;i++)
  {
      (*pbodies)[particle].X[i] = particles[particle].position[i];
  }
  accepted=1;
}
  
