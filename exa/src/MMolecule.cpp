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
#include <omp.h>

using exafmm::Bodies;
using exafmm::B_iter;
using exafmm::BoundBox;
using exafmm::Bounds;
using exafmm::BuildTree;
using exafmm::Cells;
using exafmm::Traversal;
using exafmm::UpDownPass;
using exafmm::vec3;

MMolecule::MMolecule(
    Box& box_,
    Molecules& molecules_,
    CellList& cells_,
    double particleDiameter_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_,
    NodeInfo * nodeinfo_) : Model(box_, (Particles&) molecules_, cells_ , maxInteractions_, interactionEnergy_, interactionRange_), molecules(molecules_),particleDiameter(particleDiameter_)
{
    if(nodeinfo_ != NULL) {
        args.threads = nodeinfo_->threads;
    }
    //args.threads = 8;
 
    pbodies = NULL;
    pvbodies = NULL;
    pparticles = NULL;
    // pbodies2 = new Bodies(particles.size());
    // pjbodies = new Bodies(particles.size());
    // pbuffer = NULL;
    pboundBox = new BoundBox(args.nspawn);
    pbounds = new Bounds;
    pbuildTree = new BuildTree(args.ncrit, args.nspawn);
    pcells = new Cells;
    pjcells = new Cells;
    ptraversal = new Traversal(args.nspawn, args.images);
    pupDownPass = new UpDownPass(args.theta, args.useRmax, args.useRopt);
    pvboundBox = new BoundBox(args.nspawn);
    pvbounds = new Bounds;
    pvbuildTree = new BuildTree(args.ncrit, args.nspawn);
    pvcells = new Cells;
    pvjcells = new Cells;
    pvtraversal = new Traversal(args.nspawn, args.images);
    pvupDownPass = new UpDownPass(args.theta, args.useRmax, args.useRopt);
    
    exafmm::kernel::eps2 = 0.0;
    exafmm::kernel::setup();
    exafmm::logger::verbose = args.verbose;
    exafmm::logger::printTitle("FMM Parameters");

    //args.print(exafmm::logger::stringLength, P);

    // Work out the potential shift.
    potentialShift = std::pow(1.0/interactionRange, 12) - std::pow(1/interactionRange, 6);
}

MMolecule::~MMolecule() {
    delete pbodies;
    delete pparticles;
    delete pboundBox;
    delete pbounds;
    delete pbuildTree;
    delete pcells;
    delete ptraversal;
    delete pupDownPass;
}    

void MMolecule::initBodies()
{
    long unsigned int numAtoms = molecules.size()*molecules[0].atoms.size();

    std::cout << "pre news" << std::endl;
    pbodies = new Bodies(numAtoms);
    pvbodies = new Bodies(numAtoms);
    pbuffer = new Bodies(numAtoms);
    pparticles = new Bodies(molecules.size());
    std::cout << "before doubles" << std::endl;
    unsigned int nat = args.nat;
    args.prscale = new double[nat*nat];
    args.pgscale = new double[nat*nat];
    args.pfgscale = new double[nat*nat];
    for (unsigned int i = 0; i < nat*nat; i++) { 
	args.prscale[i] = 1;
	args.pgscale[i] = 0.0001;
	args.pfgscale[i] = 0.0001;
    }
    std::cout << "pre VanDerWaals" << std::endl;
    VDW = new VanDerWaals(1.0, 5.0, cycle, nat,
			  args.prscale, args.pgscale, args.pfgscale);
    // pbodies2 = new Bodies(particles.size());
    // pjbodies = new Bodies(particles.size());
    //pbuffer = new Bodies(numAtoms);
    std::vector<double> apos(box.dimension);
    unsigned int bc = 0;
    for (unsigned int i = 0; i < molecules.size(); i++) {
	for (unsigned int j = 0; j < box.dimension; j++)
	{
	    (*pparticles)[i].X[j] = molecules[i].position[j];
	}
    	for(Molecule::atom_iter a=molecules[i].atoms.begin();
	    a < molecules[i].atoms.end();
	    a++)
	{
	    molecules[i].get_apos(*a, molecules.particleDiameter, apos);
	    (*pbodies)[bc].SRC = a->charge;
	    for (unsigned int j = 0; j < 3; j++) {
		if (j < box.dimension) {
		    (*pbodies)[bc].X[j] = apos[j];
		} else {
		    (*pbodies)[bc].X[j] = 0;
		}
	    }
	    bc++;
	}
    }
    std::cout << "pre pvbodies copy" <<std::endl;
    *pvbodies = *pbodies; 
}

void MMolecule::updateBodies()
{
    std::vector<double> apos(box.dimension);
    unsigned int bc = 0;
    for (unsigned int i = 0; i < molecules.size(); i++) {
	for (unsigned int j = 0; j < box.dimension; j++)
	{
	    (*pparticles)[i].X[j] = molecules[i].position[j];
	}
    	for(Molecule::atom_iter a=molecules[i].atoms.begin();
	    a < molecules[i].atoms.end();
	    a++)
	{
	    molecules[i].get_apos(*a, molecules.particleDiameter, apos);
	    for (unsigned int j = 0; j < 3; j++) {
		if (j < box.dimension) {
		    (*pbodies)[bc].X[j] = apos[j];
		} else {
		    (*pbodies)[bc].X[j] = 0;
		}
	    }
	    bc++;
	}
    }
    *pvbodies = *pbodies;
}

double MMolecule::computeEnergy(unsigned int particle, double position[], double orientation[])
{
    double energy = 0;      // energy counter
    unsigned int atomIndex = 0;
    FloatingPoint<double> lhs0(position[0]);
    FloatingPoint<double> rhs0((*pparticles)[particle].X[0]);
    FloatingPoint<double> lhs1(position[1]);
    FloatingPoint<double> rhs1((*pparticles)[particle].X[1]);
    // FloatingPoint<double> lhs2(position[2]);
    // FloatingPoint<double> rhs2((*pparticles)[particle].X[2]);
    if (!moved && lhs0.AlmostEquals(rhs0) &&
	lhs1.AlmostEquals(rhs1) ){ //&&
//	lhs2.AlmostEquals(rhs2)) {
	
	for (unsigned int i = 0; i < molecules.size(); i++)
	{
	    for (unsigned int j = 0; j < box.dimension; j++)
	    {
		(*pparticles)[i].X[j] = molecules[i].position[j];
	    }
	    for(Molecule::atom_iter a=molecules[i].atoms.begin();
		a < molecules[i].atoms.end();
		a++)
	    {
		energy += ((*pbodies)[atomIndex].TRG[0] * (*pbodies)[atomIndex].SRC);
		atomIndex++; 
	    }
	}

	atomIndex = 0;
	for (unsigned int i = 0; i < molecules.size(); i++)
	{
	    for (unsigned int j = 0; j < box.dimension; j++)
	    {
		(*pparticles)[i].X[j] = molecules[i].position[j];
	    }
	    for(Molecule::atom_iter a=molecules[i].atoms.begin();
		a < molecules[i].atoms.end();
		a++)
	    {
		energy += (*pvbodies)[atomIndex].TRG[0];
		atomIndex++; 
	    }
	}
	
	energy = energy * c_ftov;

	if (energy!=energy) {
	    std::cout << particle << " produced: " << energy << std::endl;
	}

    } else {
        initTargets(pbodies);
        *pbounds = pboundBox->getBounds(*pbodies);
	*pcells = pbuildTree->buildTree(*pbodies, *pbuffer, *pbounds);
	pupDownPass->upwardPass(*pcells);
	ptraversal->initListCount(*pcells);
	ptraversal->initWeight(*pcells);
	ptraversal->traverse(*pcells, *pcells, cycle, args.dual, args.mutual);
	pupDownPass->downwardPass(*pcells);
	ptraversal->writeList(*pcells,0);
	vec3 dipole = pupDownPass->getDipole(*pbodies,0);
	// vec3 localDipole = upDownPass.getDipole(bodies,0);
	// vec3 globalDipole = baseMPI.allreduceVec3(localDipole);
	pupDownPass->dipoleCorrection(*pbodies,dipole,(*pbodies).size(), cycle);

	moved=0;

	for (unsigned int i = 0; i < molecules.size(); i++)
	{
	    for (unsigned int j = 0; j < box.dimension; j++)
	    {
		(*pparticles)[i].X[j] = molecules[i].position[j];
	    }
	    for(Molecule::atom_iter a=molecules[i].atoms.begin();
		a < molecules[i].atoms.end();
		a++)
	    {
		energy += ((*pbodies)[atomIndex].TRG[0] * (*pbodies)[atomIndex].SRC);
		atomIndex++; 
	    }
	}

	initTargets(pvbodies);
        *pvbounds = pvboundBox->getBounds(*pvbodies);
	*pvcells = pvbuildTree->buildTree(*pvbodies, *pbuffer, *pvbounds);
	VDW->evaluate(*pvcells,*pvcells);
	atomIndex = 0;
	for (unsigned int i = 0; i < molecules.size(); i++)
	{
	    for (unsigned int j = 0; j < box.dimension; j++)
	    {
		(*pparticles)[i].X[j] = molecules[i].position[j];
	    }
	    for(Molecule::atom_iter a=molecules[i].atoms.begin();
		a < molecules[i].atoms.end();
		a++)
	    {
		energy += (*pvbodies)[atomIndex].TRG[0];
		atomIndex++; 
	    }
	}
	
	energy = energy * c_ftov;

	if (energy!=energy) {
	    std::cout << particle << " produced: " << energy << std::endl;
	}
    }    
    double nonPair;
    nonPair = this->nonPairwiseCallback(particle,position,orientation);  
    energy += nonPair;
    //std::cout << " energy " << std::setprecision(16) << energy << std::endl;
    return energy;
}

double MMolecule::computePairEnergy(unsigned int particle1,
				    double position1[],
				    double orientation1[],
				    unsigned int particle2,
				    double position2[],
				    double orientation2[])
{
    // the intention is to only include particles that are very close here
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
	double energy = 0.0;
        
	for(Molecule::atom_iter MCP1 = molecules[particle1].atoms.begin();
	    MCP1 < molecules[particle1].atoms.end();
	    MCP1++) {
	    for(Molecule::atom_iter MCP2 = molecules[particle2].atoms.begin();
		MCP2 < molecules[particle2].atoms.end();
		MCP2++) {
		std::vector<double> apos1(box.dimension);
		std::vector<double> apos2(box.dimension);
		Molecule::rotate2D(MCP1->position, apos1, orientation1[0]);
		Molecule::rotate2D(MCP2->position, apos2, orientation2[0]);
		for (unsigned int i=0;i<box.dimension;i++)	    
		{
		    apos1[i] += position1[i] * molecules.particleDiameter;
		    apos2[i] += position2[i] * molecules.particleDiameter;
		    sep[i] = apos1[i] - apos2[i];
		}

		
		double normSqd = 0;

		// Calculate squared norm of vector.
		for (unsigned int i=0;i<box.dimension;i++)
		    normSqd += sep[i]*sep[i];

		double r2Inv = 1.0 / normSqd;
		double r6Inv = r2Inv*r2Inv*r2Inv;
		energy += 4.0*interactionEnergy*((r6Inv*r6Inv) - r6Inv - potentialShift);
		energy += (MCP1->charge * MCP2->charge) / sqrt(normSqd); //coulomb term
	    }
	}
	return energy;
    }
    return 0;
}

double MMolecule::nonPairwiseCallback(unsigned int particle,
				      double position[],
				      double orientation[])
{
  //this is where the surface potential will go
  return 0.0;
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
  moved=1;
}

double MMolecule::getEnergy()
{
    double energy = 0;

    for (unsigned int i=0;i<molecules.size();i++)
    {
        energy += computeEnergy(i, &molecules[i].position[0], &molecules[i].orientation[0]);
	//std::cout << std::setprecision(12) << energy << std::endl;
    }
    std::cout << std::setprecision(12) << energy/(2*molecules.size()) << std::endl;
    return energy/(2*molecules.size());
}

