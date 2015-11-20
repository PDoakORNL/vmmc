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

#include "Exa.h"
#include "VMMC.h"
#include "MMolecule.hpp"
#include <omp.h>

int main(int argc, char** argv)
{
    // Simulation parameters.
    unsigned int dimension = 2;                     // dimension of simulation box
    unsigned int nMolecules = 300;                 // number of atoms
    double interactionEnergy = 2;                   // interaction energy scale (in units of kBT)
    double interactionRange = 4.5;                  // size of interaction range (in units of particle diameter)
    double density = 0.05;                          // particle density
    double baseLength;                              // base length of simulation box
    unsigned int maxInteractions = 100;             // maximum number of interactions per particle

    omp_set_num_threads(4);
    // std::cout << "threads: " << omp_get_num_threads() << " " << args.threads
    // 	      << std::endl;
    
    // Data structures.
    std::vector<Molecule> vmolecules;    // Molecule is derived from Particle
    CellList cells;                                 // cell list
    // Resize particle container.
    vmolecules.resize(nMolecules);
    Molecules molecules(&vmolecules);
    // Work out base length of simulation box (particle diameter is one).
    if (dimension == 2) baseLength = std::pow((nMolecules*M_PI)/(2.0*density), 1.0/2.0);
    else baseLength = std::pow((nMolecules*M_PI)/(6.0*density), 1.0/3.0);

    std::vector<double> boxSize;
    for (unsigned int i=0;i<dimension;i++)
        boxSize.push_back(baseLength);

    // Initialise simulation box object.
    Box box(boxSize);

    // Initialise input/output class.
    InputOutput io;

    // Create VMD bounding box.
    io.vmdScript(boxSize);

    // Initialise cell list.
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    // Initialise the MMolecule potential model.
    MMolecule mmolecule(box, molecules, cells,
        maxInteractions, interactionEnergy, interactionRange);
    // Initialise random number generator.
    MersenneTwister rng;

    // Initialise particle initialisation object.
    InitialiseMolecules initialise;

    // Generate a random particle configuration.
    initialise.random(molecules, cells, box, rng, false);

    // setup the exafmm bodies in molecule from vmmc particles
    mmolecule.initBodies();
   
    // Initialise data structures needed by the VMMC class.
    double coordinates[dimension*nMolecules];
    double orientations[dimension*nMolecules];
    bool isIsotropic[nMolecules];
    
    // Copy particle coordinates and orientations into C-style arrays.
    for (unsigned int i=0;i<nMolecules;i++)
    {
        for (unsigned int j=0;j<dimension;j++)
        {
            coordinates[dimension*i + j] = molecules[i].position[j];
            orientations[dimension*i + j] = molecules[i].orientation[j];
        }

        isIsotropic[i] = false;

    }

    // Initialise the VMMC callback functions.
    using namespace std::placeholders;
    vmmc::CallbackFunctions callbacks;

    
    callbacks.energyCallback =
        std::bind(&MMolecule::computeEnergy, &mmolecule, _1, _2, _3);
    callbacks.pairEnergyCallback =
        std::bind(&MMolecule::computePairEnergy, &mmolecule, _1, _2, _3, _4, _5, _6);
    callbacks.nonPairwiseCallback =
      std::bind(&MMolecule::nonPairwiseCallback, &mmolecule, _1, _2, _3);
    callbacks.interactionsCallback =
        std::bind(&MMolecule::computeInteractions, &mmolecule, _1, _2, _3, _4);
    callbacks.postMoveCallback =
        std::bind(&MMolecule::applyPostMoveUpdates, &mmolecule, _1, _2, _3);

    // Initialise the VMMC object.
#ifndef ISOTROPIC
    vmmc::VMMC vmmc(nMolecules, dimension, coordinates, orientations,
        0.15, 0.2, 0.5, 0.5, maxInteractions, &boxSize[0], isIsotropic, true, callbacks);
#else
    vmmc::VMMC vmmc(nMolecules, dimension, coordinates,
        0.15, 0.2, 0.5, 0.5, maxInteractions, &boxSize[0], true, callbacks);
#endif

    // Execute the simulation.
    for (unsigned int i=0;i<1000;i++)
    {
        // Increment simulation by 1000 Monte Carlo Sweeps.
        //std::vector<Particle> tempAtoms(atoms);
        vmmc+= 100 * nMolecules;
	// std::vector<Particle>::iterator i_atoms, i_tempAtoms;
	// i_atoms = atoms.begin();
	// i_tempAtoms = tempAtoms.begin();
	// while(i_tempAtoms != tempAtoms.end()) {
	//   if (i_atoms->position != i_tempAtoms->position) {
	//     std::cout << i_atoms->index << std::endl;
	//     std::cout << std::setprecision(12) << std::setw(14) << i_atoms->position[0] << i_atoms->position[1] << i_atoms->position[2] << std::endl;
	//     std::cout << i_tempAtoms->position[0] << std::endl;
	    
	//   }
	//   i_tempAtoms++;
	//   i_atoms++;
	// }
	
        // Append particle coordinates to an xyz trajectory.
        if (i == 0) io.appendXyzTrajectory(dimension, molecules, true);
        else io.appendXyzTrajectory(dimension, molecules, false);
 
        // Report.
        printf("sweeps = %9.4e, energy = %10f\n", ((double) (i+1)*100), mmolecule.getEnergy());
	printf("accepts = %6llu\n", vmmc.getAccepts());
    }

    std::cout << "\nComplete!\n";

    // We're done!
    return (EXIT_SUCCESS);
}
