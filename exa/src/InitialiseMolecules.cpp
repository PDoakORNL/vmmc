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

#include "InitialiseMolecules.h"

InitialiseMolecules::InitialiseMolecules()
{
}

void InitialiseMolecules::random(Molecules& molecules, CellList& cells, Box& box, MersenneTwister& rng, bool isSpherocylinder)
{
    if (isSpherocylinder && (box.dimension != 3))
    {
        std::cerr << "[ERROR] InitialiseMolecules: Spherocylindrical boundary only valid for three dimensional simulation box!\n";
        exit(EXIT_FAILURE);
    }

    // Copy box dimensions.
    boxSize = box.boxSize;

    for (unsigned i=0;i<molecules.size();i++)
    {
        // Current number of attempted molecule insertions.
        unsigned int nTrials = 0;

        // Whether molecule overlaps.
        bool isOverlap = true;

        // Temporary vector.
        std::vector<double> vec(box.dimension);

        // Set molecule index.
        molecules[i].index = i;

	molecules[i].charge = 0.0;

	molecules[i].atoms.resize(2);
	molecules[i].atoms[0].element = 1;
	molecules[i].atoms[0].charge = 0.69;
	molecules[i].atoms[0].position = {.46,0,0};
	molecules[i].atoms[1].element = 9;
	molecules[i].atoms[1].charge = -0.69;
	molecules[i].atoms[1].position = {-.46,0,0};
		

        // Keep trying to insert molecule until there is no overlap.
        while (isOverlap)
        {
            nTrials++;

            // Generate a random position.
            for (unsigned int j=0;j<box.dimension;j++)
                vec[j] = rng()*box.boxSize[j];

            molecules[i].position = vec;

            // Generate a random orientation.
            for (unsigned int j=0;j<box.dimension;j++)
                vec[j] = rng.normal();

            // Calculate vector norm.
            double norm = 0;
            for (unsigned int j=0;j<box.dimension;j++)
                norm += vec[j]*vec[j];
            norm = sqrt(norm);

            // Convert orientation to a unit vector.
            for (unsigned int j=0;j<box.dimension;j++)
                vec[j] /= norm;

            molecules[i].orientation = vec;

            // Calculate the molecule's cell index.
            molecules[i].cell = cells.getCell((Particle&)(molecules[i]));

            // Enforce spherocylindrical boundary.
//             if (isSpherocylinder)
//             {
//                 // Make sure molecule lies within the spherocylinder.
// #ifndef ISOTROPIC
//                 if (!outsideSpherocylinder(i, &molecules[i]->position[0], &molecules[i]->orientation[0]))
// #else
//                 if (!outsideSpherocylinder(i, &molecules[i]->position[0]))
// #endif
//                 {
//                     // See if there is any overlap between molecules.
// 			isOverlap = checkOverlap(*(particles *)molecules[i], molecules, cells, box);
//                 }
//                 else isOverlap = true;
//             }
//             else
//             {
                // See if there is any overlap between molecules.

	    isOverlap = checkOverlap(molecules[i], molecules, cells, box);
            

            // Check trial limit isn't exceeded.
            if (nTrials == MAX_TRIALS)
            {
                std::cerr << "[ERROR] InitialiseMolecules: Maximum number of trial insertions reached.\n";
                exit(EXIT_FAILURE);
            }
        }

        // Update cell list.
        cells.initCell(molecules[i].cell, (Particle&)molecules[i]);
    }
}


bool InitialiseMolecules::checkOverlap(Molecule& molecule, Molecules& molecules, CellList& cells, Box& box)
{
    unsigned int cell, neighbour;

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        cell = cells[molecule.cell].neighbours[i];

        // Check all molecules within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            neighbour = cells[cell].particles[j];

            // Make sure molecules are different.
            if (neighbour != molecule.index)
            {
                // Molecule separtion vector.
                std::vector<double> sep(box.dimension);

                // Compute separation.
                for (unsigned int k=0;k<box.dimension;k++)
                    sep[k] = molecule.position[k] - molecules[neighbour].position[k];

                // Compute minimum image.
                box.minimumImage(sep);

                double normSqd = 0;

                // Calculate squared norm of vector.
                for (unsigned int k=0;k<box.dimension;k++)
                    normSqd += sep[k]*sep[k];

		
                // Overlap if normSqd is less than molecule diameter (box is scaled in diameter units).
                if (normSqd < 1) return true;
            }
        }
    }

    // If we get this far, no overlaps.
    return false;
}
