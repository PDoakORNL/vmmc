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

#include "Model.h"

double INF = std::numeric_limits<double>::infinity();

Model::Model(Box& box_, std::vector<Particle>& particles_, CellList& cells_,
    unsigned int maxInteractions_, double interactionEnergy_, double interactionRange_) :
        box(box_),
        particles(particles_),
        cells(cells_),
        maxInteractions(maxInteractions_),
        interactionEnergy(interactionEnergy_),
        interactionRange(interactionRange_)
{
    // Work out squared cut-off distance.
    squaredCutOffDistance = interactionRange * interactionRange;
}

double Model::computeEnergy(unsigned int particle, double position[], double orientation[])
{
    // N.B. This method is somewhat redundant since the same functionality
    // could be achieved by using a combination of the computeInteractions
    // and model specific computePairEnergy methods.

    unsigned int cell;      // cell index
    unsigned int neighbour; // index of neighbouring particle
    double energy = 0;      // energy counter

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        cell = cells[particles[particle].cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            neighbour = cells[cell].particles[j];

            // Make sure the particles are different.
            if (neighbour != particle)
            {
                // Calculate model specific pair energy.
                energy += computePairEnergy(particle, position, orientation,
                          neighbour, &particles[neighbour].position[0],
                          &particles[neighbour].orientation[0]);

                // Early exit test for hard core overlaps and large finite energy repulsions.
                if (energy > 1e6) return INF;
            }
        }
    }

    return energy;
}

double Model::computePairEnergy(unsigned int, double[], double[], unsigned int, double[], double[])
{
    std::cerr << "[ERROR] Model: Virtual function Model::computePairEnergy() must be defined.\n";
    exit(EXIT_FAILURE);
}

unsigned int Model::computeInteractions(unsigned int particle,
    double position[], double orientation[], unsigned int interactions[])
{
    unsigned int cell;              // cell index
    unsigned int neighbour;         // index of neighbouring particle
    unsigned int nInteractions = 0; // interaction counter

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        cell = cells[particles[particle].cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            neighbour = cells[cell].particles[j];

            // Make sure the particles are different.
            if (neighbour != particle)
            {
                std::vector<double> sep(box.dimension);

                // Compute separation.
                for (unsigned int k=0;k<box.dimension;k++)
                    sep[k] = position[k] - particles[neighbour].position[k];

                // Enforce minimum image.
                box.minimumImage(sep);

                double normSqd = 0;

                // Calculate squared norm of vector.
                for (unsigned int k=0;k<box.dimension;k++)
                    normSqd += sep[k]*sep[k];

                // Particles interact.
                if (normSqd < squaredCutOffDistance)
                {
                    if (nInteractions == maxInteractions)
                    {
                        std::cerr << "[ERROR] Model: Maximum number of interactions exceeded!\n";
                        exit(EXIT_FAILURE);
                    }

                    interactions[nInteractions] = neighbour;
                    nInteractions++;
                }
            }
        }
    }

    return nInteractions;
}

void Model::applyPostMoveUpdates(unsigned int particle, double position[], double orientation[])
{
    // Copy coordinates/orientations.
    for (unsigned int i=0;i<box.dimension;i++)
    {
        particles[particle].position[i] = position[i];
        particles[particle].orientation[i] = orientation[i];
    }

    // Calculate the particle's cell index.
    unsigned int newCell = cells.getCell(particles[particle]);

    // Update cell lists if necessary.
    if (particles[particle].cell != newCell)
        cells.updateCell(newCell, particles[particle], particles);
}
