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

#ifndef _MMOLECULE_H
#define _MMOLECULE_H

#include "Model.h"
#include "types.h"
#include "nodeinfo.h"
#include "bound_box.h"
#include "traversal.h"
#include "up_down_pass.h"
#include "build_tree.h"

/*! \file MMolecule.h
*/

using exafmm::Bodies;
using exafmm::BoundBox;
using exafmm::Bounds;
using exafmm::BuildTree;
using exafmm::Cells;
using exafmm::Traversal;
using exafmm::UpDownPass;

//! Class defining the MONC Molecule potential.
class MMolecule : public Model
{
public:
    //! Constructor.
    /*! \param box_
            A reference to the simulation box object.

        \param particles_
            A reference to the particle list.

        \param cells_
            A reference to the cell list object.

        \param maxInteractions_
            The maximum number of interactions per particle.

        \param interactionEnergy_
            The potential energy scale (in units of kBT).

        \param interactionRange_
            The potential cut-off distance.
     */
    MMolecule(Box&, std::vector<Particle>&, CellList&, unsigned int, double, double, NodeInfo * = NULL);

  virtual double computeEnergy(unsigned int, double[], double[]);
  
    //! Calculate the pair energy between two particles.
    /*! \param particle1
            The index of the first particle.

        \param position1
            The position vector of the first particle.

        \param orientation1
            The orientation vector of the first particle.

        \param position2
            The position vector of the second particle.

        \param orientation2
            The orientation vector of the second particle.

        \return
            The pair energy between particles 1 and 2.
     */
    double computePairEnergy(unsigned int, double[], double[],
			     unsigned int, double[], double[]);


    double nonPairwiseCallback(unsigned int, double[], double[]);
    
    void applyPostMoveUpdates(unsigned int, double[], double[]);

    void updateBodies();
    void initBodies();
    bool accepted = true;
    
private:
    const exafmm::real_t cycle = 2 * M_PI;

    struct {
	int ncrit = 64;
	int dual = 1;
	int full = 0;
	int graft = 0;
	int getMatrix = 0;
	int images = 0;
	int verbose = 0;
	int mutual = 0;
	int useRopt = 0;
	int nspawn = 5000;
	int threads = 1;
	double theta =.4;
	int useRmax = 0;
    } args;

    Bodies *pbodies;
    Bodies *pbodies2;
    Bodies *pjbodies;
    Bodies *pbuffer;
    BoundBox *pboundBox;
    Bounds *pbounds;
    BuildTree *pbuildTree;
    Cells *pcells;
    Cells *pjcells;
    Traversal *ptraversal;
    UpDownPass *pupDownPass;

    double potentialShift;  //!< Shift factor to zero potential at cut-off.
    
};

#endif  /* _MMOLECULE_H */

