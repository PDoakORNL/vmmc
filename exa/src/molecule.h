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

#ifndef _MOLECULE_H
#define _MOLECULE_H

#include <vector>
#include <memory>
#include "Particle.h"

/*! \file molecule.h
    \brief A simple particle data type.
*/



//! Structure containing attributes for an individual particle.
struct Molecule : public Particle
{
public:
    //! Default constructor.
    Molecule();

    //! Constructor.
    /*! \param index
            The particle index.
     */
    Molecule(unsigned int);
    
    unsigned int species;

    class Atom {
    public:
	int element;
        double charge;
	std::vector<double> position;
    };
    typedef std::vector<Atom>::iterator atom_iter;
    std::vector<Atom> atoms;

    double selfEnergy = 0.0;

    static void rotate2D(const std::vector<double>&, std::vector<double>&, double);
    static void rotate3D(std::vector<double>&, std::vector<double>&, std::vector<double>&, double);
    void get_apos(Atom&, double, std::vector<double>&);

};


//typedef std::vector<std::shared_ptr<Molecule>> Molecules;
//typedef std::shared_ptr<Particle> PParticle;

#endif  /* _MOLECULE_H */
