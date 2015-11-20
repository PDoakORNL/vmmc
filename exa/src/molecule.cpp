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
 
#include "molecule.h"
#include <cmath>
#include <cassert>

Molecule::Molecule() 
{
}

Molecule::Molecule(unsigned int index_) : Particle(index_)
{
}

void Molecule::rotate2D(const std::vector<double>& v1, std::vector<double>& v2, double angle)
{
        double c = std::cos(angle);
        double s = std::sin(angle);

        v2[0] = (v1[0]*c - v1[1]*s);// - v1[0]; // I don't really think these subtraction terms belong since I'm always rotating about 0
        v2[1] = (v1[0]*s + v1[1]*c);// - v1[1];
}

void Molecule::rotate3D(std::vector<double>& v1, std::vector<double>& v2, std::vector<double>& v3, double angle)
    {
	    double c = std::cos(angle);
	    double s = std::sin(angle);

        double v1Dotv2 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

        v3[0] = ((v1[0] - v2[0]*v1Dotv2))*(c - 1) + (v2[2]*v1[1] - v2[1]*v1[2])*s;
        v3[1] = ((v1[1] - v2[1]*v1Dotv2))*(c - 1) + (v2[0]*v1[2] - v2[2]*v1[0])*s;
        v3[2] = ((v1[2] - v2[2]*v1Dotv2))*(c - 1) + (v2[1]*v1[0] - v2[0]*v1[1])*s;
    }

void Molecule::get_apos(Atom& atom, double particleDiameter, std::vector<double>& apos)
{
	assert(position.size()==apos.size());
	Molecule::rotate2D(atom.position, apos, orientation[0]);
	for (unsigned int i=0;i<position.size();i++)	    
	{
	    apos[i] += position[i]*particleDiameter;
	}
}
