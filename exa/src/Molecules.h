#ifndef _MOLECULES_H
#define _MOLECULES_H

//! Container class for storing a list particles. since we want to be able to have objects derived from particle and not be dead in the water

#include "Particles.h"
#include "molecule.h"

class Molecules : public Particles
{
public:
    //! Default constructor.
    Molecules();
    virtual ~Molecules();
    Molecules(std::vector<Molecule> * vpMolecules);
    virtual Molecule& operator[](unsigned int);
    virtual long unsigned int size();
    
    std::vector<Molecule> * vpMolecules;
};

#endif /* _MOLECULES_H */
