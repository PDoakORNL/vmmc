
//! Container class for storing a list particles. since we want to be able to have objects derived from particle and not be dead in the water
//! This class contains the main cell list that is manipulated by the simulation.
#ifndef _PARTICLES_H
#define _PARTICLES_H


#include "Particle.h"

class Particles
{
public:
    //! Default constructor.
    Particles();
    Particles(std::vector<Particle> *);
    virtual ~Particles();
    
    virtual Particle& operator[](unsigned int);
    virtual long unsigned int size();
    
    std::vector<Particle> * vpParticles;
};

#endif  /* _PARTICLES_H */
