#include "Particles.h"

Particles::Particles() {
}

Particles::Particles(std::vector<Particle> * particles_) : vpParticles(particles_)
{
}

Particles::~Particles()
{
}

Particle& Particles::operator[](unsigned int index) {
	return (*vpParticles)[index];
}

long unsigned int Particles::size() {
	return vpParticles->size();
}
