#include "Molecules.h"

Molecules::Molecules() {
}

Molecules::Molecules(std::vector<Molecule> * molecules_) : vpMolecules(molecules_)
{
}

Molecules::~Molecules()
{
}


Molecule& Molecules::operator[](unsigned int index) {
	return (*vpMolecules)[index];
}

long unsigned int Molecules::size() {
	return vpMolecules->size();
}
