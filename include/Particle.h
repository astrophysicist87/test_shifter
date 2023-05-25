#ifndef SHIFT_LIB_PARTICLE_H
#define SHIFT_LIB_PARTICLE_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "FourVector.h"

using namespace std;

namespace shift_lib
{

	typedef struct
	{
		int particleID;

		double m, m2;

		Vec4 x, p;
	} Particle;

	std::ostream& operator << (std::ostream& out, const Particle & a);

	typedef std::pair<Particle, Particle> particle_pair;

}

#endif
