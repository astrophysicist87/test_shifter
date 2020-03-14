#ifndef PARTICLERECORD_H
#define PARTICLERECORD_H

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "FourVector.h"

using namespace std;

namespace shift_lib
{

	typedef struct
	{
		int particleID;

		double m, m2;

		Vec4 x, p, pShift, pComp;	
	} ParticleRecord;

	std::ostream& operator << (std::ostream& out, const ParticleRecord & a);

}

#endif
