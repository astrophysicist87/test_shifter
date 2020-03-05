#include "../include/ParticleRecord.h"

std::ostream& operator << (std::ostream& out, const ParticleRecord & a)
{
    out << a.p.e() << "   "
		<< a.p.px() << "   "
		<< a.p.py() << "   "
		<< a.p.pz() << "   "
		<< a.x.e() << "   "	//t
		<< a.x.px() << "   "	//x
		<< a.x.py() << "   "
		<< a.x.pz() << std::endl;
    return out;
}
