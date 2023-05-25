#include "../include/Particle.h"


namespace shift_lib
{

	std::ostream& operator << (std::ostream& out, const Particle & a)
	{
		out << a.p.t() << "   "
			<< a.p.x() << "   "
			<< a.p.y() << "   "
			<< a.p.z() << "   "
			<< a.x.t() << "   "	//t
			<< a.x.x() << "   "	//x
			<< a.x.y() << "   "
			<< a.x.z() << std::endl;
		return out;
	}

}
