#ifndef SHIFT_LIB_RANDOM_EVENTS_H
#define SHIFT_LIB_RANDOM_EVENTS_H

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <random>

#include "./distribution_3D.h"
#include "./ParticleRecord.h"
#include "./ParameterReader.h"

using namespace std;

namespace shift_lib
{

	void generate_events(vector<ParticleRecord> & allParticles, ParameterReader * paraRdr)
	{
		//cout << "Using random number generator for toy model calculation!" << endl;

		allParticles.clear();

		double mass 	= paraRdr->getVal("mass");
		double RNG_R 	= paraRdr->getVal("RNG_R");
		double RNG_p0 	= paraRdr->getVal("RNG_p0");

		int RNG_Nev 	= paraRdr->getVal("RNG_Nev");
		int RNG_mult 	= paraRdr->getVal("RNG_mult");
		int RNG_xDir 	= paraRdr->getVal("RNG_xDir");
		int RNG_yDir 	= paraRdr->getVal("RNG_yDir");
		int RNG_zDir 	= paraRdr->getVal("RNG_zDir");

		int RNG_seed	= (int)paraRdr->getVal("RNG_seed");

		unsigned time_seed = chrono::system_clock::now().time_since_epoch().count();
		int seed_to_use = ( RNG_seed < 0 ) ? time_seed : RNG_seed;
		default_random_engine generator = default_random_engine (seed_to_use);

		normal_distribution<double> distribution(0.0, 1.0);

		// use the 3D distribution here
		param_list parameters ( { {"scale", RNG_p0},
															{"do_x", (bool)RNG_xDir},
															{"do_y", (bool)RNG_yDir},
															{"do_z", (bool)RNG_zDir} } );
		distribution_3D dist3D( "Gaussian", seed_to_use, parameters );

		//for (int iEvent = 0; iEvent < RNG_Nev; ++iEvent)
		//{
		//	EventRecord event;
		//
			for (int iParticle = 0; iParticle < RNG_mult; ++iParticle)
			{

				double tP = 0.0;	// cf. paper
				double xP = RNG_xDir ? RNG_R * distribution(generator) : 0.0;
				double yP = RNG_yDir ? RNG_R * distribution(generator) : 0.0;
				double zP = RNG_zDir ? RNG_R * distribution(generator) : 0.0;

				// double px = RNG_xDir ? RNG_p0 * distribution(generator) : 0.0;
				// double py = RNG_yDir ? RNG_p0 * distribution(generator) : 0.0;
				// double pz = RNG_zDir ? RNG_p0 * distribution(generator) : 0.0;
				// std::tie(px,py,pz)
				auto [px,py,pz] = dist3D.generate_sample();
cout << "px,py,pz = " << px << "   " << py << "   " << pz << endl;
if (1) std::terminate();
				double Ep = sqrt( mass*mass + px*px + py*py + pz*pz );

				ParticleRecord particle;
				//particle.eventID 	= iEvent;
				particle.particleID = iParticle;
				particle.p.e(Ep);
				particle.p.px(px);
				particle.p.py(py);
				particle.p.pz(pz);
				particle.x.e(tP);
				particle.x.px(xP);
				particle.x.py(yP);
				particle.x.pz(zP);

				particle.m = mass;
				particle.m2 = mass*mass;

				//event.particles.push_back( particle );
				allParticles.push_back( particle );

			}

		//	allEvents.push_back( event );
		//
		//}

		return;
	}
}

#endif
