#ifndef RANDOM_EVENTS_H
#define RANDOM_EVENTS_H

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <random>

#include "./ParticleRecord.h"
#include "./ParameterReader.h"

using namespace std;

void generate_events(vector<ParticleRecord> & allParticles, ParameterReader * paraRdr)
{
	//cout << "Using random number generator for toy model calculation!" << endl;

	allParticles.clear();

	double mass 	= paraRdr->getVal("mass");
	double RNG_R 	= paraRdr->getVal("RNG_R");
	double RNG_a 	= paraRdr->getVal("RNG_a");

	int RNG_Nev 	= paraRdr->getVal("RNG_Nev");
	int RNG_mult 	= paraRdr->getVal("RNG_mult");
	int RNG_xDir 	= paraRdr->getVal("RNG_xDir");
	int RNG_yDir 	= paraRdr->getVal("RNG_yDir");
	int RNG_zDir 	= paraRdr->getVal("RNG_zDir");

	bool RNG_seed	= (bool)paraRdr->getVal("RNG_seed");

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator;
	if ( RNG_seed )
		generator = default_random_engine (seed);

	normal_distribution<double> distribution(0.0, RNG_R);

	//for (int iEvent = 0; iEvent < RNG_Nev; ++iEvent)
	//{
	//	EventRecord event;
	//
		for (int iParticle = 0; iParticle < RNG_mult; ++iParticle)
		{
			
			double tP = 0.0;	// cf. paper
			double xP = RNG_xDir ? distribution(generator) : 0.0;
			double yP = RNG_yDir ? distribution(generator) : 0.0;
			double zP = RNG_zDir ? distribution(generator) : 0.0;

			double px = 0.2*RNG_xDir ? distribution(generator) : 0.0;
			double py = 0.2*RNG_yDir ? distribution(generator) : 0.0;
			double pz = 0.2*RNG_zDir ? distribution(generator) : 0.0;
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

			//event.particles.push_back( particle );
			allParticles.push_back( particle );

		} 

	//	allEvents.push_back( event );
	//
	//}

	return;
}

#endif
