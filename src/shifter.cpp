#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "../include/shifter.h"

using namespace std;

namespace shift_lib
{
	constexpr bool BE_VERBOSE = false;

	constexpr double MM2FM = 1e12;
	constexpr double FM2MM = 1e-12;
	constexpr double HBARC = 0.19733;

	constexpr double COMPRELERR = 1e-10;
	constexpr double COMPFACMAX = 1000.;
	constexpr int    NCOMPSTEP  = 10;

	constexpr double TINY       = 1e-6;


	void shifter::process_event( vector<Particle> & allParticles_in )
	{
		allParticles = allParticles_in;

		// Perform shifts.
		shiftEvent();

		// Return shifted results.
		allParticles_in	= allParticles;

		number_of_shifted_events++;

		return;
	}


	void shifter::print(int eventID, vector<Particle> & allParticles_in, const string & filename)
	{
		ofstream out(filename.c_str(), ios::app);
		out << eventID << "   " << allParticles_in.size() << "\n";
		for (const auto & particle: allParticles_in)
			out << particle.p.x() << "   " << particle.p.y() << "   " << particle.p.z() << "\n";
		out.close();
		return;
	}

	//--------------------------------
	vector<vector<double>> shifter::get_pairs( const vector<Particle> & particles )
	{
		vector<vector<double>> result;

		const int number_of_particles = particles.size();

		// get all values in vector first
		for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
		for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
		{
			Vec4 q = particles[i1].p - particles[i2].p;
			result.push_back( vector<double>({q.x(), q.y(), q.z()}) );
		}

		return result;
	}

	//--------------------------------------------------------------------------
	// Perform Bose-Einstein corrections on an event.
	void shifter::shiftEvent()
	{
		//------------------------------------------------
		// set random seed and distribution/generator needed below
		unsigned seed = ((int)paraRdr->getVal("RNG_seed") < 0) ?
										chrono::system_clock::now().time_since_epoch().count() :
										(int)paraRdr->getVal("RNG_seed");
		default_random_engine generator(seed);
		std::uniform_real_distribution<double> uniform(0.0, 1.0);


	std::normal_distribution<double>       normal(0.0, 1.0);
	for (const auto & p: allParticles) {normal(generator); normal(generator);}
	int RNG_xDir 	= paraRdr->getVal("RNG_xDir");
	int RNG_yDir 	= paraRdr->getVal("RNG_yDir");
	int RNG_zDir 	= paraRdr->getVal("RNG_zDir");
	const double RNG_p0 = paraRdr->getVal("RNG_p0");


		//------------------------------------------------
		// set up spectra sampler
		QuantumSampler qs( allParticles, paraRdr->getVal("sigma"),
												(int)paraRdr->getVal("RNG_xDir"),
												(int)paraRdr->getVal("RNG_yDir"),
												(int)paraRdr->getVal("RNG_zDir"), seed );

		//------------------------------------------------
		// save original version of particles
		const int number_of_particles = allParticles.size();
		vector<Particle> allParticles_Original = allParticles;

		//------------------------------------------------
		// compute current pairs (shifted pairs needed below)
		current_pairs = get_pairs(allParticles);
		shifted_pairs = current_pairs;

		//------------------------------------------------
		// set tabulated BE distances for below
		double R = paraRdr->getVal("RNG_R") / HBARC;
		BoseEinsteinDistance BEdist( "SingleScale", { {"R", R} } );
		vector<double> current_BE_distances(current_pairs.size(), 0.0);
		std::transform( current_pairs.cbegin(), current_pairs.cend(),
	                  current_BE_distances.begin(), BEdist.get_distance_v );
		vector<double> shifted_BE_distances = current_BE_distances;

		// if (true)
		// {
		// cout << "R: "<< R << endl;
		// cout << "BEdistances:" << endl;
		// int total_n = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*current_BE_distances.size())));
		// int index = 0;
		// for (int i = 0; i < total_n; i++)
		// for (int j = i+1; j < total_n; j++)
		// {
		// 	const auto & q = current_pairs.at(index);
		// 	cout << index << "  " << i << "  " << j << "  " << q[0] << "  " << q[1] << "  " << q[2]
		// 	     << "  " << current_BE_distances.at(index) << "\n";
		// 	index++;
		// }
		// std::terminate();
		// }

		//------------------------------------------------
		// compute probability of initial configuration
		const double precision = 1e-2;
		const bool assume_sparse_pair_matrix = false;
		param_list parameters ( { {"precision", precision},
															{"n_particles", (int)allParticles.size()},
															{"assume_sparse", assume_sparse_pair_matrix} } );

// cout << "Initializing cp..." << endl;
// cout << "current_pairs.size() = " << current_pairs.size() << endl;
// cout << "current_pairs[0] = " << "  " << current_pairs[0][0] << "  " << current_pairs[0][1]
// 			<< "  " << current_pairs[0][2] << endl;
// cout << "Direct = " << BEdist.get_distance(current_pairs[0]) << endl;
// cout << "current_BE_distances[0] = " << current_BE_distances[0] << endl;
		ConfigurationProbability cp( SHIFT_MODE, parameters, paraRdr, "Short" );
		double P1 = cp.get_probability( allParticles, current_BE_distances, -1 );
																		// last -1 means place all particles in clusters
		cout << "current_BE_distances.size() = " << current_BE_distances.size() << endl;
// cout << "Finished initializing cp." << endl;
// cout << "CHECK: " << P1 << "  " << P1FP << endl;

		//------------------------------------------------
		// begin looping over particles
		Stopwatch swTotal;
		int nLoops = paraRdr->getVal("shifter_nLoops");
		for (int iLoop = 0; iLoop < nLoops; iLoop++)	// loop over all particles fixed number of times
		{
			swTotal.Reset();
			swTotal.Start();
			cerr << "Loop #" << iLoop << ": \n";
			for (int iParticle = 0; iParticle < number_of_particles; iParticle++) // loop over particles, re-sample one at a time
			{
				cerr << "  --> Particle #" << iParticle << "\n";

				// generate a shifted momentum
				double x1 = allParticles[iParticle].p.x();
				double y1 = allParticles[iParticle].p.y();
				double z1 = allParticles[iParticle].p.z();

				// QuantumSampler fixes distribution in terms of original particles
				// auto [x2, y2, z2] = qs.sample(iParticle/*, std::make_tuple(x1, y1, z1)*/);

				// use simple Gaussian for checks
				double x2 = RNG_xDir ? RNG_p0 * normal(generator) : x1;
				double y2 = RNG_yDir ? RNG_p0 * normal(generator) : y1;
				double z2 = RNG_zDir ? RNG_p0 * normal(generator) : z1;

				// set configuration containing shifted particle
				vector<Particle> allParticles_with_shift = allParticles;
				allParticles_with_shift[iParticle].p.x(x2);
				allParticles_with_shift[iParticle].p.y(y2);
				allParticles_with_shift[iParticle].p.z(z2);

				// get shifted pairs
				get_shifted_pairs( /*shifted_pairs,*/ shifted_BE_distances,
													 allParticles_with_shift, iParticle, BEdist );

				// get probability of shifted configuration
				double P2 = cp.get_probability( allParticles_with_shift,
																				shifted_BE_distances, iParticle );


			// std::cout << "CHECK: " << SHIFT_MODE << "   "
			// 					<< x1 << "   " << y1 << "   " << z1 << "   "
			// 					<< x2 << "   " << y2 << "   " << z2 << "   "
			// 					<< P2/P1 << std::endl << setprecision(4);
// cout << "--------------------------------------------------------------------------------------------------------------------" << endl;
// 			for (const auto & p: allParticles) std::cout << setw(10) << right << p.p.z() << "   ";
// 			std::cout << std::endl;
// 			for (const auto & p: allParticles_with_shift) std::cout << setw(10) << right << p.p.z() << "   ";
// 			std::cout << std::endl;
// double tmp = uniform(generator);
// cout << "Acceptance probability: " << P2/P1 << " vs. " << tmp << endl;


				// choose new configuration (shifted or original)
				double log_alpha = std::min(0.0, log(P2/P1));
				bool shift_this_particle = (log(uniform(generator)) <= log_alpha);
				// bool shift_this_particle = (log(tmp) <= log_alpha);

				// re-set shifted quantities
				if (shift_this_particle)
				{
					// current_pairs = shifted_pairs;		// re-set current configuration
					get_shifted_pairs( /*current_pairs,*/ current_BE_distances,
														 allParticles_with_shift, iParticle, BEdist );
					P1 = P2;													// re-set probability of current configuration
					allParticles[iParticle].p.x(x2);	// re-set momentum of particle
					allParticles[iParticle].p.y(y2);	// re-set momentum of particle
					allParticles[iParticle].p.z(z2);	// re-set momentum of particle
				}
				else	// overwrite shifted pairs with current (unshifted)
				{
					// current particles vector used to overwrite shifted pairs and distances
					get_shifted_pairs( /*shifted_pairs,*/ shifted_BE_distances,
														 allParticles, iParticle, BEdist );

					// unshifted particles and distances overwrite current internal state
					// of ConfigurationProbability cp object, if any
					cp.revert_state( allParticles, current_BE_distances );
				}

			}
			swTotal.Stop();
			cout << "finished in " << swTotal.printTime() << " s.\n";
		}

		// Done.
		return;
	}


	//--------------------------------
	void shifter::get_shifted_pairs( //vector<vector<double>> & pairs,
																	 vector<double> & BE_distances,
																   const vector<Particle> & particles,
															     const int shifted_particle_index,
 																	 const BoseEinsteinDistance & BEdist )
	{
		const int np = particles.size();
		auto & spi = shifted_particle_index;

		// indexes upper-triangular list of pairs given indices (i,j)
		auto UTindexer = [](int i, int j, int n){return -1 + j - i*(3 + i - 2*n)/2;};

		for (int i1 = 0; i1 < spi; ++i1)
		{
			// int iPair = UTindexer(i1, spi, np);
			// auto q = { particles[i1].p.x() - particles[spi].p.x(),
			// 					 particles[i1].p.y() - particles[spi].p.y(),
			// 					 particles[i1].p.z() - particles[spi].p.z() };
			// BE_distances[iPair] = BEdist.get_distance(vector<double>(q));
			BE_distances[UTindexer(i1, spi, np)]
				= BEdist.get_distance({ particles[i1].p.x() - particles[spi].p.x(),
								 particles[i1].p.y() - particles[spi].p.y(),
								 particles[i1].p.z() - particles[spi].p.z() });
			// pairs[iPair].assign(q);
		}

		for (int i1 = spi+1; i1 < np; ++i1)
		{
			// int iPair = UTindexer(spi, i1, np);
			// auto q = { particles[spi].p.x() - particles[i1].p.x(),
			// 					 particles[spi].p.y() - particles[i1].p.y(),
			// 					 particles[spi].p.z() - particles[i1].p.z() };
			BE_distances[UTindexer(spi, i1, np)]
				= BEdist.get_distance({ particles[spi].p.x() - particles[i1].p.x(),
								 particles[spi].p.y() - particles[i1].p.y(),
								 particles[spi].p.z() - particles[i1].p.z() });
			// pairs[iPair].assign(q);
		}

		return;
	}

}	//End of namespace

//End of file
