#ifndef SHIFT_LIB_SHIFTER_H
#define SHIFT_LIB_SHIFTER_H

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <omp.h>

#include "ParameterReader.h"
#include "Particle.h"
#include "BoseEinsteinDistance.h"
#include "probability.h"
#include "QuantumSampler.h"
#include "stopwatch.h"

using namespace std;

namespace shift_lib
{

	template <class T>
	inline T sgn(T v) { return T(v > T(0)) - T(v < T(0)); }

	//============================================================================
	class shifter
	{
		//--------------------------------------------------------------------------
		private:
			//------------------------------------------------------------------------
			const string SHIFT_MODE;

			ParameterReader * paraRdr;

			vector<Particle> allParticles;
			vector<Particle> allParticles_Shifted;

			// miscellaneous
			string path;
			ostream & out;
			ostream & err;

			int number_of_shifted_events = 0;

			// need these objects for finding combinations
			vector<int> indices;

			vector<vector<double>> current_pairs, shifted_pairs;

			//------------------------------------------------------------------------
			void process_event( vector<Particle> & allParticles_in );

			void shiftEvent();

			bool setSortedPairs( const vector<Particle> & particles_to_sort );

			vector<vector<double>> get_pairs(
				const vector<Particle> & particles );

			void get_shifted_pairs(
				vector<vector<double>> & pairs,
				vector<double> & BE_distances,
			  const vector<Particle> & particles,
		    const int shifted_particle_index,
				const BoseEinsteinDistance & BEdist );

		//--------------------------------------------------------------------------
		public:
			//------------------------------------------------------------------------
			// Constructors, destructors, and initializers
			shifter( ParameterReader * paraRdr_in,
				vector<Particle> & allParticles_in,
				std::string & shift_mode,
				ostream & out_stream = std::cout,
				ostream & err_stream = std::cerr )
				:
				paraRdr{ paraRdr_in },
				SHIFT_MODE{ shift_mode },
				out{ out_stream },
				err{ err_stream }
				{
					process_event( allParticles_in );
				};

			~shifter(){}

			void print(int eventID, vector<Particle> & allParticles_in, const string & filename);
	};

}

// End of file
#endif
