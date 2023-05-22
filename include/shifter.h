#ifndef SHIFT_LIB_SHIFTER_H
#define SHIFT_LIB_SHIFTER_H

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include <omp.h>

#include "ParameterReader.h"
#include "ParticleRecord.h"

using namespace std;

namespace shift_lib
	{

	template <class T>
	inline T
	sgn(T v) {
		return T(v > T(0)) - T(v < T(0));
	}


	class shifter
	{
		private:

			const string SHIFT_MODE;

			ParameterReader * paraRdr;

			vector<ParticleRecord> allParticles;
			vector<ParticleRecord> allParticles_Shifted;

			// miscellaneous
			string path;
			ostream & out;
			ostream & err;

			bool include_pair_density;

			vector<double> pairShifts;
			vector<bool> this_pair_shifted;

			// The pair density
			vector<double> denBar;

			// Various ways to sort pairs
			vector< pair< double, pair <int,int> > > pairs_sorted_by_qz;
			vector< pair< double, pair <int,int> > > pairs_sorted_by_abs_qz;

			int number_of_shifted_events = 0;

			// need these objects for finding combinations
			std::mt19937 comb_gen;
			vector<int> indices;

			vector<vector<double>> current_pairs, shifted_pairs;

		public:

			// Constructors, destructors, and initializers
			shifter( ParameterReader * paraRdr_in,
						vector<ParticleRecord> & allParticles_in,
						std::string & shift_mode,
						ostream & out_stream = std::cout,
						ostream & err_stream = std::cerr )
						:
						SHIFT_MODE(shift_mode),
						out(out_stream),
						err(err_stream),
						comb_gen(std::mt19937{std::random_device{}()})
						{ initialize_all( paraRdr_in, allParticles_in ); };


			void initialize_all( ParameterReader * paraRdr_in, vector<ParticleRecord> & allParticles_in );

			void update_records( vector<ParticleRecord> & allParticles_in );

			void process_event( vector<ParticleRecord> & allParticles_in );

			void print(int eventID, vector<ParticleRecord> & allParticles_in, const string & filename);

			~shifter();

			// void shiftEvent();
			void shiftEvent_efficient();

			bool setSortedPairs( const vector<ParticleRecord> & particles_to_sort );

			vector<vector<double>> get_pairs(
				const vector<ParticleRecord> & particles );

			vector<vector<double>> get_shifted_pairs(
				const vector<vector<double>> & pairs,
				const vector<ParticleRecord> & particles,
				const int shifted_particle_index );

			void get_shifted_pairs_FAST(
				vector<vector<double>> & pairs,
				vector<double> & BE_distances,
			  const vector<ParticleRecord> & particles,
		    const int shifted_particle_index,
        const double R );

			double standard_deviation( const vector<double> & v );

			double get_probability( const double R, const vector<vector<double>> & qVec, const vector<double> & BE_distances );
			double get_probability( const double R, const vector<vector<double>> & qVec, const string & CHOSEN_SHIFT_MODE );


			double get_RMSscale( const vector<ParticleRecord> & particles );

			inline vector<long> dec2binarr(long n, long dim)
			{
			    // note: res[dim] will save the sum res[0]+...+res[dim-1]
			    // long* res = (long*)calloc(dim + 1, sizeof(long));
			    vector<long> res(dim+1);
			    long pos = dim - 1;

			    // note: this will crash if dim < log_2(n)...
			    while (n > 0)
			    {
			        res[pos] = n % 2;
			        res[dim] += res[pos];
			        n = n / 2; // integer division
			        pos--;
			    }

			    return res;
			}

			double permanent(const vector<double> & A, long n);

	};

}

// End of file
#endif
