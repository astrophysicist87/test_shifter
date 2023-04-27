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

			vector<double> get_pairs( const vector<ParticleRecord> & particles );
			vector<double> get_shifted_pairs( const vector<double> & pairs,
			 																	const vector<ParticleRecord> & particles,
																				const int shifted_particle_index );

			void get_combinations(int N, int K, vector<vector<int>> & combinations);

			double standard_deviation( const vector<double> & v );

			double get_probability( const double R, const vector<double> & pair_qzs );

			double get_RMSscale( const vector<ParticleRecord> & particles );


			// long long BinomialCoefficient(const long long n, const long long k)
			// {
			//   std::vector<long long> aSolutions(k);
			//   aSolutions[0] = n - k + 1;
			//
			//   for (long long i = 1; i < k; ++i)
			//     aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
			//
			//   return aSolutions[k - 1];
			// }

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
