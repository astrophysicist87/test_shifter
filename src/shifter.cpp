#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <unordered_set>

#include "../include/shifter.h"
#include "../include/ParameterReader.h"

using namespace std;

namespace shift_lib
{
	// constexpr const char * SHIFT_MODE = "TRIAL3";

	constexpr bool BE_VERBOSE = false;

	constexpr double MM2FM = 1e12;
	constexpr double FM2MM = 1e-12;
	constexpr double HBARC = 0.19733;

	constexpr double COMPRELERR = 1e-10;
	constexpr double COMPFACMAX = 1000.;
	constexpr int    NCOMPSTEP  = 10;

	constexpr double TINY       = 1e-6;


	void shifter::initialize_all( ParameterReader * paraRdr_in,
		vector<ParticleRecord> & allParticles_in )
	{
		include_pair_density = false;

		// Load parameters.
		paraRdr			= paraRdr_in;

		// cout << "SHIFT_MODE = " << SHIFT_MODE << endl;

		// process (i.e., shift and output) this event
		process_event( allParticles_in );

		return;
	}

	shifter::~shifter()
	{
		//clear everything

		return;
	}

	void shifter::update_records( vector<ParticleRecord> & allParticles_in )
	{
		allParticles.clear();
		pairs_sorted_by_qz.clear();
		pairs_sorted_by_abs_qz.clear();

		// process (i.e., shift) this event
		process_event( allParticles_in );

		return;
	}

	void shifter::process_event( vector<ParticleRecord> & allParticles_in )
	{
		allParticles = allParticles_in;

		// Perform shifts.
		shiftEvent_efficient();

		// Return shifted results.
		allParticles_in	= allParticles;

		number_of_shifted_events++;

		return;
	}


	void shifter::print(int eventID, vector<ParticleRecord> & allParticles_in, const string & filename)
	{
		ofstream out(filename.c_str(), ios::app);
		out << eventID << "   " << allParticles_in.size() << "\n";
		for (const auto & particle: allParticles_in)
			out << particle.p.px() << "   " << particle.p.py() << "   " << particle.p.pz() << "\n";
		out.close();
		return;
	}


	double shifter::standard_deviation( const vector<double> & v )
	{
		const double n = v.size();
		double x2 = 0.0, runsum = 0.0, variance = 0.0;
		for (const auto e: v)
		{
			variance -= e*runsum;
			runsum += e;
			x2 += e*e;
		}
		variance *= 2.0/(n*(n-1.0));
		variance += x2/n;
		return sqrt(variance);
	}

double shifter::permanent(const vector<double> & A, long n) // expects n by n matrix encoded as vector
{
    double sum = 0.0;
    double rowsumprod = 0.0, rowsum = 0.0;
    vector<long> chi(n + 1);
    double C = (double)pow((double)2, n);

    // loop all 2^n submatrices of A
    for (long k = 1; k < C; k++)
    {
        rowsumprod = 1.0;
        chi = dec2binarr(k, n); // characteristic vector

        // loop columns of submatrix #k
        for (long m = 0; m < n; m++)
        {
            rowsum = 0.0;

            // loop rows and compute rowsum
            for (long p = 0; p < n; p++)
                rowsum += chi[p] * A[m * n + p];

            // update product of rowsums
            rowsumprod *= rowsum;

            // (optional -- use for sparse matrices)
            if (rowsumprod < TINY) break;
        }

        sum += (double)pow((double)-1, n - chi[n]) * rowsumprod;
    }

    return sum;
}


//------------------------------------------------------------------------------
double shifter::get_probability( const double R, const vector<vector<double>> & qVec )
{
	if ( SHIFT_MODE == "Exact" )
	{
		const int n = qVec.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		vector<double> A(np*np);
		{
			int index = 0;
			for (int i = 0; i < np; i++)
			{
				A[i*np+i] = 1.0;
				for (int j = i+1; j < np; j++)
				{
					auto q = qVec[index];
					double q2 = inner_product(q.cbegin(), q.cend(),
																		q.cbegin(),
																		0.0);
					double tmp = exp(-0.25*q2*R*R);
					if (tmp < TINY) tmp = 0.0;	// make matrix as sparse as possible
					A[i*np+j] = tmp;
					A[j*np+i] = tmp; // matrix is symmetric
					index++;
				}
			}
		}
		return permanent(A, np);
	}
	//--------------------------------------------------------------------------
	// FULL PRODUCT
	else if ( SHIFT_MODE == "FullProduct" )
	{
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		for (const auto & q: qVec)
		{
			double q2 = inner_product(q.cbegin(), q.cend(),
																q.cbegin(),
																0.0);
			result *= 1.0 + normalization*exp(-0.5*q2*R*R);
		}
		return result;
	}
	//--------------------------------------------------------------------------
	// EXPERIMENTAL
	else if ( SHIFT_MODE == "TRIAL3" )
	{
		// use only np-1 independent pairs, and cycle over which gets omitted
		const int n = qVec.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		int i = -1;
		double factor = 0.0;
		for (int i1 = 0; i1 < np - 1; ++i1)
		for (int i2 = i1 + 1; i2 < np; ++i2)
		{
			i++;
			bool include_this_pair = (i2 == i1+1) || (i1 == 0 && i2 == np-1);
			if (!include_this_pair) continue;
			auto q = qVec[i];
			double q2 = inner_product(q.cbegin(), q.cend(),
																q.cbegin(),
																0.0);
			double term = 1.0 + 0.5*np*normalization*exp(-0.5*q2*R*R);
			result *= term;
			factor += 1.0/term;
		}
		return factor*result/np;
	}
	//--------------------------------------------------------------------------
	else if ( SHIFT_MODE == "TRIAL5" )
	{
		// use adjacent particle pairs and next-to-neighbor pairs
		const int n = qVec.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double total = 0.0;
		double normalization = paraRdr->getVal("shifter_norm");
		int maxsep = np/2;
		for (int step = 1; step <= maxsep; step++) // sum over independent pairs (modulo step)
		{
			int i = -1;
			double result = 1.0;
			for (int i1 = 0; i1 < np - 1; ++i1)
			for (int i2 = i1 + 1; i2 < np; ++i2)
			{
				i++;
				int di = std::abs(i2-i1);
				bool include_this_pair = (std::min(di, np-di) == step);
				if (!include_this_pair) continue;
				auto q = qVec[i];
				double q2 = inner_product(q.cbegin(), q.cend(),
																	q.cbegin(),
																	0.0);
				result *= 1.0 + 0.5*(np-1.0)*normalization*exp(-0.5*q2*R*R);
			}
			total += result;
		}
		return total/maxsep;
	}
	//--------------------------------------------------------------------------
	else if ( SHIFT_MODE == "RMSscale" )
	{
		const int n = qVec.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		int i = -1;
		double factor = 0.0;
		for (int i1 = 0; i1 < np - 1; ++i1)
		for (int i2 = i1 + 1; i2 < np; ++i2)
		{
			i++;
			bool include_this_pair = (i2 == i1+1) || (i1 == 0 && i2 == np-1);
			if (!include_this_pair) continue;
			auto q = qVec[i];
			double q2 = inner_product(q.cbegin(), q.cend(),
																q.cbegin(),
																0.0);
			double term = 1.0 + 0.5*np*normalization*exp(-0.5*q2*R*R);
			result *= term;
			factor += 1.0/term;
		}
		return factor*result/np;
	}
	//--------------------------------------------------------------------------
	else
	{
		cerr << "This mode (" << SHIFT_MODE << ") not supported!" << endl;
		abort();
	}
}

//--------------------------------
vector<vector<double>> shifter::get_pairs( const vector<ParticleRecord> & particles )
{
	vector<vector<double>> result;

	const int number_of_particles = particles.size();

	// get all values in vector first
	for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
	for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
	{
		// Vec4 q = particles.at(i1).p - particles.at(i2).p;
		Vec4 q = particles[i1].p - particles[i2].p;
		result.push_back( vector<double>({q.px(), q.py(), q.pz()}) );
	}

	return result;
}


//--------------------------------------------------------------------------
// Perform Bose-Einstein corrections on an event.

void shifter::shiftEvent_efficient()
{
	// set up distributions
	int RNG_seed	= (int)paraRdr->getVal("RNG_seed");

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator;
	if ( RNG_seed < 0 )
		generator = default_random_engine (seed);
	else
		generator = default_random_engine (RNG_seed);

	std::uniform_real_distribution<double> uniform(0.0, 1.0);
	std::normal_distribution<double>       normal(0.0, 1.0);

	// this RNG is for shuffling particle momenta
	auto rng = std::default_random_engine { std::random_device{}() };

	// avoid replicating previously generated values (cycle through twice to be safe)
	for (const auto & p: allParticles) {normal(generator); normal(generator);}

	const double R = paraRdr->getVal("RNG_R") / HBARC;
	const double RMSscale = get_RMSscale(allParticles) / HBARC;

	const double RNG_p0 = paraRdr->getVal("RNG_p0");
	const int number_of_particles = allParticles.size();
	vector<ParticleRecord> allParticles_Original = allParticles;

	vector<vector<double>> current_pairs = get_pairs(allParticles);

	// need comparator for sorting particles by momentum
	auto particleSort = []( const ParticleRecord & p1,
                          const ParticleRecord & p2 )
                        { return (p1.p.px() < p2.p.px())
															&& (p1.p.py() < p2.p.py())
															&& (p1.p.pz() < p2.p.pz()); };

	// sort by pz
	std::sort(allParticles.begin(), allParticles.end(), particleSort);

	// get probability of current configuration
	double P1 = (SHIFT_MODE == "RMSscale") ?
							get_probability( RMSscale, current_pairs ) :
							get_probability( R, current_pairs );

	int RNG_xDir 	= paraRdr->getVal("RNG_xDir");
	int RNG_yDir 	= paraRdr->getVal("RNG_yDir");
	int RNG_zDir 	= paraRdr->getVal("RNG_zDir");

	constexpr bool check_number_of_shifted_particles = false;
	int iLoop = 0;
	int nLoops = paraRdr->getVal("shifter_nLoops");
	for (iLoop = 0; iLoop < nLoops; iLoop++)
	{
		int number_of_shifted_particles = 0;

		// loop over particles, re-sample one at a time
		for (int iParticle = 0; iParticle < number_of_particles; iParticle++)
		{
			// generate a shifted momentum
			double x1 = allParticles[iParticle].p.px();
			double y1 = allParticles[iParticle].p.py();
			double z1 = allParticles[iParticle].p.pz();
			double x2 = RNG_xDir ? RNG_p0 * normal(generator) : 0.0;	// corresponds to choice of parameters in random_events.h
			double y2 = RNG_yDir ? RNG_p0 * normal(generator) : 0.0;	// corresponds to choice of parameters in random_events.h
			double z2 = RNG_zDir ? RNG_p0 * normal(generator) : 0.0;	// corresponds to choice of parameters in random_events.h

			// compute shifted configuration
			vector<ParticleRecord> allParticles_with_shift = allParticles;
			allParticles_with_shift[iParticle].p.px(x2);
			allParticles_with_shift[iParticle].p.py(y2);
			allParticles_with_shift[iParticle].p.pz(z2);
			vector<vector<double>> shifted_pairs = get_shifted_pairs( current_pairs,
																			allParticles_with_shift, iParticle );

			// get probability of shifted configuration
			double P2 = (SHIFT_MODE == "RMSscale") ?
									get_probability( RMSscale, shifted_pairs ) :
									get_probability( R, shifted_pairs );

			std::cout << "CHECK: " << SHIFT_MODE << "   "
								<< x1 << "   " << y1 << "   " << z1 << "   "
								<< x2 << "   " << y2 << "   " << z2 << "   "
								<< P2/P1 << std::endl;
			// for (const auto & p: allParticles) std::cout << p.p.pz() << "   ";
			// std::cout << std::endl;
			// for (const auto & p: allParticles_with_shift) std::cout << p.p.pz() << "   ";
			// std::cout << std::endl;

			// choose new configuration (shifted or original)
			double log_alpha = std::min(0.0, log(P2/P1));
			bool shift_this_particle = (log(uniform(generator)) <= log_alpha);

			// re-set shifted quantities
			if (shift_this_particle)
			{
				current_pairs = shifted_pairs;		// re-set current configuration
				P1 = P2;													// re-set probability of current configuration
				number_of_shifted_particles++;
				allParticles[iParticle].p.px(x2);	// re-set momentum of particle
				allParticles[iParticle].p.py(y2);	// re-set momentum of particle
				allParticles[iParticle].p.pz(z2);	// re-set momentum of particle

				// sort by pz
				// std::sort(allParticles.begin(), allParticles.end(), particleSort);
			}
		}
		if ( check_number_of_shifted_particles && number_of_shifted_particles == 0 ) break;
		std::cout << std::endl << "Before: ";
		for (const auto & p: allParticles_Original) std::cout << p.p.pz() << "   ";
		std::cout << std::endl << "After: ";
		for (const auto & p: allParticles) std::cout << p.p.pz() << "   ";
		std::cout << std::endl;
		if (true) std::terminate();
	}

	// Done.
	return;

}

		//--------------------------------
		vector<vector<double>> shifter::get_shifted_pairs( const vector<vector<double>> & pairs,
		 																					 const vector<ParticleRecord> & particles,
																							 const int shifted_particle_index )
		{
			vector<vector<double>> result = pairs;
			// auto indexer = [nCols](size_t i, size_t j) { return i*nCols+j; };

			const int number_of_particles = particles.size();

			int iPair = 0;
			for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
			for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
			{
				if (i1==shifted_particle_index || i2==shifted_particle_index)
					result[iPair] = vector<double>({
														particles[i1].p.px() - particles[i2].p.px(),
														particles[i1].p.py() - particles[i2].p.py(),
														particles[i1].p.pz() - particles[i2].p.pz()});
				iPair++;
			}

			return result;
		}

	//--------------------------------
	double shifter::get_RMSscale( const vector<ParticleRecord> & particles )
	{
		double result = 0.0;
		const int number_of_particles = particles.size();
		const int number_of_pairs = number_of_particles*(number_of_particles-1)/2;

		for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
		for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
		{
			const double dz = particles[i1].x.pz() - particles[i2].x.pz();
			result += dz*dz;
		}

		return sqrt(0.5*result/number_of_pairs);
	}

}	//End of namespace

//End of file
