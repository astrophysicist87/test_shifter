#ifndef PERMANENT_H
#define PERMANENT_H

#include <cassert>
#include <functional>
#include <vector>

#include "FourVector.h"
#include "Particle.h"
#include "stopwatch.h"

using namespace std;
using shift_lib::Vec4;
using shift_lib::Particle;

class MatrixPermanent
{
  typedef struct
	{
		int clusterID;
		double permanent;
		std::vector<Particle> v;
	} Cluster;

  typedef struct
	{
		int ID1, ID2;
    Particle p1, p2;
		std::vector<double> q;
	} Pair;

  //============================================================================
  private:
    //--------------------------------------------------------------------------
    static constexpr bool VERBOSE = false;
    static constexpr bool APPROXIMATE_LARGE_N = false;
    static constexpr double R = 5.0/0.19733;
    static constexpr int CUTOFF = 25;
    static constexpr int MAX_ROWSUM = 10000;
    bool ASSUME_SPARSE = false;
    double TINY = 1e-3;
    string PREFIX = "";
    int SAVE_shifted_particle_index = -1;

    vector<double> Cvec;

    // use these for recycling quantities
    bool clusters_set = false;
    vector<int> cluster_of_particle;
    vector<Cluster> clusters;
    vector<bool> place_particle_in_cluster;
    vector<double> permanent_by_cluster;

    //--------------------------------------------------------------------------
    void print_clusters( const vector<Cluster> & clusters_to_print, bool print_particles = false );
    void print_matrix( const vector<double> & A, long n );
    //--------------------------------------------------------------------------
    vector<Pair> get_pairs( const vector<Particle> & particles );
    double permanent_RNW( const vector<double> & A, const long long n );
    //--------------------------------------------------------------------------
    vector<double> get_A( const vector<Pair> & pairs, const int np,
                          const vector<double> & BE_distance );
    //--------------------------------------------------------------------------
    void set_clusters_with_merging( const vector<Particle> & particles,
                                    const vector<double> & BE_distance );


    //--------------------------------------------------------------------------
    double compute_permanent_from_cluster( const vector<Particle> & clusterList,
                                           const vector<double> & BE_distance );
    double get_full_product_from_pairs( const vector<Pair> & pairs );
    //--------------------------------------------------------------------------
    double permanent_by_decomposition( const vector<Particle> & particles,
                                       const vector<double> & BE_distance );
    void remove_shifted_cluster( int shifted_particle_index );
    //--------------------------------------------------------------------------
    double sparse_permanent( const vector<double> & A, long n );
    vector<double> take_minor( const vector<double> & A, long i0, long j0, long n );
    double permanent_by_expansion( const vector<double> & A, long n );
    vector<long> get_rowsums( const vector<double> & A, long n );


  //============================================================================
  public:
    MatrixPermanent(){}
    MatrixPermanent( const int n, double TOLERANCE, bool ASSUME_SPARSE_IN, const string & PREFIX_IN )
    : TINY{TOLERANCE},
      ASSUME_SPARSE{ASSUME_SPARSE_IN},
      PREFIX{PREFIX_IN}
    {
      Cvec.resize(n+1);
      for (long i = 0; i <= n; i++)
        Cvec[i] = (double)pow((double)2, i);
      cluster_of_particle.resize(n, -1);
      place_particle_in_cluster.resize(n, true);
    }


    double evaluate_approximate_permanent( const vector<Particle> & particles,
                        const vector<double> & BE_distance,
                        const int shifted_particle_index = -1 );

    double evaluate_exact_permanent( const vector<Particle> & particles,
                        const vector<double> & BE_distance,
                        const int shifted_particle_index = -1 );

    double evaluate_full_product( const vector<Particle> & particles,
                        const vector<double> & BE_distance,
                        const int shifted_particle_index = -1 );

    // resets internal state of MatrixPermanent object if proposed shift rejected
    void revert_state( const vector<Particle> & particles,
                       const vector<double> & BE_distance );

};


#endif
