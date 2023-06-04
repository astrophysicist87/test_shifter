#ifndef PERMANENT_H
#define PERMANENT_H

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
  template <typename T>
    using vv = std::vector<std::vector<T>>;

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
    static constexpr int CUTOFF = 1;
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
    /*{
    	std::cout << "\n\n";
      std::cout << "--------------------------------------------------------------------------------\n"
                << "PRINTING CLUSTERS:\n"
                << "Clusters (size = " << clusters_to_print.size() << "):\n" << setprecision(4);
      int iCluster = 0;
    	for (const auto & cluster: clusters_to_print)
    	{
        std::cout << PREFIX << "  - Cluster " << iCluster << " (size = " << cluster.v.size() << "):";
    		if (print_particles)
          for (const auto & node: cluster.v)
    			 std::cout << "  " << node.particleID;
    		std::cout << "\n";
        iCluster++;
    	}
      std::cout << "--------------------------------------------------------------------------------\n\n";
    }*/

    //--------------------------------------------------------------------------
    vector<Pair> get_pairs( const vector<Particle> & particles )
    {
    	vector<Pair> result;

    	const int n = particles.size();

    	// get all values in vector first
    	for (int i1 = 0; i1 < n - 1; ++i1)
    	for (int i2 = i1 + 1; i2 < n; ++i2)
      {
    		Vec4 q = particles[i1].p - particles[i2].p;
    		result.push_back( Pair({ particles[i1].particleID,
                                 particles[i2].particleID,
                                 particles[i1], particles[i2],
                                 vector<double>({q.x(), q.y(), q.z()}) }) );
      }

    	return result;
    }


    //--------------------------------------------------------------------------
    // inline vector<long> dec2binarr(long n, long dim)
    inline void dec2binarr(long n, long dim, vector<long> & res)
    {
        // note: res[dim] will save the sum res[0]+...+res[dim-1]
        // long* res = (long*)calloc(dim + 1, sizeof(long));
        // vector<long> res(dim+1);
        long pos = dim - 1;

        // note: this will crash if dim < log_2(n)...
        while (n > 0)
        {
            res[pos] = n % 2;
            res[dim] += res[pos];
            n /= 2; // integer division
            pos--;
        }

        return;
    }

    void print_matrix( const vector<double> & A, long n );
    /*{
      std::cout << "\nn = " << n << "\n";
      std::cout << "\nA(almost exact) =\n" << fixed << setprecision(4);
      for (int i = 0; i < n; ++i)
      {
        for (int j = 0; j < n; ++j)
          std::cout << " " << A[i*n+j];
        cout << "\n";
      }
      cout << "\n";
    }*/


    //--------------------------------------------------------------------------
    // expects n by n matrix encoded as vector
    double permanent( const vector<double> & A, long n )
    {
      if (VERBOSE)
        print_matrix(A, n);

      double sum = 0.0;
      double rowsumprod = 0.0, rowsum = 0.0;
      vector<long> chi(n + 1);
      // double C = (double)pow((double)2, n);
      double C = Cvec[n];

      // loop all 2^n submatrices of A
      for (long k = 1; k < C; k++)
      {
        rowsumprod = 1.0;
        std::fill( chi.begin(), chi.end(), 0 );
        dec2binarr(k, n, chi); // characteristic vector

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
          if (ASSUME_SPARSE && rowsumprod < TINY) break;
        }

        sum += (double)pow((double)-1, n - chi[n]) * rowsumprod;
      }

      return sum;
    }


    //--------------------------------------------------------------------------
    vector<double> get_A( const vector<Pair> & pairs, const int np,
                          const vector<double> & BE_distance )
    {
			auto UTindexer = [](int i, int j, int n){return -1 + j - i*(3 + i - 2*n)/2;};
      int index = 0;
      vector<double> A(np*np);
      int total_n = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*BE_distance.size())));
      for (int i = 0; i < np; i++)
      {
        A[i*np+i] = 1.0;
        for (int j = i+1; j < np; j++)
        {
          auto & pair = pairs[index];
          const auto & q = pair.q;
          double q2 = inner_product(q.cbegin(), q.cend(),
                                    q.cbegin(),
                                    0.0);
          double tmp = exp(-0.25*q2*R*R);
// if (abs(BE_distance[UTindexer(pair.ID1, pair.ID2, np)]-tmp) > 1e-10)
// {
//   cout << "q: " << q[0] << "  " << q[1] << "  " << q[2] << "\n\n";
//   cout << "Particle 1: " << pair.p1 << "\n";
//   cout << "Particle 2: " << pair.p2 << "\n";
//   cout << tmp << "  " << UTindexer(pair.ID1, pair.ID2, total_n)
//       << "  " << BE_distance[UTindexer(pair.ID1, pair.ID2, total_n)] << "  "
//       << abs(tmp-BE_distance[UTindexer(pair.ID1, pair.ID2, total_n)]) << "\n\n";
//   cout << "CHECK: " << pair.ID1 << "  " << pair.ID2 << "  " << BE_distance[UTindexer(pair.ID1, pair.ID2, total_n)] << endl;
//   std::terminate();
// }
          if (tmp < TINY) tmp = 0.0;	// make matrix as sparse as possible
          A[i*np+j] = tmp;
          A[j*np+i] = tmp; // matrix is symmetric
          index++;
        }
      }
      return A;
    }

//--------------------------------------------------------------------------
    inline double get_q2( const Vec4 & q )
    {
      return q.x()*q.x() + q.y()*q.y() + q.z()*q.z();
    }

    //--------------------------------------------------------------------------
    inline double get_BE_distance( const Particle & p1,
                                   const Particle & p2 )
    {
    	return exp(-0.25*get_q2(p1.p - p2.p)*R*R);
    }


    //--------------------------------------------------------------------------
    void set_clusters_with_merging( const vector<Particle> & particles,
                                    const vector<double> & BE_distance )
    {
			auto UTindexer = [](int i, int j, int n){return -1 + j - i*(3 + i - 2*n)/2;};

      int iParticle = 0;
      int number_of_particles = particles.size();
      for (const auto & particle: particles)
      {

        // if this particle already belongs to a cluster, continue
        if ( !place_particle_in_cluster[ particle.particleID ] )
          continue;

        int iCluster = 0;
        vector<int> indices;
        for (const auto & cluster: clusters)  // loop over all clusters
        {
          // check if particle is close to any node (particle) in this cluster
          for (const auto & node: cluster.v)
          {
            auto [i,j] = std::minmax(node.particleID, particle.particleID);

            // if ( get_BE_distance(particle, node) > TINY )
            if ( BE_distance[UTindexer( i, j, number_of_particles )] > TINY )
            {
              indices.push_back(iCluster);
              break;  // no need to check other nodes in this cluster
            }
          }
          ++iCluster;
        }

        // merge particle with any identified clusters
        if (indices.empty())          // add particle to new cluster
          clusters.push_back( Cluster({-1, -1.0, vector<Particle>({particle})}) );
        else if (indices.size()==1)   // add particle to unique matching cluster
        {
          auto & this_cluster = clusters[ indices[0] ];
          this_cluster.v.push_back( particle );
          this_cluster.permanent = -1.0;
        }
        else // otherwise add to appropriate clusters and merge
        {
          vector<Particle> merged_clusters{particle};
          for (auto it = indices.rbegin(); it != indices.rend(); ++it)
          {
            merged_clusters.insert( merged_clusters.end(),
                                    clusters[*it].v.begin(),
                                    clusters[*it].v.end() );
            clusters.erase (clusters.begin() + *it);
          }
          clusters.push_back( Cluster({-1, -1.0, merged_clusters }) );
        }
      }

      // if (VERBOSE)
      //   print_clusters(clusters, true);

      // set vector to track which cluster each particle belongs to
      // set also indices of cluster
      clusters_set = true;
      std::fill( cluster_of_particle.begin(),
                 cluster_of_particle.end(), -1 );
      {
        int iCluster = 0;
        for (auto & cluster: clusters)  // loop over all clusters
        {
          cluster.clusterID = iCluster;
          for (const auto & node: cluster.v)
            cluster_of_particle[ node.particleID ] = iCluster;

          // compute permanent of this cluster
          if (cluster.permanent < 0.0)
            cluster.permanent = compute_permanent_from_cluster(cluster.v, BE_distance);
          iCluster++;
        }
      }


      return;
    }


    //--------------------------------------------------------------------------
    double compute_permanent_from_cluster(
             const vector<Particle> & clusterList,
             const vector<double> & BE_distance )
    {
      //------------------------------------------------------------------------
      // find pairs
      vector<Pair> pairs = get_pairs( clusterList );

      //------------------------------------------------------------------------
      // set matrix
      vector<double> A = get_A(pairs, clusterList.size(), BE_distance);

      if (APPROXIMATE_LARGE_N && clusterList.size()>20)
      {
        // full product
        double full_product = 1.0;
        double sum1 = 1.0;
        for (const auto & pair: pairs)
        {
          double q2 = inner_product( pair.q.cbegin(), pair.q.cend(),
                                     pair.q.cbegin(), 0.0 );
          full_product *= 1.0 + exp(-0.5*q2*R*R);
          sum1 += exp(-0.5*q2*R*R);
        }

        Stopwatch sw;
        double time1 = 0.0, time2 = 0.0;
        // if (VERBOSE)
//           cout << "COMPARE: " << clusterList.size() << "  " << full_product
//                 << "  " << sum1;
//           sw.Start();
//           cout << "  " << permanent( A, clusterList.size() );
//           sw.Stop();
//           time1 = sw.printTime();
//           sw.Reset();
//           sw.Start();
//           cout << "  " << sparse_permanent( A, clusterList.size() ) << endl;
//           sw.Stop();
//           time2 = sw.printTime();
// cout << "  -> Took: " << time1 << " vs. " << time2 << endl;
// cout << "------------------------------------------------------------------------" << endl;

        return full_product;
      }

      //------------------------------------------------------------------------
      // compute and return permanent
      if (VERBOSE)
      {
        double tmp = permanent( A, clusterList.size() );
        cout << "--------------------------------------------------" << endl;
        cout << "Permanent of above A: " << tmp << endl;
        cout << "--------------------------------------------------" << endl;
        return tmp;
      }
      else
      {
// cout << "COMPARE: " << setprecision(12) << permanent( A, clusterList.size() ) << "  " << sparse_permanent( A, clusterList.size() ) << "\n";
        return sparse_permanent( A, clusterList.size() );
      }
    }


    //--------------------------------------------------------------------------
    double permanent_by_decomposition( const vector<Particle> & particles,
                                       const vector<double> & BE_distance )
    {
      // set any clusters which have changed
      set_clusters_with_merging( particles, BE_distance );

      // total permanent is product of sub-permanents
      double decomposed_permanent = 1.0;
      for (const auto & cluster: clusters)
        decomposed_permanent *= cluster.permanent;

      return decomposed_permanent;
    }


    void remove_shifted_cluster( int shifted_particle_index )
    {
      // identify cluster which contains shifted particle
      int cluster_to_remove = cluster_of_particle[ shifted_particle_index ];

      // mark each particle in cluster as needing to be re-placed into a new cluster
      // (must explicitly ignore those which are unchanged!)
      std::fill( place_particle_in_cluster.begin(),
                 place_particle_in_cluster.end(), false );
      for ( const auto & node: clusters[ cluster_to_remove ].v )
        place_particle_in_cluster[ node.particleID ] = true;

      // remove this cluster
      clusters.erase( clusters.begin() + cluster_to_remove );
    }




double sparse_permanent( const vector<double> & A, long n )
{
  // if matrix is small enough, use Ryser-Niejenhuis-Wilf algorithm
  if (n <= CUTOFF)
    return permanent( A, n );
    // otherwise, decompose into minors and try again
  else
  {
    // compute rowsums of each particle
    vector<long> rowsums = get_rowsums( A, n );

    // sort particles by increasing rowsums
    std::vector<long> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&](int i, int j) -> bool { return rowsums[i] < rowsums[j]; });

// print_matrix(A, n);
//
// cout << "\n\nRowsums:";
// for (int i = 0; i < n; ++i) cout << "  " << rowsums[i];
// cout << endl;

    // re-compute sorted matrix
    vector<double> A_sorted(n*n);
    for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      A_sorted[i*n+j] = A[indices[i]*n+indices[j]];

// print_matrix(A_sorted, n);

// rowsums = get_rowsums( A_sorted, n );

// cout << "\n\nRowsums (sorted):";
// for (int i = 0; i < n; ++i) cout << "  " << rowsums[i];
// cout << endl;
//
// std::terminate();

    // expand low-rowsum particles in minors until RNW is viable
    return permanent_by_expansion( A_sorted, n );
  }
}

vector<double> take_minor( const vector<double> & A, long i0, long j0, long n )
{
  vector<double> A_minor((n-1)*(n-1));
  long index = 0;
  for (int i = 0; i < n; ++i)
  for (int j = 0; j < n; ++j)
    if (i!=i0 && j!=j0)
      A_minor[index++] = A[i*n+j];
  return A_minor;
}


double permanent_by_expansion( const vector<double> & A, long n )
{
  // if matrix is small enough, use Ryser-Niejenhuis-Wilf algorithm
  if (n <= CUTOFF)
    return permanent( A, n );
  // otherwise, decompose into minors and try again
  else
  {
    double result = 0.0;

    if (n == 2)
      return A[0*n+0]*A[1*n+1]+A[0*n+1]*A[1*n+0];

    // loop over non-zero column entries in 0th row
    constexpr long i = 0;  // 0th row
    for (long j = 0; j < n; ++j)
    {
      if (abs(A[i*n+j]) > TINY)
        result += A[i*n+j] * permanent_by_expansion( take_minor(A, i, j, n), n-1 );
    }
    return result;
  }
}


vector<long> get_rowsums( const vector<double> & A, long n )
{
  vector<long> rowsums(n);
  for (int i = 0; i < n; ++i)
  for (int j = 0; j < n; ++j)
    if (A[i*n+j] > TINY)
      rowsums[i]++;
  return rowsums;
}

  //============================================================================
  public:
    MatrixPermanent(){}
    ~MatrixPermanent(){}
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

    inline double evaluate( const vector<Particle> & particles,
                            const vector<double> & BE_distance,
                            const int shifted_particle_index = -1 )
    {
      // allows to make revert_state consistent if needed
      SAVE_shifted_particle_index = shifted_particle_index;

      // if no/negative particle index passed, get permanent of all particles
      if ( shifted_particle_index < 0 )
      {
        clusters.clear();
        std::fill( place_particle_in_cluster.begin(),
                   place_particle_in_cluster.end(), true );
      }
      else
        remove_shifted_cluster( shifted_particle_index );

      // permanent_by_decomposition as usual either way
      // (clusters get re-set inside)
      return permanent_by_decomposition(particles, BE_distance);
    }

    // resets internal state of MatrixPermanent object if proposed shift rejected
    void revert_state( const vector<Particle> & particles,
                       const vector<double> & BE_distance )
    {
      // negative index means nothing to do,
      // clusters were completely reset and nothing was saved
      if (SAVE_shifted_particle_index < 0)
        return;

      // remove the cluster containing the shifted particle
      remove_shifted_cluster( SAVE_shifted_particle_index );

      // re-distribute removed cluster amongst remaining clusters
      set_clusters_with_merging( particles, BE_distance );

      return;
    }

};


#endif
