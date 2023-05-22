#ifndef PERMANENT_H
#define PERMANENT_H

#include <functional>
#include <vector>

#include "FourVector.h"
#include "ParticleRecord.h"

using namespace std;
using shift_lib::Vec4;
using shift_lib::ParticleRecord;

// forward declaration
// class MatrixPermanent;
// constexpr bool MatrixPermanent::VERBOSE = false;
// constexpr double MatrixPermanent::R = 5.0/0.19733;


class MatrixPermanent
{
  template <typename T>
  using vv = std::vector<std::vector<T>>;

  //============================================================================
  private:
    //--------------------------------------------------------------------------
    const bool VERBOSE = false;
    const double R = 5.0/0.19733;
    bool ASSUME_SPARSE = false;
    double TINY = 1e-3;

    vector<double> Cvec;

    //--------------------------------------------------------------------------
    void print_clusters( const vv<ParticleRecord> & clusters )
    {
    	std::cout << "\n\nClusters (size = " << clusters.size() << "):\n";
      int iCluster = 0;
    	for (const auto & cluster: clusters)
    	{
        std::cout << "Cluster " <<iCluster  << " (size = " << cluster.size() << "):";
    		for (const auto & node: cluster)
    			std::cout << "  " << node;
    		std::cout << "\n";
        iCluster++;
    	}
    }

    //--------------------------------------------------------------------------
    vv<double> get_pairs( const vector<ParticleRecord> & particles )
    {
    	vv<double> result;

    	const int n = particles.size();

    	// get all values in vector first
    	for (int i1 = 0; i1 < n - 1; ++i1)
    	for (int i2 = i1 + 1; i2 < n; ++i2)
      {
    		Vec4 q = particles[i1].p - particles[i2].p;
    		result.push_back( vector<double>({q.px(), q.py(), q.pz()}) );
      }

    	return result;
    }


    //--------------------------------------------------------------------------
    // inline void dec2binarr(long n, long dim, vector<long> & res)
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
            n /= 2; // integer division
            pos--;
        }

        return res;
    }

    //--------------------------------------------------------------------------
    // expects n by n matrix encoded as vector
    double permanent( const vector<double> & A, long n )
    {
      if (VERBOSE)
      {
        std::cout << "\nA(almost exact) =\n" << fixed << setprecision(6);
        for (int i = 0; i < n; ++i)
        {
          for (int j = 0; j < n; ++j)
            std::cout << " " << A[i*n+j];
          cout << "\n";
        }
        cout << "\n";
      }

      double sum = 0.0;
      double rowsumprod = 0.0, rowsum = 0.0;
      vector<long> chi(n + 1);
      double C = (double)pow((double)2, n);
      // double C = Cvec[n];

      // loop all 2^n submatrices of A
      for (long k = 1; k < C; k++)
      {
          rowsumprod = 1.0;
          chi = dec2binarr(k, n); // characteristic vector
          // dec2binarr(k, n, chi); // characteristic vector

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
    vector<double> get_A( const vv<double> & pairs, const int np,
                          const vector<double> & BE_distance )
    {
      int index = 0;
      vector<double> A(np*np);
      for (int i = 0; i < np; i++)
      {
        A[i*np+i] = 1.0;
        for (int j = i+1; j < np; j++)
        {
          const auto & q = pairs[index];
          double q2 = inner_product(q.cbegin(), q.cend(),
                                    q.cbegin(),
                                    0.0);
          double tmp = exp(-0.25*q2*R*R);
          // double tmp = BE_distance[index];
// cout << "check: " << tmp << "   " << BE_distance[index] << endl;
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
      return q.px()*q.px() + q.py()*q.py() + q.pz()*q.pz();
    }

    //--------------------------------------------------------------------------
    inline double get_BE_distance( const ParticleRecord & p1,
                                   const ParticleRecord & p2 )
    {
    	return exp(-0.25*get_q2(p1.p - p2.p)*R*R);
    }

    //--------------------------------------------------------------------------
    vv<ParticleRecord> get_clusters_with_merging(
                         const vector<ParticleRecord> & particles,
                         const vector<double> & BE_distance )
    {
			auto UTindexer = [](int i, int j, int n){return -1 + j - i*(3 + i - 2*n)/2;};

    	vv<ParticleRecord> clusters;
      int number_of_particles = particles.size();
      for (const auto & particle: particles)
      {
        int iCluster = 0;
        vector<int> indices;
        for (auto & cluster: clusters)  // loop over all clusters
        {
          // check if particle is close to any node (particle) in this cluster
          for (const auto & node: cluster)
          {
            if ( get_BE_distance(particle, node) > TINY )
            // if ( BE_distance[UTindexer( node.particleID,  // NOTE THE ORDER!!
            //                             particle.particleID,
            //                             number_of_particles )] > TINY )
            {
              indices.push_back(iCluster);
              break;  // no need to check other nodes in this cluster
            }
          }
          iCluster++;
        }

        if (indices.empty())          // add particle to new cluster
          clusters.push_back( vector<ParticleRecord>({particle}) );
        else if (indices.size()==1)   // add particle to unique matching cluster
          clusters[ indices[0] ].push_back( particle );
        else // otherwise add to appropriate clusters and merge
        {
          vector<ParticleRecord> merged_clusters{particle};
          for (auto it = indices.rbegin(); it != indices.rend(); ++it)
          {
            merged_clusters.insert( merged_clusters.end(),
                                    clusters[*it].begin(),
                                    clusters[*it].end() );
            clusters.erase (clusters.begin() + *it);
          }
          clusters.push_back( merged_clusters );
        }

        // iParticle++;
      }

      // print_clusters(clusters);

      return clusters;
    }

    //--------------------------------------------------------------------------
    double compute_permanent_from_cluster(
             const vector<ParticleRecord> & cluster,
             const vector<double> & BE_distance )
    {
      //------------------------------------------------------------------------
      // find pairs
      vv<double> pairs = get_pairs( cluster );

      //------------------------------------------------------------------------
      // set matrix
      vector<double> A = get_A(pairs, cluster.size(), BE_distance);

      //------------------------------------------------------------------------
      // compute and return permanent
      return permanent( A, cluster.size() );
    }


    //--------------------------------------------------------------------------
    double permanent_by_decomposition( const vector<ParticleRecord> & particles,
                                       const vector<double> & BE_distance )
    {
      double decomposed_permament = 1.0;
      for (const auto & cluster: get_clusters_with_merging(particles, BE_distance))
        decomposed_permament *= compute_permanent_from_cluster(cluster, BE_distance);
      return decomposed_permament;
    }

  //============================================================================
  public:
    MatrixPermanent( const int n, const double TOLERANCE = 1e-6,
                     const bool ASSUME_SPARSE_IN = false)
    : TINY{TOLERANCE},
      ASSUME_SPARSE{ASSUME_SPARSE_IN}
    {
      Cvec.resize(n+1);
      for (long i = 0; i < n; i++)
        Cvec[i] = (double)pow((double)2, i);
    }
    ~MatrixPermanent(){}
    inline double evaluate( const vector<ParticleRecord> & particles,
                            const vector<double> & BE_distance )
                  { return permanent_by_decomposition(particles, BE_distance); }

};


#endif
