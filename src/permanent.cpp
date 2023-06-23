#include <algorithm>
#include <climits>
#include <functional>
#include <numeric>
#include <vector>

#include "../include/permanent.h"
#include "../include/sparse_matrix_permanent.h"

using namespace std;
using shift_lib::Vec4;
using shift_lib::Particle;

//--------------------------------------------------------------------------
void MatrixPermanent::print_clusters( const vector<Cluster> & clusters_to_print, bool print_particles )
{
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
}

//------------------------------------------------------------------------------
void MatrixPermanent::print_matrix( const vector<double> & A, long n )
{
  std::cout << "\nn = " << n << "\n";
  std::cout << "\nA(almost exact) =\n" << fixed << setprecision(4);
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
      std::cout << " " << A[i*n+j];
    cout << "\n";
  }
  cout << "\n";
}


//--------------------------------------------------------------------------
vector<MatrixPermanent::Pair> MatrixPermanent::get_pairs( const vector<Particle> & particles )
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
// expects n by n matrix encoded as vector
double MatrixPermanent::permanent_RNW( const vector<double> & A, const long long n )
{
	assert((A.size() == n*n) && "A must be an n x n matrix!");
  if (VERBOSE)
    print_matrix(A, n);

	if (n >= 30) cout << "n = " << n << endl;

  // loop all 2^n submatrices of A
  double sum = 0.0;
  unsigned long long count = 0;
  unsigned long long C = (1ULL << n); // bitshift equals integer pow() for base2

  vector<bool> chi(n);
  vector<double> rowsums(n, 0.0);
  for ( unsigned long long k = 0; k < C - 1; ++k )
  {
    // order submatrices by gray code, identify which bit changes
    unsigned long long mask = 1ULL, index = 0;
		while (k & mask)
		{
			mask <<= 1ULL;
			++index;
		}

    // flip the bit and store the change in sign
    chi[index] = !chi[index];
    double sign = static_cast<double>(chi[index]) - static_cast<double>(!chi[index]);
    count += llround(sign); // store the current number of 1s

    // evaluate this term in Ryser's formula
    double rowsumprod = 1.0;
    for ( unsigned long long m = 0; m < n; ++m )
    {
      rowsums[m] += sign * A[(m + 1) * n - index - 1];
      rowsumprod *= rowsums[m];
    }

		// if (n>20) cout << "CHECK: " << k << "  " << C-1 << "  " << sum << "\n";

    sum += rowsumprod * ( ((n - count) % 2) ? -1.0 : 1.0 );
  }
	// if (n>20) std::terminate();
  return sum;
}



//--------------------------------------------------------------------------
vector<double> MatrixPermanent::get_A( const vector<Pair> & pairs, const int np,
                      const vector<double> & BE_distance )
{
	// linear indexer for upper-triangular matrix
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
			auto [minID, maxID] = std::minmax(pair.ID1, pair.ID2);
			double tmp = BE_distance[UTindexer(minID, maxID, total_n)];
      if (tmp < TINY) tmp = 0.0;	// make matrix as sparse as possible
      A[i*np+j] = tmp;
      A[j*np+i] = tmp; // matrix is symmetric
      index++;
    }
  }
  return A;
}


//--------------------------------------------------------------------------
void MatrixPermanent::set_clusters_with_merging( const vector<Particle> & particles,
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
    // print_clusters(clusters, true);

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
double MatrixPermanent::compute_permanent_from_cluster(
         const vector<Particle> & clusterList,
         const vector<double> & BE_distance )
{
  //------------------------------------------------------------------------
  // find pairs
  vector<Pair> pairs = get_pairs( clusterList );

  //------------------------------------------------------------------------
  // approximate this cluster with full product, if necessary
	const int n = clusterList.size();
  if (APPROXIMATE_LARGE_N && n >= CUTOFF)
    return get_full_product_from_pairs( pairs );

  //------------------------------------------------------------------------
  // set matrix
  vector<double> A = get_A(pairs, n, BE_distance);

  //------------------------------------------------------------------------
  // compute and return permanent
	if (n < CUTOFF)
  	return permanent_RNW( A, n );
	else
	{
		// compute rowsums of each particle
		vector<long> rowsums = get_rowsums( A, n );

		// sort particles by increasing rowsums
		std::vector<long> indices(n);
		std::iota(indices.begin(), indices.end(), 0);
		std::sort(indices.begin(), indices.end(),
							[&](int i, int j) -> bool { return rowsums[i] < rowsums[j]; });

		// re-compute sorted matrix
		vector<double> A_sorted(n*n);
		for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			A_sorted[i*n+j] = A[indices[i]*n+indices[j]];

		A = std::move(A_sorted);

		vector<Sparse::Element> A_sparse;
		for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (A[i*n+j] > TINY)
				A_sparse.push_back( Sparse::Element( {i, j, A[i*n+j]} ) );
		return Sparse::permanent( Sparse::Matrix(A_sparse, n) );
	}
}

//------------------------------------------------------------------------------
double MatrixPermanent::get_full_product_from_pairs( const vector<Pair> & pairs,
                                   const vector<double> & BE_distance )
{
  // full product
  double full_product = 1.0;
	for (const auto & BEd: BE_distance)
		full_product *= 1.0 + BEd*BEd;

  return full_product;
}


//--------------------------------------------------------------------------
double MatrixPermanent::permanent_by_decomposition( const vector<Particle> & particles,
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


//--------------------------------------------------------------------------
void MatrixPermanent::remove_shifted_cluster( int shifted_particle_index )
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


//----------------------------------------------------------------------------
vector<long> MatrixPermanent::get_rowsums( const vector<double> & A, long n )
{
  vector<long> rowsums(n);
  for (int i = 0; i < n; ++i)
  for (int j = 0; j < n; ++j)
    if (A[i*n+j] > TINY)
      rowsums[i]++;
  return rowsums;
}


//----------------------------------------------------------------------------
double MatrixPermanent::evaluate_approximate_permanent(
                          const vector<Particle> & particles,
                          const vector<double> & BE_distance,
                          const int shifted_particle_index )
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


//----------------------------------------------------------------------------
double MatrixPermanent::evaluate_exact_permanent(
                          const vector<Particle> & particles,
                          const vector<double> & BE_distance,
                          const int shifted_particle_index )
{
  return compute_permanent_from_cluster( particles, BE_distance );
}

//----------------------------------------------------------------------------
double MatrixPermanent::evaluate_full_product(
                          const vector<Particle> & particles,
                          const vector<double> & BE_distance,
                          const int shifted_particle_index )
{
  return get_full_product_from_pairs( get_pairs( particles ) );
}


//----------------------------------------------------------------------------
// resets internal state of MatrixPermanent object if proposed shift rejected
void MatrixPermanent::revert_state( const vector<Particle> & particles,
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
