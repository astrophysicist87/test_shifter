#ifndef PROBABILITY_H
#define PROBABILITY_H

#include "param_list.h"
#include "ParameterReader.h"
#include "permanent.h"

using shift_lib::ParameterReader;

//==============================================================================
class ConfigurationProbability
{
  //----------------------------------------------------------------------------
  private:
    double R { 0.0 };
    MatrixPermanent mp;
		ostream & out;
		ostream & err;
    ParameterReader * paraRdr;

    static constexpr double TINY = 1e-6;

    //--------------------------------------------------------------------------
    double get_probability_Exact( const vector<Particle> & particles,
                                  // const vector<vector<double>> & qVec,
                                  const vector<double> & BE_distances,
                                  const int shifted_particle_index )
    {
      return mp.evaluate_exact_permanent( particles, BE_distances,
                                          shifted_particle_index );
  	}

  	//--------------------------------------------------------------------------
    // sorted in MP as this gives a mediocre approximation of permanent
    double get_probability_FullProduct( const vector<Particle> & particles,
                                        // const vector<vector<double>> & qVec,
                                        const vector<double> & BE_distances,
                                        const int shifted_particle_index )
  	{
      return mp.evaluate_full_product( particles, BE_distances,
                                       shifted_particle_index );
  	}

  	//--------------------------------------------------------------------------
    double get_probability_Speed( const vector<Particle> & particles,
                                  const vector<double> & BE_distances,
                                  const int shifted_particle_index )
  	{
  		// use only np-1 independent pairs, and cycle over which gets omitted
      const int np = particles.size();
  		double normalization = paraRdr->getVal("shifter_norm");

  		auto square = [](double x){return x*x;};
  		auto UTindexer = [](int i, int j, int n){return -1 + j - i*(3 + i - 2*n)/2;};
  		double init = 1.0 + 0.5 * np * normalization
                          * square( BE_distances[ UTindexer(0, np-1, np) ] );
  		double result = init;
  		double factor = 1.0/init;

  		// update with other pairs
  		for (int i1 = 0; i1 < np - 1; ++i1)
  		{
  			auto i = UTindexer(i1, i1+1, np);
  			double term = 1.0 + 0.5 * np * normalization * square( BE_distances[i] );
  			result *= term;
  			factor += 1.0/term;
  		}

  		return factor*result/np;
  	}

  	//--------------------------------------------------------------------------
    double get_probability_AlmostExact( const vector<Particle> & particles,
                                        // const vector<vector<double>> & qVec,
                                        const vector<double> & BE_distances,
                                        const int shifted_particle_index )
    {
      return mp.evaluate_approximate_permanent( particles, BE_distances,
                                                shifted_particle_index );
    }


  //----------------------------------------------------------------------------
  public:
    ConfigurationProbability( const string & mode,
                      				param_list & parameters,
                              ParameterReader * paraRdr_in,
                              const string & prefix,
                      				ostream & out_stream = std::cout,
                      				ostream & err_stream = std::cerr )
		: paraRdr{paraRdr_in},
      out{out_stream},
			err{err_stream}
		{
			if (mode == "AlmostExact")
			{
				bool assume_sparse = std::get<bool>(parameters.at("assume_sparse"));
				int n_particles    = std::get<int>(parameters.at("n_particles"));
				double precision   = std::get<double>(parameters.at("precision"));

        // set function for computing probability of given configuration
				get_probability = [this]( const vector<Particle> & particles,
                                  const vector<double> & BE_distances,
                                  const int shifted_particle_index )
                          { return get_probability_AlmostExact(
                                      particles, BE_distances,
                                      shifted_particle_index ); };

        // initialize MatrixPermanent mp object and pass to revert_state lambda
        mp = MatrixPermanent(n_particles, precision, assume_sparse, prefix);
        revert_state = [this]( const vector<Particle> & particles,
                               const vector<double> & BE_distances )
                       { mp.revert_state( particles, BE_distances ); };
			}
      else if (mode == "Exact")
      {
        bool assume_sparse = std::get<bool>(parameters.at("assume_sparse"));
        int n_particles    = std::get<int>(parameters.at("n_particles"));
        double precision   = std::get<double>(parameters.at("precision"));

        // set function for computing probability of given configuration
        get_probability = [this]( const vector<Particle> & particles,
                                  const vector<double> & BE_distances,
                                  const int shifted_particle_index )
                          { return get_probability_Exact(
                                      particles, BE_distances,
                                      shifted_particle_index ); };

        // initialize MatrixPermanent mp object and pass to revert_state lambda
        mp = MatrixPermanent(n_particles, precision, assume_sparse, prefix);
        revert_state = [this]( const vector<Particle> & particles,
                               const vector<double> & BE_distances )
                       { mp.revert_state( particles, BE_distances ); };
      }
      else if (mode == "FullProduct")
      {
        bool assume_sparse = std::get<bool>(parameters.at("assume_sparse"));
        int n_particles    = std::get<int>(parameters.at("n_particles"));
        double precision   = std::get<double>(parameters.at("precision"));

        // set function for computing probability of given configuration
        get_probability = [this]( const vector<Particle> & particles,
                                  const vector<double> & BE_distances,
                                  const int shifted_particle_index )
                          { return get_probability_FullProduct(
                                    particles, BE_distances,
                                    shifted_particle_index ); };

        // initialize MatrixPermanent mp object and pass to revert_state lambda
        mp = MatrixPermanent(n_particles, precision, assume_sparse, prefix);
        revert_state = [this]( const vector<Particle> & particles,
                               const vector<double> & BE_distances )
                       { mp.revert_state( particles, BE_distances ); };
      }
      else if (mode == "Speed")
      {
        // set function for computing probability of given configuration
        get_probability = [this]( const vector<Particle> & particles,
                                  const vector<double> & BE_distances,
                                  const int shifted_particle_index )
                          { return get_probability_FullProduct(
                                    particles, BE_distances,
                                    shifted_particle_index ); };

        // revert_state does nothing in this case
        revert_state = [this]( const vector<Particle> & particles,
                               const vector<double> & BE_distances ) { ; };
      }
      else
			{
				err << "Invalid mode = " << mode << endl;
				std::terminate();
			}
		}
    // ~ConfigurationProbability(){}

    double get_probability( const vector<Particle> &,
                          const vector<double> &,
                          const int, const string & chosen_mode )
    {
      if (chosen_mode == "AlmostExact")
        return get_probability_AlmostExact( particles, BE_distances, shifted_particle_index );
      else if (chosen_mode == "Exact")
        return get_probability_Exact( particles, BE_distances, shifted_particle_index );
      else if (chosen_mode == "FullProduct")
        return get_probability_FullProduct( particles, BE_distances, shifted_particle_index );
      else if (chosen_mode == "Speed")
        return get_probability_Speed( particles, BE_distances, shifted_particle_index );
      else
			{
				err << "Invalid mode = " << mode << endl;
				std::terminate();
			}
    }

    std::function<double( const vector<Particle> &,
                          const vector<double> &,
                          const int )> get_probability;

    std::function<void( const vector<Particle> &,
                        const vector<double> & )> revert_state;
};

#endif
