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

  	double permanent(const vector<double> & A, long n) // expects n by n matrix encoded as vector
  	{
  		// std::cout << "\nA(exact) =\n" << fixed << setprecision(6);
  		// for (int i = 0; i < n; ++i)
  		// {
  		// 	for (int j = 0; j < n; ++j)
  		// 		std::cout << "  " << A[i*n+j];
  		// 	cout << "\n";
  		// }
  		// cout << "\n";

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


    //--------------------------------------------------------------------------
    double get_probability_Exact( const vector<Particle> & particles,
                                  const vector<vector<double>> & qVec,
                                  const vector<double> & BE_distances )
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
  // cout << "Check q: " << q[0] << "  " << q[1] << "  " << q[2] << "  " << q2 << "  " << tmp << "\n";
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
    double get_probability_FullProduct( const vector<Particle> & particles,
                                        const vector<vector<double>> & qVec,
                                        const vector<double> & BE_distances )
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
    double get_probability_Speed( const vector<Particle> & particles,
                                  const vector<vector<double>> & qVec,
                                  const vector<double> & BE_distances )
  	{
  		// use only np-1 independent pairs, and cycle over which gets omitted
  		const int n = qVec.size();
  		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
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
                                        const vector<vector<double>> & qVec,
                                        const vector<double> & BE_distances )
    {
      return mp.evaluate(particles, BE_distances);
    }


  //----------------------------------------------------------------------------
  public:
    ConfigurationProbability( const string & mode,
                      				param_list & parameters,
                              ParameterReader * paraRdr_in,
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
				get_probability = [this]( const vector<Particle> & particles,
                                  const vector<vector<double>> & qVec,
                                  const vector<double> & BE_distances )
                          { return get_probability_AlmostExact(
                                    particles, qVec, BE_distances); };
        mp = MatrixPermanent(n_particles, precision, assume_sparse);
			}
			else
			{
				err << "Invalid mode = " << mode << endl;
				std::terminate();
			}
		}
    ~ConfigurationProbability(){}

    std::function<double( const vector<Particle> &,
                          const vector<vector<double>> &,
                          const vector<double> &)> get_probability;

};

#endif
