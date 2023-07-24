#include <algorithm>
#include <complex>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "ParameterReader.h"
#include "stopwatch.h"

using namespace std;

typedef struct
{
  int ID;
  vector<double> p, x;
} Particle;

constexpr double TINY = 1e-6;
constexpr double HBARC = 0.19733;
//==============================================================================
class QuantumCorrelator
{
  //----------------------------------------------------------------------------
  private:
    //--------------------------------------------------------------------------
    int RNG_xDir { 0 };
    int RNG_yDir { 0 };
    int RNG_zDir { 0 };
    int dimension { 0 };
    double sigma { 1.0 };  // fm
    double s1p_normalization { 0.0 };
    vector<double> K;
    vector<Particle> & particles;

    //--------------------------------------------------------------------------
    inline double vec_sq(const auto & v)
    {
      return inner_product(v.begin(), v.end(), v.begin(), 0.0);
    }
    inline double fourvec_dot(const vector<double> & a, const vector<double> & b)
    {
      return a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];
    }

    inline double s1p_prob(const vector<double> & pi, const vector<double> & K)
    {
      return s1p_normalization * exp(-sigma*sigma*vec_sq({ pi[1] - K[1],
                                                           pi[2] - K[2],
                                                           pi[3] - K[3] }) );
    }

  //----------------------------------------------------------------------------
  public:
    QuantumCorrelator( vector<Particle> & particles_in,
                    double sigma_in, int RNG_xDir_in, int RNG_yDir_in, int RNG_zDir_in )
    : particles{ particles_in }
    {
      sigma = sigma_in/HBARC;  // convert to inverse GeV
      RNG_xDir = RNG_xDir_in;
      RNG_yDir = RNG_yDir_in;
      RNG_zDir = RNG_zDir_in;
      dimension = RNG_xDir + RNG_yDir + RNG_zDir;
      s1p_normalization = pow(sqrt(M_PI)*sigma, dimension);
      K = vector<double>({0.5*(P1[1]+P2[1]), 0.5*(P1[2]+P2[2]), 0.5*(P1[3]+P2[3])});  // 3-vector
      q = vector<double>({P2[0]-P1[0], P2[1]-P1[1], P2[2]-P1[2], P2[3]-P1[3]});       // 4-vector
    }
    ~QuantumCorrelator(){}

    double get_spectra(const vector<double> & P)
    {
      double spectra = 0.0;
      for (const auto & p: particles) spectra += s1p_prob(p.p, P);
      return spectra;
    }
    double get_correlator(const vector<double> & P1, const vector<double> & P2)
    {
      // References:
      //   - Phys.Rept. 319 (1999) 145-230  <<--- using this version
      //   - Phys.Rev.C 56 (1997) R614-R618
      using namespace std::complex_literals;
      auto sq = [](const auto & x){ return x*x; };

      double spectra1 = 0.0, spectra2 = 0.0, Tc_num = 0.0, Tc_den = 0.0;
      complex<double> FTspectra = 0.0;
      double quantum_prefactor = exp(-0.5*sigma*sigma*vec_sq(q));
      for (const auto & p: particles)
      {
        double s1K = s1p_prob(p.p, K);
        double s1P1 = s1p_prob(p.p, P1);
        double s1P2 = s1p_prob(p.p, P2);
        spectra1 += s1P1;
        spectra2 += s1P2;
        Tc_num += s1K*s1K;
        Tc_den += s1P1*s1P2;
        FTspectra += s1K*exp((1i)*fourvec_dot(q,p.x));
      }
      double numerator = sq(abs(FTspectra)) - Tc_num;
      double denominator = spectra1*spectra2 - Tc_den;
      return 1.0 + quantum_prefactor * numerator / denominator; // == C(q,K)
    }
};
