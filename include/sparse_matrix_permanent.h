#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "stopwatch.h"

using namespace std;

namespace Sparse
{
  using ll      = long long;
  using Element = std::tuple<ll, ll, long double>;

  class Matrix
  {
  //==============================================================================
  private:
    ll matrix_size = 0;
    vector<Element> v;
    vector<ll> rowsums, colsums;

  //==============================================================================
  public:
    Matrix(const vector<Element> & v_in, ll matrix_size_in)
    : v{v_in}, matrix_size{matrix_size_in}
    {
      // cout << "-------------------------------------------------------------------------------\n";
      // cout << "Created Matrix object containing following elements:";
      // for (const auto & e: v) cout << " (" << get<0>(e) << "," << get<1>(e) << ")";
      // cout << "\n";
      // cout << "get_matrix_size() = " << get_matrix_size() << "\n";
      // cout << "-------------------------------------------------------------------------------\n";

      // count number of non-zero elements in each row and column
      rowsums.resize(get_matrix_size(), 0);
      colsums.resize(get_matrix_size(), 0);
      for (const auto & e: v)
      {
        // cout << "get<0>(e) = " << get<0>(e) << endl;
        // cout << "get<1>(e) = " << get<1>(e) << endl;
        // cout << "----------------------------------------" << endl;
        ++rowsums.at(get<0>(e));
        ++colsums.at(get<1>(e));
      }
    }

    //----------------------------------------------------------------------------
    Matrix get_minor(ll i, ll j) const
    {
      // select all elements not in i-th row or j-th column
      auto isNotInRowOrColumn = [i,j] (Element e)
                                { return !(get<0>(e) == i) && !(get<1>(e) == j); };
      const ll n = get_nonzero_size();
      vector<Element> vCopy(v.size());
      auto it = std::copy_if(v.begin(), v.end(), vCopy.begin(), isNotInRowOrColumn );
      vCopy.resize(std::distance(vCopy.begin(), it));

      // shift any affected indices to make iterated indexing easier
      for (auto & e: vCopy)
      {
        if (get<0>(e) > i) --get<0>(e);
        if (get<1>(e) > j) --get<1>(e);
      }
      return Matrix(vCopy, matrix_size - 1);
    }

    inline ll get_nonzero_size() const { return v.size(); }
    inline ll get_matrix_size() const { return matrix_size; }

    friend long double permanent(const Matrix & sm);

  };


  //----------------------------------------------------------------------------
  long double permanent(const Matrix & sm)
  {
    if (sm.get_nonzero_size() == 0)
      return 0.0;
    else
    {
      // check for zero rows or columns  ==>> permanent automatically zero
      for (ll index = 0; index < sm.get_matrix_size(); ++index)
        if (sm.rowsums.at(index) == 0 || sm.colsums.at(index) == 0)
          return 0.0;

      // otherwise, check if only one element remaining
      if (sm.get_nonzero_size() == 1 && get<0>(sm.v[0]) == 0)
        return get<2>(sm.v[0]);
      else
      {
        double sum = 0.0;
        for (const auto & e: sm.v)
          if (get<0>(e) == 0)
          {
            // cout << "===============================================================\n";
            // cout << "Getting minor for e = (" << get<0>(e) << "," << get<1>(e) << ")\n";
            sum += get<2>(e) * permanent( sm.get_minor( get<0>(e), get<1>(e) ) );
            // cout << "Check: sum = " << sum << "\n";
            // cout << "Got minor for e = (" << get<0>(e) << "," << get<1>(e) << ")\n";
            // cout << "===============================================================\n";
          }
        return sum;
      }
    }
  }
}
