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

//--------------------------------------------------------------------------
// expects n by n matrix encoded as vector
long double permanent_RNW( const vector<long double> & A, const long long n )
{
	assert((A.size() == n*n) && "A must be an n x n matrix!");

  // loop all 2^n submatrices of A
  long double sum = 0.0;
  unsigned long long count = 0;
  unsigned long long C = (1ULL << n); // bitshift equals integer pow() for base2

  vector<bool> chi(n);
  vector<long double> rowsums(n, 0.0);
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
    long double sign = static_cast<long double>(chi[index]) - static_cast<long double>(!chi[index]);
    count += llround(sign); // store the current number of 1s

    // evaluate this term in Ryser's formula
    long double rowsumprod = 1.0;
    for ( unsigned long long m = 0; m < n; ++m )
    {
      rowsums[m] += sign * A[(m + 1) * n - index - 1];
      rowsumprod *= rowsums[m];
    }

    sum += rowsumprod * ( ((n - count) % 2) ? -1.0 : 1.0 );
  }
  return sum;
}


//==============================================================================
namespace Sparse
{
  using ll      = long long;
  typedef struct
  {
    ll row, col;
    long double value;
  } Element;

  class Matrix
  {
  //==============================================================================
  private:
    ll matrix_size = 0;
    vector<Element> v;
    vector<ll> rowsums, colsums;

  //==============================================================================
  public:
    Matrix(vector<Element> v_in, ll matrix_size_in)
    : v{std::move(v_in)}, matrix_size{matrix_size_in}
    {
      const ll & n = get_matrix_size();
      // cout << "v.size() = " << v.size() << endl;

      // count number of non-zero elements in each row and column
      rowsums.resize(n, 0);
      colsums.resize(n, 0);
      // cout << "Nonzero elements:";
      for (const auto & e: v)
      {
        // cout << "  (" << e.row << "," << e.col << ")";
        ++rowsums[e.row];
        // ++colsums[e.col];
      }
      // cout << endl;
      // cout << "rowsums:";
      // for (const auto r: rowsums) cout << "  " << r;
      // cout << endl;
      // cout << "colsums:";
      // for (const auto c: colsums) cout << "  " << c;
      // cout << endl;

      // sort rows by increasing rowsums
      std::vector<ll> indices(n);
      std::iota(indices.begin(), indices.end(), 0);
      // std::sort(indices.begin(), indices.end(),
      //           [&](ll i, ll j) -> bool { return rowsums[i] < rowsums[j]; });
      std::partial_sort(indices.begin(), indices.begin()+1, indices.end(),  // only need first row in place
                [&](ll i, ll j) -> bool { return rowsums[i] < rowsums[j]; });

      std::vector<ll> indices2(n);
      for (ll i = 0; i < n; ++i)
        indices2[indices[i]] = i; // inverse permutation
      indices = std::move(indices2);

      // indices = vector<ll>({3,1,4,0,2});

      // cout << "Reordering:\n";
      // for (ll i = 0; i < indices.size(); ++i)
      //   cout << "  " << i << " --> " << indices.at(i) << endl;
      // cout << endl;

      // re-compute sorted matrix
      // cout << "Sorting:\n";
      for (auto & e: v)
      {
        e.row = indices[e.row];
        e.col = indices[e.col];
      }
      std::sort(v.begin(), v.end(),
                [&](const Element & a, const Element & b) -> bool
                   { return a.row * get_matrix_size() + a.col
                          < b.row * get_matrix_size() + b.col; });

      // print_full_matrix( get_full_matrix() );

      std::fill(rowsums.begin(), rowsums.end(), 0);
      // std::fill(colsums.begin(), colsums.end(), 0);
      for (const auto & e: v)
      {
        ++rowsums[e.row];
        ++colsums[e.col];
      }

      // cout << "sorted rowsums:";
      // for (const auto r: rowsums) cout << "  " << r;
      // cout << "\n";
      // cout << "sorted colsums:";
      // for (const auto c: colsums) cout << "  " << c;
      // cout << "\n";
      // if (true) std::terminate();
    }

    //----------------------------------------------------------------------------
    Matrix get_minor(/*ll i,*/ ll j) const
    {
      // select all elements not in i-th row or j-th column
      // auto isNotInRowOrColumn = [i,j] (const Element & e)
      //                           { return !(e.row == i) && !(e.col == j); };
      auto isNotInRowOrColumn = [j] (const Element & e)
                                { return !(e.row == 0) && !(e.col == j); };
      vector<Element> vCopy(v.size());
      auto it = std::copy_if(v.begin(), v.end(), vCopy.begin(), isNotInRowOrColumn );
      vCopy.resize(std::distance(vCopy.begin(), it));

      // shift any affected indices to make iterated indexing easier
      for (auto & e: vCopy)
      {
        /*if (e.row > 0)*/ --e.row; // all rows in minor have i > 0
        if (e.col > j) --e.col;
      }
      return Matrix(vCopy, matrix_size - 1);
    }

    vector<long double> get_full_matrix() const
    {
      vector<long double> full_matrix(get_matrix_size()*get_matrix_size());
      for (auto & e: v)
        full_matrix[ e.row * get_matrix_size() + e.col ] = e.value;
      return full_matrix;
    }

    void print_full_matrix( const vector<long double> & full_matrix )
    {
      const ll n = get_matrix_size();
      cout << setprecision(12) << "{";
      for (long long i = 0; i < n; i++)
      {
        cout << "{";
        for (long long j = 0; j < n; j++)
          cout << "  " << full_matrix[i*n+j] << ",";
        cout << "},";
      }
      cout << "}" << endl;
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
        if (sm.rowsums[index] == 0 || sm.colsums[index] == 0)
          return 0.0;

      // otherwise, check if only one element remaining
      if (sm.get_nonzero_size() == 1 && sm.v[0].row == 0)
        return sm.v[0].value;
      else
      {
        if (sm.get_matrix_size() > 7)
        {
          double sum = 0.0;
          for (const auto & e: sm.v)
          {
            if (e.row > 0)
              break;
            sum += e.value * permanent( sm.get_minor( /*e.row,*/ e.col ) ); // e.row == 0
          }
          return sum;
        }
        else
          return permanent_RNW( sm.get_full_matrix(), sm.get_matrix_size() );
      }
    }
  }
}
