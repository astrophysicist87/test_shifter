#ifndef SYMMETRIC_PAIR_TABLE_H
#define SYMMETRIC_PAIR_TABLE_H

#include <vector>

template<typename T>
class SymmetricPairTable<T>
{
  private:
    long n {0};
    std::vector<T> v;

    inline long UTindexer(const long i, const long j)
    {
      return -1 + j - i*(3 + i - 2*n)/2;
    }

  public:
    SymmetricPairTable(const vector<T> & v_in)
    : n{v_in.size()},
      v{v_in}
    {
    }
    inline T& operator() (const long i, const long j)
    {
      return v[(i<j) ? UTindexer(i, j): UTindexer(j, i)];
    }
    inline const T& operator() (const long i, const long j) const
    {
      return v[(i<j) ? UTindexer(i, j): UTindexer(j, i)];
    }
};


#endif
