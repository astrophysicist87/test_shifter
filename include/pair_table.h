#ifndef PAIR_TABLE_H
#define PAIR_TABLE_H

#include <vector>

template<typename T>
class PairTable<T>
{
  private:
    std::vector<long> v;

    inline long UTindexer(long i, long j, long n)
    {
      return -1 + j - i*(3 + i - 2*n)/2;
    }

  public:
    
};


#endif
