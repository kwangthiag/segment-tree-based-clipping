#include "util/util.h"

/* Function to check if x is power of 2*/
bool Util :: isPowerOfTwo(int n)
{
  if (n == 0)
    return 0;
  while (n != 1)
  {
    if (n%2 != 0)
      return 0;
    n = n/2;
  }
  return 1;
}

int Util :: getNPowOf2(int num)
{
    bool isNPowerOf2 = false;
    while(!isNPowerOf2)
    {
        if(isPowerOfTwo(num))
          isNPowerOf2 = true;
        else  
          num++;
    }
    
    return num;
}