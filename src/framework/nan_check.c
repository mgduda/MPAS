#include <math.h>

int nan_check_(double * x, int * n)
{
   int i;
   double * d;

   d = x;
   for (i=0; i<(*n); i++) {
      if (isnan(*d) || isinf(*d)) return (i+1);
      d++;
   }
   
   return 0;
}
