#include <math.h>

#ifdef SINGLE_PRECISION
#define MPAS_C_RKIND float
#else
#define MPAS_C_RKIND double
#endif

#ifdef UNDERSCORE
#define nan_check nan_check_
#else
#ifdef DOUBLEUNDERSCORE
#define nan_check nan_check__
#endif
#endif

int nan_check_( MPAS_C_RKIND * x, int * n )
{
   int i;
   MPAS_C_RKIND * d;

   d = x;
   for (i=0; i<(*n); i++) {
      if (isnan(*d) || isinf(*d)) return (i+1);
      d++;
   }
   
   return 0;
}
