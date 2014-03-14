#include <stdlib.h>
#include <sys/time.h>

struct timeval start_time[10];
struct timeval end_time[10];

void start_timer_(int * n)
{
   gettimeofday(&(start_time[*n]), NULL);
}

void stop_timer_(int * n, int * secs, int * usec)
{
   gettimeofday(&(end_time[*n]), NULL);

   *secs = (int)(end_time[*n].tv_sec - start_time[*n].tv_sec);
   *usec = (int)(end_time[*n].tv_usec - start_time[*n].tv_usec);

   if (*usec < 0)  {
      *secs   -= 1;
      *usec += 1000000;
   }
}
