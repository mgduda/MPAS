#include <stdlib.h>
#include <sys/time.h>

struct timeval start_time[10];
struct timeval end_time[10];

void start_timer(int id)
{
	gettimeofday(&(start_time[id]), NULL);
}

void stop_timer(int id, int *secs, int *usec)
{
	gettimeofday(&(end_time[id]), NULL);

	*secs = (int)(end_time[id].tv_sec - start_time[id].tv_sec);
	*usec = (int)(end_time[id].tv_usec - start_time[id].tv_usec);

	if (*usec < 0)  {
		*secs   -= 1;
		*usec += 1000000;
	}
}
