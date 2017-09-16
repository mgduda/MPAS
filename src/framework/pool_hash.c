#include <stdlib.h>

#if (defined(UNDERSCORE) || defined(FORTRANUNDERSCORE))
#define pool_hash pool_hash_
#else
#if (defined(DOUBLEUNDERSCORE) || defined(FORTRANDOUBLEUNDERSCORE))
#define pool_hash pool_hash__
#endif
#endif

void pool_hash(int* hash, char* key, int* len)
{
	int i;
	unsigned int whash;

	whash = 0;

	for (i=0; i<(*len); i++) {
		whash += (unsigned int)key[i];
	} 

	*hash = (int)(whash & 0x7fffffff);
}
