#include <stdlib.h>
#include "Header1new.h"
#define NR_END 1
int  *vectornew(long nl,long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
         int  *v;

	v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror2new("allocation failure in dvector()");
	return v-nl+NR_END;
}
