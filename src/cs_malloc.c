#include "cs.h"
#ifdef MATLAB_MEX_FILE
#define malloc mxMalloc
#define free mxFree
#define realloc mxRealloc
#define calloc mxCalloc
#endif

/* wrapper for malloc */
void *cs_malloc (int n, size_t size)
{
    return (CS_OVERFLOW (n,size) ? NULL : malloc (CS_MAX (n,1) * size));
}

/* wrapper for calloc */
void *cs_calloc (int n, size_t size)
{
    return (CS_OVERFLOW (n,size) ? NULL : calloc (CS_MAX (n,1), size));
}

/* wrapper for free */
void *cs_free (void *p)
{
    if (p) free (p) ;	    /* free p if it is not already NULL */
    return (NULL) ;	    /* return NULL to simplify the use of cs_free */
}

/* wrapper for realloc */
void *cs_realloc (void *p, int n, size_t size, int *ok)
{
    void *p2 ;
    *ok = !CS_OVERFLOW (n,size) ;	    /* guard against int overflow */
    if (!(*ok)) return (p) ;		    /* p unchanged if n too large */
    p2 = realloc (p, CS_MAX (n,1) * size) ; /* realloc the block */
    *ok = (p2 != NULL) ;
    return ((*ok) ? p2 : p) ;		    /* return original p if failure */
}
