#include "cs.h"
/* print a sparse matrix */
int cs_print (const cs *A, int brief)
{
    int p, j, m, n, nzmax, nz, *Ap, *Ai ;
    double *Ax ;
    if (!A) { Rprintf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    Rprintf ("CSparse Version %d.%d.%d, %s.  %s\n", CS_VER, CS_SUBVER,
	CS_SUBSUB, CS_DATE, CS_COPYRIGHT) ;
    if (nz < 0)
    {
	Rprintf ("%d-by-%d, nzmax: %d nnz: %d, 1-norm: %g\n", m, n, nzmax,
		Ap [n], cs_norm (A)) ;
	for (j = 0 ; j < n ; j++)
	{
	    Rprintf ("    col %d : locations %d to %d\n", j, Ap [j], Ap [j+1]-1);
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		Rprintf ("      %d : %g\n", Ai [p], Ax ? Ax [p] : 1) ;
		if (brief && p > 20) { Rprintf ("  ...\n") ; return (1) ; }
	    }
	}
    }
    else
    {
	Rprintf ("triplet: %d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, nz) ;
	for (p = 0 ; p < nz ; p++)
	{
	    Rprintf ("    %d %d : %g\n", Ai [p], Ap [p], Ax ? Ax [p] : 1) ;
	    if (brief && p > 20) { Rprintf ("  ...\n") ; return (1) ; }
	}
    }
    return (1) ;
}
