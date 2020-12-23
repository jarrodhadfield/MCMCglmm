#include "MCMCglmmcc.h"

extern "C"{  

void pkkR(
      int *k,
      double *prob,
      double *size,
      double *p
){         


  p[0] = pkk(prob, size[0], k[0]);

}
}
