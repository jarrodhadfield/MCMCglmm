#include "MCMCglmm.h"

double dcutpoints(const cs *liab, double *yP, int *observed, int start,int finish, double *oldcutpoints, double *newcutpoints, int stcutpoints, int ncutpoints, double sdcp, double sdl)
{
    int i,j,w;
    double llik = 0.0;

    for (j = 2 ; j < (ncutpoints-2); j++){
        llik += log(pnorm(oldcutpoints[stcutpoints+j+1]-oldcutpoints[j], 0.0, sdcp, TRUE,FALSE)-pnorm(newcutpoints[stcutpoints+j-1]-oldcutpoints[j], 0.0, sdcp, TRUE,FALSE));
        llik -= log(pnorm(newcutpoints[stcutpoints+j+1]-newcutpoints[j], 0.0, sdcp, TRUE,FALSE)-pnorm(oldcutpoints[stcutpoints+j-1]-newcutpoints[j], 0.0, sdcp, TRUE,FALSE));
    }

    llik += log(1.0-pnorm(newcutpoints[stcutpoints+ncutpoints-3]-oldcutpoints[stcutpoints+ncutpoints-2], 0.0, sdcp, TRUE,FALSE));
    llik -= log(1.0-pnorm(oldcutpoints[stcutpoints+ncutpoints-3]-newcutpoints[stcutpoints+ncutpoints-2], 0.0, sdcp, TRUE,FALSE));

    for (i = start ; i < finish; i++){
        w = yP[i];
        if(w>1 && observed[i]==1){
          if(w==(ncutpoints-1)){
            llik += log(1.0-pnorm(newcutpoints[stcutpoints+w-1], liab->x[i], sdl, TRUE,FALSE));
            llik -= log(1.0-pnorm(oldcutpoints[stcutpoints+w-1], liab->x[i], sdl, TRUE,FALSE));
          }else{
            llik += log(pnorm(newcutpoints[stcutpoints+w], liab->x[i], sdl, TRUE,FALSE)-pnorm(newcutpoints[stcutpoints+w-1], liab->x[i], sdl, TRUE,FALSE));
            llik -= log(pnorm(oldcutpoints[stcutpoints+w], liab->x[i], sdl, TRUE,FALSE)-pnorm(oldcutpoints[stcutpoints+w-1], liab->x[i], sdl, TRUE,FALSE));
          }
        }
    }
    return llik;
}
