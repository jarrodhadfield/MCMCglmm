#include "MCMCglmm.h"

double rtnorm(double mu, double sd, double lower, double upper)
{

    double z,pz,u,slower,supper,tr,alpha;
    int sample=1;
 
    if(lower>=upper){
      return((lower+upper)/2);
    }
    if(lower < -1e+32 || upper > 1e+32){
     
      if(lower < -1e+32 && upper > 1e+32){
         z = rnorm(mu, sd);
         return(z);
      }else{
        if(upper > 1e+32){
          tr = (lower-mu)/sd;
        }else{
          tr = (mu-upper)/sd;
        }
        if(tr<0){                          // if sampling >0.5 of a normal density possibly quicker just to sample and reject
          while(sample==1){
            z = rnorm(0.0,1.0);
            if(z>tr){
              sample = 0;
            }
          }
        }else{
          alpha = (tr+sqrt((tr*tr)+4.0))/2.0;
          while(sample==1){
            z = rexp(1.0/alpha)+tr;
            pz = -((alpha-z)*(alpha-z)/2.0);
            u = -rexp(1.0);
            if(u<=pz){
              sample = 0;
            }
          }
        }
      }
    }else{

      slower = (lower-mu)/sd;
      supper = (upper-mu)/sd;

      tr = pnorm(supper, 0.0,1.0, TRUE, FALSE)-pnorm(slower, 0.0,1.0, TRUE, FALSE);

      if(tr>0.5){                   // if sampling >0.5 of a normal density possibly quicker just to sample and reject
        while(sample==1){
          z = rnorm(0.0,1.0);
          if(z>slower && z<supper){
            sample = 0;
          }
        }
      }else{
        while(sample==1){
          z = runif(slower,supper);

          if(slower<=0.0 && 0.0<=supper){
              pz = -z*z/2.0;
          }else{
            if(supper<0.0){
              pz = (supper*supper-z*z)/2.0;
            }else{
              pz = (slower*slower-z*z)/2.0;
            }
          }
          u = -rexp(1.0);
          if(u<pz){
            sample=0;
          }
        }
      }
    }
    if(lower < -1e+32){
      return(mu-z*sd);
    }else{
      return(z*sd+mu);
    }
}



