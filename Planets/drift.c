#include "../andrea.h"

static double rate = 0.0;
static double dexp = 0.0;
static double q = 0.0;
static double t_max = 0.0;
static double mach = 0.0;
static double eps_frac = 0.0;
static double r_init = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2;

   rate  = theDomain->theParList.Drift_Rate;
   dexp  = theDomain->theParList.Drift_Exp;
   q     = theDomain->theParList.Mass_Ratio;
   t_max = theDomain->theParList.t_max;
   mach  = theDomain->theParList.Disk_Mach;
   eps_frac = theDomain->theParList.eps_frac;
   r_init = theDomain->theParList.r_init;

}

int planet_motion_analytic( void ){
   return(1);
}

double drift_pos( double R , double t ){
   // planet initial position is set by drift rate and sim time
   return( r_init * pow( 1. + R*(t - t_max)/dexp , dexp) );
}

void initializePlanets( struct planet * thePlanets ){

   // r is total separation
   double r = drift_pos(rate , 0.0);
   double mu = q/(1.0 + q);
   //units: G*M_total = 1
   double om = pow(r,-1.5); 

   thePlanets[0].M     = (1.0 - mu); 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = om; 
   thePlanets[0].r     = r*mu; 
   thePlanets[0].phi   = 0.0; 
   thePlanets[0].eps   = 0.0;

   thePlanets[1].M     = mu; 
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = om; 
   thePlanets[1].r     = r*(1.0 - mu); 
   thePlanets[1].phi   = 0.0; 
   //amd: smoothing length is some fraction of scale height
   thePlanets[1].eps   = eps_frac/mach;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){

   double mu = q/(1.0 + q);
   double r = drift_pos(rate , t+dt);
   
   thePlanets[0].r     = r*mu;
   thePlanets[0].phi += thePlanets[0].omega*dt;
   thePlanets[1].r     = r*(1.0 - mu);
   thePlanets[1].phi  += thePlanets[1].omega*dt;
   // update smoothing length as planet drifts within a Mach profile disk
   thePlanets[1].eps   = eps_frac * r/mach;

}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

