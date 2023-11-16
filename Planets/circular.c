
#include "../andrea.h"

static double q = 1.0;
static double mach = 0.0;
static double eps_frac = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2; 
   q = theDomain->theParList.Mass_Ratio;
   mach = theDomain->theParList.Disk_Mach;
   eps_frac = theDomain->theParList.eps_frac;

}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   // r is total separation
   double r = 1.0;
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
   thePlanets[1].eps   = eps_frac/mach;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   thePlanets[0].phi += thePlanets[0].omega*dt;
   thePlanets[1].phi += thePlanets[1].omega*dt;
}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

