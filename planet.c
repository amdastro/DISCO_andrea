#include "andrea.h"

//amd
static double nu = 0.0;
static double t_sink_factor = 0.0;
static double r_sink = 0.0;


double PHI_ORDER = 2.0;

double get_dp( double , double );

double phigrav( double M , double r , double eps ){
   double n = PHI_ORDER;
   return( M/pow( pow(r,n) + pow(eps,n) , 1./n ) ) ;
}

double fgrav( double M , double r , double eps ){
   double n = PHI_ORDER;
   return( M*pow(r,n-1.)/pow( pow(r,n) + pow(eps,n) ,1.+1./n) );
}

void adjust_gas( struct planet * pl , double * x , double * prim , double gam ){

   double r   = x[0];
   double phi = x[1];

   double rp = pl->r;
   double pp = pl->phi;
   double cosp = cos(phi);
   double sinp = sin(phi);
   double dx = r*cosp-rp*cos(pp);
   double dy = r*sinp-rp*sin(pp);
   double script_r = sqrt(dx*dx+dy*dy);

   double pot = phigrav( pl->M , script_r , pl->eps );

   double c2 = gam*prim[PPP]/prim[RHO];
   double factor = 1. + (gam-1.)*pot/c2;

   prim[RHO] *= factor;
   prim[PPP] *= pow( factor , gam );

}

void planetaryForce( struct planet * pl , double r , double phi , double z , double * fr , double * fp , double * fz , int mode ){

   z = 0.0;

   double rp = pl->r;
   double pp = pl->phi;
   double cosp = cos(phi);
   double sinp = sin(phi);
   double dx = r*cosp-rp*cos(pp);
   double dy = r*sinp-rp*sin(pp);
   double script_r = sqrt(dx*dx+dy*dy+z*z);
   double script_r_perp = sqrt(dx*dx+dy*dy);

   double f1 = -fgrav( pl->M , script_r , pl->eps );

   double cosa = dx/script_r;
   double sina = dy/script_r;

   double cosap = cosa*cosp+sina*sinp;
   double sinap = sina*cosp-cosa*sinp;

   if( mode==1 ){
      cosap = cosa*cos(pp)+sina*sin(pp);
      sinap = sina*cos(pp)-cosa*sin(pp);
   }
/*
   double rH = rp*pow( pl->M/3.,1./3.);
   double pd = 0.8; 
   double fd = 1./(1.+exp(-( script_r/rH-pd)/(pd/10.)));
*/

   double sint = script_r_perp/script_r;
   double cost = z/script_r;

   *fr = cosap*f1*sint; //*fd;
   *fp = sinap*f1*sint; //*fd;
   *fz = f1*cost;

}

void planet_src( struct planet * pl , double * prim , double * cons , double * xp , double * xm , double dVdt ){

   double rp = xp[0];
   double rm = xm[0];
   double rho = prim[RHO];
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   double omega = prim[UPP];
   
   double r = 0.5*(rp+rm);
   double vp  = r*omega;
   double dphi = get_dp(xp[1],xm[1]);
   double phi = xm[1] + 0.5*dphi;
   double z = .5*(xp[2]+xm[2]);

   double Fr,Fp,Fz;
   planetaryForce( pl , r , phi , z , &Fr , &Fp , &Fz , 0 );

   cons[SRR] += rho*Fr*dVdt;
   cons[SZZ] += rho*Fz*dVdt;
   cons[LLL] += rho*Fp*r*dVdt;
   cons[TAU] += rho*( Fr*vr + Fz*vz + Fp*vp )*dVdt;

}

//amd
void planet_sink(struct planet * pl , struct domain * theDomain , double * prim , double * cons , double * xp , double * xm , double dVdt ){
   // This is assuming we have nu viscosity
   // In the future can read in the alpha flag to make this more general.
   // Then if alpha is specified, can calculate the corresponding nu.
   nu   = theDomain->theParList.viscosity;
   t_sink_factor  = theDomain->theParList.t_sink_factor;
   r_sink  = theDomain->theParList.r_sink;

   // Get planet data
   double r_p = pl->r;
   double p_p = pl->phi;

   // This is cell data? What is rp and rm?
   // Not sure if we need all these
   double rp = xp[0];
   double rm = xm[0];
   double rho = prim[RHO];

   // some measure of the radius of the cell?
   // Is this distance of cell to planet?
   double r = 0.5*(rp+rm);

   // t_sink = some factor times the viscous time
   double t_visc = 2./3. * r*r/nu;
   double t_sink = t_sink_factor * t_visc;

   double drho_dt_sink = 0.0;

   // If r < r_sink, then drho_dt source term is calculated
   if (r < r_sink){ 
      drho_dt_sink = rho / t_sink;
   }

   // update conservative variables
   cons[RHO=] -= drho_dt_sink*dVdt;
   // Yike also updated the following, since these all contain a factor of rho
   //cons[DDD] -= drho_dt_sink *dt/rho * cons[DDD];
   //cons[TAU] -= drho_dt_sink *dt/rho * cons[TAU];
   //cons[SRR] -= drho_dt_sink *dt/rho * cons[SRR];
   //cons[LLL] -= drho_dt_sink *dt/rho * cons[LLL];
   //cons[SZZ] -= drho_dt_sink *dt/rho * cons[SZZ];

}


void planet_RK_copy( struct planet * pl ){
   pl->RK_r     = pl->r;
   pl->RK_phi   = pl->phi;
   pl->RK_M     = pl->M;
   pl->RK_omega = pl->omega;
   pl->RK_vr    = pl->vr;
}

void planet_RK_adjust( struct planet * pl , double RK ){
   pl->r     = (1.-RK)*pl->r     + RK*pl->RK_r;
   pl->phi   = (1.-RK)*pl->phi   + RK*pl->RK_phi;
   pl->M     = (1.-RK)*pl->M     + RK*pl->RK_M;
   pl->omega = (1.-RK)*pl->omega + RK*pl->RK_omega;
   pl->vr    = (1.-RK)*pl->vr    + RK*pl->RK_vr;
}
