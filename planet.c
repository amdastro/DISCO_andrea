#include "andrea.h"


double PHI_ORDER = 2.0;

// amd
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
void get_drho_dt(struct planet * pl , struct domain * theDomain , double r , double phi , double rho , double pres , double * drho_dt_sink){
   // From planet data:
   double r_p = pl->r;
   double p_p = pl->phi;
   //double om_p = pl->omega;
   double m_p = pl->M;

   // Get the distance of the cell to the planet:
   double dx = r*cos(phi)-r_p*cos(p_p);
   double dy = r*sin(phi)-r_p*sin(p_p);
   double script_r = sqrt(dx*dx+dy*dy);

   // This is set up only for the secondary planet to accrete
   // So mtotal = q + 1

   double visc_flag = theDomain->theParList.alpha_flag;
   double t_sink_factor  = theDomain->theParList.t_sink_factor;
   double r_sink  = theDomain->theParList.r_sink;

   

   double nu = theDomain->theParList.viscosity;

   if (visc_flag){
      double alpha = theDomain->theParList.viscosity;
      double nu = pow( 1.0/(alpha*pres/rho)*sqrt(pow(script_r,-3)*m_p/(m_p+1.0) ),-1);
   }
 
   // Set the accretion timescale to the local viscous timescale
   double t_visc = 2./3. * script_r*script_r/nu;
   double t_sink = t_sink_factor * t_visc;

   //can read in dt from planet_sink, but then need to read it in to report.c too
   // don't let sink timescale get too short, things could get unstable
   //if (t_sink < 10.* dt){
   //    t_sink = 10.*dt;
   //}

   // Or we can make the timescale some number of orbital times.
   // This is dependent on the planet's orbital frequency, which 
   // will change when it migrates -- think about implications of that! 
   //double t_sink = t_sink_factor / om_p;

   //*drho_dt_sink = 0.0;

   // If r < r_sink, then drho_dt source term is calculated
   //if (script_r < r_sink){
   //   *drho_dt_sink = rho / t_sink;
   //}
   
   *drho_dt_sink = rho / t_sink * exp(- pow( script_r / r_sink , 4.)); 

}

void planet_sink(struct planet * pl , struct domain * theDomain , double * prim , double * cons , double * xp , double * xm , double dVdt ){

   // This is cell data
   double rp = xp[0];
   double rm = xm[0];
   // Radius of the cell
   double r = 0.5*(rp+rm);
   double rho = prim[RHO];
   double pres  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];

   // Phi of cell
   double dphi = get_dp(xp[1],xm[1]);
   double phi = xm[1] + 0.5*dphi;

   double drho_dt_sink;

   // Call get_drho_dt with planet data and cell r, phi
   // rho, pres, viscosity

   get_drho_dt( pl , theDomain , r , phi , rho , pres , &drho_dt_sink );

   // Update all the conservative variables since they have factors of rho
   // better to use primitive variables on the right hand side, since those 
   // aren't being updated throughout a time step   cons[DDD] -= drho_dt_sink*dVdt;
   cons[SRR] -= drho_dt_sink * vr * dVdt;
   cons[TAU] -= .5 * drho_dt_sink * dVdt;
   cons[LLL] -= r * drho_dt_sink * vp * dVdt;
   cons[SZZ] -= drho_dt_sink * vz * dVdt;

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
