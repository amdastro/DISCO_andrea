#include "andrea.h"

void planetaryForce( struct planet * , double , double , double , double * , double * , double * , int );
void get_drho_dt(struct planet * , double , double , double , double * );


double get_dV( double * , double * );

void report( struct domain * theDomain ){

   double t = theDomain->t;
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Ng = theDomain->Ng;
   int rank = theDomain->rank;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   MPI_Comm grid_comm = theDomain->theComm;

   double gamma_law = theDomain->theParList.Adiabatic_Index;

   struct planet * thePlanets = theDomain->thePlanets;
   int Npl = theDomain->Npl;

   double r_p = 0.0;
   double p_p = 0.0;
   double rprimary = thePlanets[0].r;
   double pprimary = thePlanets[0].phi;
   if( Npl > 1 ){
     r_p = thePlanets[1].r;
     p_p = thePlanets[1].phi;   
   }

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int sink_flag = theDomain->theParList.sink_flag;


   int jmin = Ng;
   int jmax = Nr-Ng;
   if( dim_rank[0]==0             ) jmin = 0;
   if( dim_rank[0]==dim_size[0]-1 ) jmax = Nr;
   int kmin = Ng;
   int kmax = Nz-Ng;
   if( dim_rank[1]==0             ) kmin = 0;
   if( dim_rank[1]==dim_size[1]-1 ) kmax = Nz;

   int j,k,i;
   //double Br2     = 0.0;
   //double B2      = 0.0;
   double Power  = 0.0;
   double Torque_p = 0.0;
   double Torque_primary = 0.0;
   double Fr=0.0;
   double PsiR = 0.0;
   double PsiI = 0.0;
   double Vol = 0.0;
   double rho_min = HUGE_VAL;
   double rhoavg_min = HUGE_VAL;
   double rho_peak = HUGE_VAL;
   double rho_rp_avg = HUGE_VAL;
   double Mass = 0.0;
   double Mdot = 0.0;

   // define mdot and accretion torque on secondary
   double MdotP = 0.0;
   double JdotP = 0.0;

   //double BrBp = 0.0;
   double PdV  = 0.0;

   double S_R = 0.0;
   double S_0 = 0.0;

   double drho_dt_sink = 0.0;

   // amd: Torque inside Hill radius
   double Torq_hill = 0.0;
   //double T_cut[10];
   //double P_cut[10];
   //for( j=0 ; j<10 ; ++j ){ T_cut[j]=0.;  P_cut[j]=0.; }

   for( j=jmin ; j<jmax ; ++j ){
      double r = .5*(r_jph[j]+r_jph[j-1]);
      double rho0 = 1.0;//pow( r , -1.5 );
      double rho_avg = 0.0;
      double Vol_avg = 0.0;
      double rho_peak = 0.0;
      double rho_rp_avg = 0.0;
      for( k=kmin ; k<kmax ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phi = c->piph - .5*c->dphi;
            double Pp  = c->prim[PPP];
            double rho = c->prim[RHO];
            double pres= c->prim[PPP];
            //amd:
            double omega = c->prim[UPP];

            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV( xp , xm );
	    //bug fix for 2D quantities:
	    if( Nz == 1 ) dV /= z_kph[0]-z_kph[-1];

            PsiR += rho*dV*cos(phi);
            PsiI += rho*dV*sin(phi);

            //if( NUM_Q > BRR ){
            //   double Br = c->prim[BRR];
            //   double Bp = c->prim[BPP];
            //   double Bz = c->prim[BZZ];
            //   Br2 += .5*Br*Br*dV;
            //   BrBp += Br*Bp*dV;
            //   B2  += .5*(Br*Br+Bp*Bp+Bz*Bz)*dV;
            //}
            PdV += Pp*dV;
            Mdot += 2.*M_PI*r*rho*dV*c->prim[URR];
            Vol += dV;

            //find rho_peak
            if( rho_peak > rho/rho0 ) rho_peak = rho/rho0

            if( rho_min > rho/rho0 ) rho_min = rho/rho0;
            rho_avg += rho*dV;
            Vol_avg += dV;
            Mass += rho*dV;

            S_R += pow(rho,4.)*r*dV;
            S_0 += pow(rho,4.)*dV;

            // Calculate quantities on the secondary
            if( Npl > 1 ){
               double fr,fp,fz,fp2;
               // rprimary = primary position, r_p = secondary/planet position
               // careful with notation consistency between primary and secondary! yikes!
               double r_p = thePlanets[1].r;
	            double rprimary = thePlanets[0].r;
	            double pprimary = thePlanets[0].phi;
               double om = thePlanets[1].omega;
               double vr = thePlanets[1].vr;
               double eps = thePlanets[1].eps;
               //double mp = thePlanets[1].M;
               planetaryForce( thePlanets+0 , r , phi , 0.0 , &fr , &fp , &fz , 1 );
               planetaryForce( thePlanets+1 , r , phi , 0.0 , &fr , &fp2  , &fz , 1 );

               // amd:
               // Call planetary sink here to return drho_dt
               // If there were multiple BHs, this would need to
               // loop over each planet
               if( sink_flag ){ 
                  get_drho_dt( thePlanets+1 , r , phi , rho , &drho_dt_sink );
               }
               // Accretion rate
               MdotP += drho_dt_sink*dV;

               // Accretion torque
               double vp = r*omega;
               double dphi2 = phi - thePlanets[1].phi;
               double vp2 = vp * cos(dphi2) + vr*sin(dphi2);
               double dj = r_p * (vp2 - r_p*om);

               JdotP += drho_dt_sink * dj * dV;
               //

               //double cosp = cos(phi);
               //double sinp = sin(phi);
               //double dx = r*cosp-r_p*cos(pp);
               //double dy = r*sinp-r_p*sin(pp);
               //double script_r = sqrt(dx*dx+dy*dy);
               //double rH = pow( thePlanets[1].M/3. , 1./3. );

               Torque_p -= (rho-1.0)*r_p*fp2*dV;
               Power  -= (rho-1.0)*( r_p*om*fp + vr*fr )*dV;
               Torque_primary -= (rho-1.0)*rprimary*fp*dV;
               Fr -= (rho-1.0)*fr*dV;

               double pp = thePlanets[1].phi;
               double rhill = pow(thePlanets[1].M/3., 1./3. ) * r_p;
               
               double scriptr2 = r_p*r_p + r*r - 2.*r*r_p*cos(phi-pp);
               if( scriptr2 < rhill*rhill ){
                  Torq_hill -= (rho-1.0)*r_p*fp*dV;
               }
               // sum up density within (smoothing) annulus of secondary
               if( r > (r_p - eps) && r < (r_p + eps) ){
                  rho_rp_avg+=rho*dV;
               }

/*
               int n_cut;
               for( n_cut=0 ; n_cut<10 ; ++n_cut ){
                  double r_cut = 0.3*((double)(n_cut+1.)/10.);
                  double scriptr2 = r_p*r_p + r*r - 2.*r*r_p*cos(phi-pp);
                  if( scriptr2 > r_cut*r_cut ){
                     T_cut[n_cut] -= (rho-1.0)*r_p*fp*dV;
                     P_cut[n_cut] -= (rho-1.0)*( r_p*om*fp + vr*fr )*dV;
                  }
               }
*/

            }
         }
      }
      rho_avg /= Vol_avg;
      rho_rp_avg /= Vel_avg;
      if( rhoavg_min > rho_avg/rho0 ) rhoavg_min = rho_avg/rho0;
   }

   //MPI_Allreduce( MPI_IN_PLACE , &Br2     , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   //MPI_Allreduce( MPI_IN_PLACE , &B2      , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   //MPI_Allreduce( MPI_IN_PLACE , &BrBp    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &PdV     , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Vol     , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Torque_primary  , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Torque_p , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &MdotP   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &JdotP   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Power   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Fr      , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &PsiR    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &PsiI    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &rho_min    , 1 , MPI_DOUBLE , MPI_MIN , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &rhoavg_min , 1 , MPI_DOUBLE , MPI_MIN , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Mass    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &S_R     , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &S_0     , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Mdot    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Torq_hill   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &rho_peak   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &rho_rp_avg   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );

   Mdot /= Vol;
   S_R /= S_0;

   //double aM = BrBp/PdV;
   //double bM = PdV/B2;

   if( rank==0 ){
      FILE * rFile = fopen("report.dat","a");
      fprintf(rFile,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                t,Torque_p,Torq_hill,Torque_primary,r_p,p_p,MdotP,JdotP,Fr,rho_peak,rho_rp_avg,rprimary,pprimary);
      //fprintf(rFile,"%e %e %e ",t,Torque,Power);
      //for( j=0 ; j<10 ; ++j ) fprintf(rFile,"%e %e ",T_cut[j],P_cut[j]);
      //fprintf(rFile,"\n");
      fclose(rFile);
   }

}
