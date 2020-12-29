/*Rotation*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "input.h"
#include "opac.h"
#include "atmos.h"
#include "constant.h"
#include "include.h"
#include "nrutil.h"

// C HARADA -- update for 2stream //
#include "two_stream.h"


/* --- Global variables ------------------------------------------ */

extern struct Atmos atmos;
extern struct Opac opac;

/* --- Function prototypes --------------------------------------- */

void Locate(int n, double *array, double value, int *ilow);
double Planck(double T, double lambda);
double lint2D(double x1, double x2, double y1, double y2, double z1,
              double z2, double z3, double z4, double x, double y);
double lint3D(double x1, double x2, double y1, double y2, double z1,
              double z2, double f1, double f2, double f3, double f4,
              double f5, double f6, double f7, double f8, double x,
              double y, double z);
void Angles3d(double ds[], double theta[], double dtheta[], double lat);
double Radius(double R_pl, double ds[]);
double lint(double xa, double ya, double xb, double yb, double x);

// C HARADA -- update for 2stream //
double two_stream(int NLAYER, int kmin, double *w0_array, double *g0_array, \
                  const double *temperature_array, const double *tau_array, \
                  double NU, double NU_BIN, double incident_frac, double *dtau_array);

/* ------- begin ---------------- RT_Emit -------------------- */

/* Rest frame or Doppler shifted emission spectra for the disk (not including the limb) with clouds turned on or off */

int RT_Emit_3D(double PHASE)
{
    double ***tau_tr_east, ***tau_tr_west, **theta, **dtheta, ***tau_em, ***dtau_em, ***phi_lon_solid, ***theta_lat_solid, ***temperature_3d,
    ***aero_lw_kappa_1, ***aero_lw_tau_1,
    ***aero_lw_kappa_2, ***aero_lw_tau_2,
    ***aero_lw_kappa_3, ***aero_lw_tau_3,
    ***aero_lw_kappa_4, ***aero_lw_tau_4,
    ***kappa_nu_array, ***pressure_array;
    double **intensity, **bad_interp_array, *flux_st, *flux_pl, *flux_tr, *ds, ***dl, **phi,
    *phi_plus_e, *phi_plus_w, *phi_minus_e, *phi_minus_w, **dphi, *theta_lon, *theta_lat, *phi_lon;
    double R, test, a, b, *lat_rad, kappa_nu_plus_e, kappa_nu_plus_w,
    kappa_nu_minus_e, kappa_nu_minus_w, t_lon_plus_e, t_lon_plus_w,
    t_lon_minus_e, t_lon_minus_w, p_lon_plus_e, p_lon_plus_w,
    p_lon_minus_e, p_lon_minus_w, temperature, pressure, kappa_nu, *lon_rad,
    aero_lw_kappa_interp_1,
    aero_lw_kappa_interp_2,
    aero_lw_kappa_interp_3,
    aero_lw_kappa_interp_4;
    double **I_top, *I_bot, **dkappa_nu;
    int i, j, k, l, m, n, o, c, g, h, ii;
    double dphid, thetad, dthetad;
    FILE *file;
    double solid;
    double average, running_sum, num_points, good_intensity;
    double u_vel, v_vel, w_vel, v_los, delta_lam, omega;
    
    // C HARADA -- update for 2stream //
    int kmin, good_l, good_m, good_val;
    double incident_frac;
    double ***pi0_tot, ***asym_tot;
    
    
    char OUTPUT_FILE[200];
    sprintf(OUTPUT_FILE, "%s%06.2f.dat", OUTPUT_PREFIX, PHASE);
    file = fopen(OUTPUT_FILE, "w");
    
    /*Allocate memory*/
    tau_em = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        tau_em[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            tau_em[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    dtau_em = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        dtau_em[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            dtau_em[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    temperature_3d = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        temperature_3d[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            temperature_3d[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    phi_lon_solid = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        phi_lon_solid[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            phi_lon_solid[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    theta_lat_solid = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        theta_lat_solid[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            theta_lat_solid[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    dl = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        dl[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            dl[l][m] = malloc(NTAU*sizeof(double));
        }
    }


    kappa_nu_array = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        kappa_nu_array[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            kappa_nu_array[l][m] = malloc(NTAU*sizeof(double));
        }
    }

    pressure_array = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        pressure_array[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            pressure_array[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    
    // C HARADA -- update for 2stream //
    
    /* allocate memory for scattering parameters */
    
    pi0_tot = malloc(NLAT*sizeof(double)); // total single scattering albedo
    for(l=0; l<NLAT; l++)
    {
        pi0_tot[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            pi0_tot[l][m] = malloc(NTAU*sizeof(double));
        }
    }

    asym_tot = malloc(NLAT*sizeof(double)); // total asymmetry parameter
    for(l=0; l<NLAT; l++)
    {
        asym_tot[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            asym_tot[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    
    /* allocate memory for aero taus and kappas if clouds on */
    
    if(CLOUDS==1){

        /* MgSiO3 */
        aero_lw_kappa_1 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_kappa_1[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_kappa_1[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_lw_tau_1 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_tau_1[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_tau_1[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        
        /* Fe */
        aero_lw_kappa_2 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_kappa_2[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_kappa_2[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_lw_tau_2 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_tau_2[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_tau_2[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        
        /* Al2O3 */
        aero_lw_kappa_3 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_kappa_3[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_kappa_3[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_lw_tau_3 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_tau_3[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_tau_3[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        
        /* MnS */
        aero_lw_kappa_4 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_kappa_4[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_kappa_4[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_lw_tau_4 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_tau_4[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_tau_4[l][m] = malloc(NTAU*sizeof(double));
            }
        }
            printf("clouds: ON\n");
    }
    
    else{
        printf("clouds: OFF\n");
    }
    
    if(DOPPLER==1){
        printf("doppler: ON\n");
    }
    else{
        printf("doppler: OFF\n");
    }
    
    theta = dmatrix(0, NLAT-1, 0, NLON-1);
    dtheta = dmatrix(0, NLAT-1, 0, NLON-1);
    dphi = dmatrix(0, NLAT-1, 0, NLON-1);
    phi = dmatrix(0, NLAT-1, 0, NLON-1);
    
    intensity = dmatrix(0, NLAT-1, 0, NLON-1);
    bad_interp_array = dmatrix(0, NLAT-1, 0, NLON-1);
    I_top = dmatrix(0, NLAT-1, 0, NLON-1);
    
    flux_pl = dvector(0, NLAMBDA-1);
    ds = dvector(0, NTAU-1);
    
    lat_rad = dvector(0, NLAT-1);
    lon_rad = dvector(0, NLON-1);
    
    /*   Calculate the angular rotational speed omega */
    
    omega = 2.0*PI / (P_ROT*24.0*60.0*60.0);
    
    /*Calculate ds*/
    
    for(j=NTAU-1; j>=0; j--)
    {
        ds[j] = atmos.alt[j-1] - atmos.alt[j];
    }
    ds[0] = ds[1];
    
    /* calculate new aerosol taus and kappas, corrected for wavelength (if clouds on) */
    
    if(CLOUDS==1){
        for(l=0; l<NLAT; l++){
            for(m=0; m<NLON; m++){
                for(j=0; j<NTAU; j++){
                    
                    /* scattering efficiency correction from 5um to 2.3um
                     (PI0, G0, QE calculated from Mie Scattering code of Mischenko, used in Roman 2018) */
                    aero_lw_tau_1[l][m][j] = atmos.aero_sw_tau_1[l][m][j] * (QE_MgSiO3 / 0.01);
                    aero_lw_tau_2[l][m][j] = atmos.aero_sw_tau_2[l][m][j] * (QE_Fe / 0.16);
                    aero_lw_tau_3[l][m][j] = atmos.aero_sw_tau_3[l][m][j] * (QE_Al2O3 / 0.02);
                    aero_lw_tau_4[l][m][j] = atmos.aero_sw_tau_4[l][m][j] * (QE_MnS / 0.02);
                    
                }
            }
        }
        
        for(l=0; l<NLAT; l++){
            for(m=0; m<NLON; m++){
                for(j=0; j<NTAU; j++){
                    
                    aero_lw_kappa_1[l][m][j] = aero_lw_tau_1[l][m][j] / ds[j];
                    aero_lw_kappa_2[l][m][j] = aero_lw_tau_2[l][m][j] / ds[j];
                    aero_lw_kappa_3[l][m][j] = aero_lw_tau_3[l][m][j] / ds[j];
                    aero_lw_kappa_4[l][m][j] = aero_lw_tau_4[l][m][j] / ds[j];
                }
            }
        }
    }
    
    /*Geometry*/
    
    /*Calculating dl longitude and latitude along line-of-sight*/
    for(l=0;l<NLAT;l++)
    {
        
        lat_rad[l] = atmos.lat[l] * PI/180.0; /*theta, latitude*/
        
        for(m=0;m<NLON;m++)
        {
            if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
            {
                
                lon_rad[m] = atmos.lon[m] * PI/180.0; /*phi, longitude*/
                
                b = R_PLANET;
                for(j=NTAU-1; j>=0; j--)
                {
                    b += ds[j];
                    
                    dl[l][m][j] = pow(SQ(b) - SQ(R_PLANET * sin(lat_rad[l])) - SQ(R_PLANET * cos(lat_rad[l]) * cos(-PI/2. + lon_rad[m] + PHASE*PI/180.0)), 0.5);
                    
                    phi_lon_solid[l][m][j] = 450. + acos((R_PLANET * cos(-PI/2. + lon_rad[m] + PHASE*PI/180.0) * cos(lat_rad[l]))/(pow(SQ(b)-SQ(R_PLANET * sin(lat_rad[l])), 0.5))) * 180.0/PI;
                    
                    theta_lat_solid[l][m][j] = asin((R_PLANET * sin(lat_rad[l]))/b) * 180.0/PI;
                }
                for(j=0; j<NTAU; j++)
                {
                    if(j!=NTAU-1)
                    {
                        dl[l][m][j] -= dl[l][m][j+1];
                    }
                    else
                    {
                        dl[l][m][j] = dl[l][m][j] - R_PLANET * cos(lat_rad[l]) * sin(lon_rad[m] + PHASE*PI/180.0 - PI/2.);
                    }
                }
            }
        }
    }
    
    /*Calculating solid angle along NLAT*NLON */
    solid = 0.0;
    for(l=0; l<NLAT; l++)
    {
        for(m=0; m<NLON; m++)
        {
            if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
            {
                if(l<NLAT-1 && l>0)
                {
                    theta[l][m] = theta_lat_solid[l][m][0] * PI/180.0;
                    dtheta[l][m] = 0.5*(theta_lat_solid[l+1][m][0] - theta_lat_solid[l-1][m][0])* PI/180.0;
                }
                theta[0][m] = theta_lat_solid[0][m][0]* PI/180.0;
                theta[NLAT-1][m] = theta_lat_solid[NLAT-1][m][0]* PI/180.0;
                
                dtheta[0][m] = ( 0.5*(theta_lat_solid[1][m][0]+theta_lat_solid[0][m][0]) + 90)* PI/180.0;
                dtheta[NLAT-1][m] = (90 - 0.5*(theta_lat_solid[NLAT-1][m][0]+theta_lat_solid[NLAT-2][m][0]))* PI/180.0;
                
                // printf("%e\n", theta_lat_solid[l][m][0]);
                
                if(m>0 && m<NLON-1)
                {
                    if(atmos.lon[m]>450.0-PHASE && atmos.lon[m]<630.0-PHASE)
                    {
                        dphi[l][m] = 0.5*(phi_lon_solid[l][m+1][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 450.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m+1][0] - phi_lon_solid[l][m][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 630.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                }
                
                else if(m==0)
                {
                    if(atmos.lon[m]>450.0-PHASE && atmos.lon[m]<630.0-PHASE)
                    {
                        dphi[l][m] = 0.5*(phi_lon_solid[l][m+1][0] - phi_lon_solid[l][NLON-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 450.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m+1][0] - phi_lon_solid[l][m][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 630.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m][0] - phi_lon_solid[l][NLON-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                }
                
                else if(m==NLON-1)
                {
                    if(atmos.lon[m]>450.0-PHASE && atmos.lon[m]<630.0-PHASE)
                    {
                        dphi[l][m] = 0.5*(phi_lon_solid[l][0][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 450.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][0][0] - phi_lon_solid[l][m][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 630.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                }
                // Sometimes these go negative and should not lol
                // Love Isaac
                dphi[l][m] = fabs(dphi[l][m]);
                dtheta[l][m] = fabs(dtheta[l][m]);
                
                solid += SQ(cos(theta[l][m]))*cos(phi[l][m]-PI)*dtheta[l][m]*dphi[l][m];
            }
        }
    }
    printf("solid %f\n", solid);
    
    for(i=0; i<NLAMBDA; i++)
    {
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                for(j=0; j<NTAU; j++)
                {
                    tau_em[l][m][j]=0.0;
                }
            }
        }

        /*Optical depth*/
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
                {
                    for(j=0; j<NTAU; j++)
                    {
                        Locate(NLAT, atmos.lat, theta_lat_solid[l][m][j], &o);
                        Locate(NLON, atmos.lon, phi_lon_solid[l][m][j]-PHASE, &c);
                        
                        pressure = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.P_3d[o][c][j], atmos.P_3d[o][c+1][j], atmos.P_3d[o+1][c][j], atmos.P_3d[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                        
                        if(atmos.T_3d[o][c][j] < 100.0 || atmos.T_3d[o][c+1][j] < 100.0)
                        {
                            temperature = 0.0;
                            temperature_3d[l][m][j] = temperature;
                        }
                        else
                        {
                            temperature = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.T_3d[o][c][j], atmos.T_3d[o][c+1][j], atmos.T_3d[o+1][c][j], atmos.T_3d[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            temperature_3d[l][m][j] = temperature;

                            
                        }
                        
                        /* interpolate aero kappas anf calculate total g0 and pi0 values, if clouds on */
                        if(CLOUDS==1){
                            {
                                aero_lw_kappa_interp_1 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_lw_kappa_1[o][c][j], aero_lw_kappa_1[o][c+1][j], aero_lw_kappa_1[o+1][c][j], aero_lw_kappa_1[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            }
                            {
                                aero_lw_kappa_interp_2 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_lw_kappa_2[o][c][j], aero_lw_kappa_2[o][c+1][j], aero_lw_kappa_2[o+1][c][j], aero_lw_kappa_2[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            }
                            {
                                aero_lw_kappa_interp_3 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_lw_kappa_3[o][c][j], aero_lw_kappa_3[o][c+1][j], aero_lw_kappa_3[o+1][c][j], aero_lw_kappa_3[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            }
                            {
                                aero_lw_kappa_interp_4 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_lw_kappa_4[o][c][j], aero_lw_kappa_4[o][c+1][j], aero_lw_kappa_4[o+1][c][j], aero_lw_kappa_4[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            }
                            
                            
                            // C HARADA -- update for 2stream //
                            if (temperature < 300 || dtau_em[l][m][j] < 1e-50)
                            {
                                pi0_tot[l][m][j] = 0.0;
                                asym_tot[l][m][j] = 0.0;
                                
                            }
                            
                            else
                            {
                                double weight_1 = aero_lw_kappa_interp_1 * dl[l][m][j] / dtau_em[l][m][j];
                                double weight_2 = aero_lw_kappa_interp_2 * dl[l][m][j] / dtau_em[l][m][j];
                                double weight_3 = aero_lw_kappa_interp_3 * dl[l][m][j] / dtau_em[l][m][j];
                                double weight_4 = aero_lw_kappa_interp_4 * dl[l][m][j] / dtau_em[l][m][j];
                                                                
                                pi0_tot[l][m][j] =  (weight_1 * PI0_MgSiO3 + \
                                                    weight_2 * PI0_Fe + \
                                                    weight_3 * PI0_Al2O3 + \
                                                    weight_4 * PI0_MnS);
                                
                                asym_tot[l][m][j] = (weight_1 * G0_MgSiO3 + \
                                                    weight_2 * G0_Fe + \
                                                    weight_3 * G0_Al2O3 + \
                                                    weight_4 * G0_MnS);
                                
                                
                                if (pi0_tot[l][m][j] < 0.0 || pi0_tot[l][m][j] > 1.0)
                                {
                                   printf("pi0 total out of range at %i, %i, %i \n", l, m, j);
                                }
                                if (asym_tot[l][m][j] < -1.0 || asym_tot[l][m][j] > 1.0)
                                {
                                   printf("g0 total out of range at %i, %i, %i \n", l, m, j);
                                }
                                
                                
                            }

                        }
                        
                        /* if clouds are turned off, need to set scattering params to zero */
                        else
                        {
                            pi0_tot[l][m][j] = 0.0;
                            asym_tot[l][m][j] = 0.0;
                        }
                        
                    
                        
                        Locate(NTEMP, opac.T, temperature, &g);
                        Locate(NPRESSURE, opac.P, pressure, &h);
                        
                        /* Add doppler shift to signal, if turned on */
                        
                        if(DOPPLER==1){
                            u_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ew[o][c][j], atmos.vel_ew[o][c+1][j], atmos.vel_ew[o+1][c][j], atmos.vel_ew[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            
                            v_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ns[o][c][j], atmos.vel_ns[o][c+1][j], atmos.vel_ns[o+1][c][j], atmos.vel_ns[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            
                            w_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ve[o][c][j], atmos.vel_ve[o][c+1][j], atmos.vel_ve[o+1][c][j], atmos.vel_ve[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            
                            v_los = u_vel*sin(phi_lon_solid[l][m][j]*PI/180.0) + v_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*sin(theta_lat_solid[l][m][j]*PI/180.0) - w_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + omega*(R_PLANET + atmos.alt[j])*sin(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + R_VEL*cos((90.0-PHASE)*PI/180.0); /*Everything*/
                            
                            /* v_los = u_vel*sin(phi_lon_solid[l][m][j]*PI/180.0) + v_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*sin(theta_lat_solid[l][m][j]*PI/180.0); */ /*ewns*/
                            
                            /* v_los =  - w_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0); */ /*vertical*/
                            
                            /* v_los = u_vel*sin(phi_lon_solid[l][m][j]*PI/180.0) + v_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*sin(theta_lat_solid[l][m][j]*PI/180.0) - w_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0); */ /*wind*/
                            
                            /* v_los = omega*(R_PLANET + atmos.alt[j])*sin(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + R_VEL*cos((90.0-PHASE)*PI/180.0); */ /*rotation*/
                            
                            
                            
                            delta_lam = atmos.lambda[i]*v_los/CLIGHT;
                            Locate(NLAMBDA, atmos.lambda, atmos.lambda[i]+delta_lam, &ii);
                            
                            if(temperature < 100.0 || atmos.lambda[i]+delta_lam >= atmos.lambda[NLAMBDA-1] || atmos.lambda[i]+delta_lam <= atmos.lambda[0])
                            {
                                kappa_nu = 0.0;
                            }
                            else
                            {
                                kappa_nu = lint3D(opac.T[g], opac.T[g+1], opac.P[h],
                                                  opac.P[h+1],
                                                  atmos.lambda[ii], atmos.lambda[ii+1],
                                                  opac.kappa[ii][h][g],
                                                  opac.kappa[ii][h][g+1],
                                                  opac.kappa[ii][h+1][g],
                                                  opac.kappa[ii][h+1][g+1],
                                                  opac.kappa[ii+1][h][g],
                                                  opac.kappa[ii+1][h][g+1],
                                                  opac.kappa[ii+1][h+1][g],
                                                  opac.kappa[ii+1][h+1][g+1],
                                                  temperature, pressure, atmos.lambda[i]+delta_lam);
                                
                                /* add aerosol opacities if clouds turned on */
                                if (CLOUDS==1) {
                                    kappa_nu += aero_lw_kappa_interp_1 + aero_lw_kappa_interp_2 + aero_lw_kappa_interp_3 + aero_lw_kappa_interp_4;
                                }
                                
                            }
                        }


                        /* Wind Only */
                        else if(DOPPLER==2){
                            u_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ew[o][c][j], atmos.vel_ew[o][c+1][j], atmos.vel_ew[o+1][c][j], atmos.vel_ew[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ns[o][c][j], atmos.vel_ns[o][c+1][j], atmos.vel_ns[o+1][c][j], atmos.vel_ns[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            w_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ve[o][c][j], atmos.vel_ve[o][c+1][j], atmos.vel_ve[o+1][c][j], atmos.vel_ve[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_los = u_vel*sin(phi_lon_solid[l][m][j]*PI/180.0) + v_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*sin(theta_lat_solid[l][m][j]*PI/180.0) - w_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + R_VEL*cos((90.0-PHASE)*PI/180.0);

                            delta_lam = atmos.lambda[i]*v_los/CLIGHT;
                            Locate(NLAMBDA, atmos.lambda, atmos.lambda[i]+delta_lam, &ii);

                            if(temperature < 250.0 || atmos.lambda[i]+delta_lam >= atmos.lambda[NLAMBDA-1] || atmos.lambda[i]+delta_lam <= atmos.lambda[0])
                            {
                                kappa_nu = 0.0;
                            }
                            else
                            {
                                kappa_nu = lint3D(opac.T[g], opac.T[g+1], opac.P[h],
                                                  opac.P[h+1], atmos.lambda[ii], atmos.lambda[ii+1],
                                                  opac.kappa[ii][h][g],
                                                  opac.kappa[ii][h][g+1],
                                                  opac.kappa[ii][h+1][g],
                                                  opac.kappa[ii][h+1][g+1],
                                                  opac.kappa[ii+1][h][g],
                                                  opac.kappa[ii+1][h][g+1],
                                                  opac.kappa[ii+1][h+1][g],
                                                  opac.kappa[ii+1][h+1][g+1],
                                                  temperature, pressure, atmos.lambda[i]+delta_lam);
                            }
                        }

                        /* Rotation Only */
                        else if(DOPPLER==3){
                            u_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ew[o][c][j], atmos.vel_ew[o][c+1][j], atmos.vel_ew[o+1][c][j], atmos.vel_ew[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ns[o][c][j], atmos.vel_ns[o][c+1][j], atmos.vel_ns[o+1][c][j], atmos.vel_ns[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            w_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ve[o][c][j], atmos.vel_ve[o][c+1][j], atmos.vel_ve[o+1][c][j], atmos.vel_ve[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_los = (cos(INPUT_INCLINATION)*(omega*(R_PLANET + atmos.alt[j])*sin(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + R_VEL*cos((90.0-PHASE)*PI/180.0)));

                            delta_lam = atmos.lambda[i]*v_los/CLIGHT;
                            Locate(NLAMBDA, atmos.lambda, atmos.lambda[i]+delta_lam, &ii);
                            
                            if(temperature < 250.0 || atmos.lambda[i]+delta_lam >= atmos.lambda[NLAMBDA-1] || atmos.lambda[i]+delta_lam <= atmos.lambda[0])
                            {
                                kappa_nu = 0.0;
                            }
                            else
                            {
                                kappa_nu = lint3D(opac.T[g], opac.T[g+1], opac.P[h],
                                                  opac.P[h+1], atmos.lambda[ii], atmos.lambda[ii+1],
                                                  opac.kappa[ii][h][g],
                                                  opac.kappa[ii][h][g+1],
                                                  opac.kappa[ii][h+1][g],
                                                  opac.kappa[ii][h+1][g+1],
                                                  opac.kappa[ii+1][h][g],
                                                  opac.kappa[ii+1][h][g+1],
                                                  opac.kappa[ii+1][h+1][g],
                                                  opac.kappa[ii+1][h+1][g+1],
                                                  temperature, pressure, atmos.lambda[i]+delta_lam);
                            }
                        }


                        
                        else{
                            if(temperature < 100.0)
                            {
                                kappa_nu = 0.0;
                            }
                            else
                            {
                                kappa_nu = lint2D(opac.T[g], opac.T[g+1], opac.P[h],
                                                  opac.P[h+1], opac.kappa[i][h][g],
                                                  opac.kappa[i][h][g+1],
                                                  opac.kappa[i][h+1][g],
                                                  opac.kappa[i][h+1][g+1],
                                                  temperature, pressure);
                                
                                /* add aerosol opacities if clouds turned on */
                                if (CLOUDS==1) {
                                    kappa_nu += aero_lw_kappa_interp_1 + aero_lw_kappa_interp_2 + aero_lw_kappa_interp_3 + aero_lw_kappa_interp_4;
                                }
                            }
                        }
                        
                        if(kappa_nu < 0.0)
                        {
                            printf("WARNING: kappa_nu < 0 at %i, %i, %i \n", l, m, j);
                            
                            if (j > 0)
                            {
                                kappa_nu = dtau_em[l][m][j-1] / dl[l][m][j-1];
                            }
                            else
                            {
                                kappa_nu = 0.0;
                            }
                        }

                        dtau_em[l][m][j] = kappa_nu * dl[l][m][j];
                        kappa_nu_array[l][m][j] = kappa_nu;
                        pressure_array[l][m][j] = pressure;
                    }         
                }
            }
        }
        
        
        
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
                {
                    for(j=0; j<NTAU; j++)
                    {
                        if(j==0)
                        {
                            tau_em[l][m][j]=dtau_em[l][m][j];
                        }
                        else if (j!=0)
                        {
                            tau_em[l][m][j] = tau_em[l][m][j-1] + dtau_em[l][m][j];
                        }
                    }
                }
            }
        }
        

        
        //Calculate the intensity of emergent rays at each latitude and longitude
        
        // C HARADA -- update for 2stream //
        running_sum = 0.0;
        num_points = 0.0;
        average = 0.0;
        good_intensity = 0.0;
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
                {                    
                    // Calculate the incident fraction of starlight for 2stream //
                    incident_frac = cos(lon_rad[m]) * cos(lat_rad[l]);
                    if (incident_frac < 0)
                    {
                        incident_frac = 0.0;
                    }
                    
                    // Find min vert level for 2stream //
                    kmin = 0;
                    for (j = 0; j<NTAU; j++)
                    {
                        if (dtau_em[l][m][j] < 1e-90 || tau_em[l][m][j] < 1e-90 || temperature_3d[l][m][j] < 500)
                        {
                          kmin = j+1;
                        }
                    }

                    // SHHHHH
                    kmin = kmin + 5;
                    
                    if (kmin >= NTAU)
                    {
                        intensity[l][m] = 0;
                    }
                    
                    else
                    {
                        intensity[l][m] = two_stream(NTAU, kmin, pi0_tot[l][m], asym_tot[l][m], temperature_3d[l][m], tau_em[l][m], \
                                                     CLIGHT / atmos.lambda[i], CLIGHT / atmos.lambda[i] - CLIGHT / atmos.lambda[i+1], 
                                                     incident_frac, dtau_em[l][m]);
                    }
                }
            }
        }
        

        /*
        // ~~~ THIS IS THE OLD RT ROUTINE ~~~ //
        
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
                {
                    intensity[l][m] = Planck(atmos.T_3d[l][m][NTAU-1], atmos.lambda[i]) * exp(-tau_em[l][m][NTAU-1]);
                    //intensity[l][m] = Planck(1000, atmos.lambda[i]) * exp(-tau_em[l][m][NTAU-1]);


                    int j_min = 0;
                    for (j=0; j<NTAU; j++)
                    {
                        if (temperature_3d[l][m][j] < 100)
                        {
                            j_min += 1;
                        }
                    }
                    I_top[l][m] = Planck(temperature_3d[l][m][j_min], atmos.lambda[i]);
                    //I_top[l][m] = Planck(1000, atmos.lambda[i]);
                }
            }
        }



        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
                {
                    for(j=0; j<NTAU; j++)
                    {
                        intensity[l][m] += Planck(temperature_3d[l][m][j], atmos.lambda[i]) * exp(-tau_em[l][m][j]) * dtau_em[l][m][j];
                        //intensity[l][m] += Planck(1000, atmos.lambda[i]) * exp(-tau_em[l][m][j]) * dtau_em[l][m][j];
                    }
                }
            }
        }
        */
 
        
        /*Calculate the total flux received by us*/
        //FILE *fptr = fopen("/home/imalsky/Desktop/test.txt", "w"); 
        flux_pl[i] = 0.0;
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
                {
                    //fprintf(fptr, "%d, %d, %.8e\n", l, m, intensity[l][m]);   
                    flux_pl[i] += intensity[l][m] * SQ(cos(theta[l][m])) * cos(phi[l][m]-PI) * dtheta[l][m] * dphi[l][m];
                }
            }
        }
        //fclose(fptr); 
        
        if(i % 100 == 0)
        {
            printf("%d out of %d lines (phase: %06.2f)\n", i, NLAMBDA, PHASE);
        }
        
        fprintf(file, "%10.8le\t%le\n", atmos.lambda[i], flux_pl[i] * PI/solid);
        
    }
    
    fclose(file);
    
    return 0;
}
