/*------------ file ------- readTP.c -----------------------------

Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified: October 20, 2009

------------------------------------------------------------------ */

/* Reads in the temperature - pressure profile from the file 
   T_P_FILE defined in input.h
------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>

#include "input.h"
#include "atmos.h"
#include "nrutil.h"

/* --- Global variables ------------------------------------------ */

extern struct Atmos atmos;

/* --- Function prototypes --------------------------------------- */

/* ------- begin --------------------- ReadTP.c ------------------ */

void ReadTP_3D()
{
  double dum;
  double num;
  int i, j, k;
  FILE *file;

  /* Allocate Memory */

  atmos.lat = dvector(0, NLAT-1);
  atmos.lon = dvector(0, NLON-1);
  atmos.alt = dvector(0, NTAU-1);

  atmos.P_3d = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.P_3d[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.P_3d[i][j] = malloc(NTAU*sizeof(double));
    }
  }
  atmos.T_3d = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.T_3d[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.T_3d[i][j] = malloc(NTAU*sizeof(double));
    }
  }
  atmos.vel_ew = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.vel_ew[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.vel_ew[i][j] = malloc(NTAU*sizeof(double));
    }
  }
  atmos.vel_ns = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.vel_ns[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.vel_ns[i][j] = malloc(NTAU*sizeof(double));
    }
  }
  atmos.vel_ve = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.vel_ve[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.vel_ve[i][j] = malloc(NTAU*sizeof(double));
    }
  }
  /* allocate memory for aero properties */

  /* MnS */
  atmos.aero_sw_tau_1 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_sw_tau_1[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_sw_tau_1[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_1 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_1[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_1[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_pi0_1 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_1[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_1[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  /* Al2O3 */
  atmos.aero_sw_tau_2 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_sw_tau_2[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_sw_tau_2[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_2 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_2[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_2[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_pi0_2 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_2[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_2[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  /* Fe */
  atmos.aero_sw_tau_3 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_sw_tau_3[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_sw_tau_3[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_3 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_3[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_3[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_pi0_3 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_3[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_3[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  /* MgSiO3 */
  atmos.aero_sw_tau_4 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_sw_tau_4[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_sw_tau_4[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_4 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_4[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_4[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.sw_pi0_4 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_4[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_4[i][j] = malloc(NTAU*sizeof(double));
      }
  }


  atmos.incident_frac = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.incident_frac[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.incident_frac[i][j] = malloc(NTAU*sizeof(double));
    }
  }
    
    
    file = fopen(T_P_3D_FILE, "r");
    if(file == NULL){
        printf("\nreadt_p.c:\nError opening file: No such file or directory\n\n");
        exit(1);
    }
    
    
    
    for (i=0; i<NLAT; i++){
        for(j=0; j<NLON; j++){
            for(k=0; k<NTAU; k++){
                if (fscanf(file, "%le %le %le %le %le %le %le %le %le %le\
                       %le %le %le %le %le %le %le %le %le %le %le %le",
                       &atmos.lat[i], &atmos.lon[j], &num, &atmos.alt[k], &atmos.P_3d[i][j][k],
                       &atmos.T_3d[i][j][k], &atmos.vel_ew[i][j][k], &atmos.vel_ns[i][j][k], &atmos.vel_ve[i][j][k],
                       &atmos.aero_sw_tau_1[i][j][k], &atmos.sw_asym_1[i][j][k], &atmos.sw_pi0_1[i][j][k],
                       &atmos.aero_sw_tau_2[i][j][k], &atmos.sw_asym_2[i][j][k], &atmos.sw_pi0_2[i][j][k],
                       &atmos.aero_sw_tau_3[i][j][k], &atmos.sw_asym_3[i][j][k], &atmos.sw_pi0_3[i][j][k],
                       &atmos.aero_sw_tau_4[i][j][k], &atmos.sw_asym_4[i][j][k], &atmos.sw_pi0_4[i][j][k],
                       &atmos.incident_frac[i][j][k]));
                
            }
            
            
            for(k=NTAU; k<NTAU; k++){
                if (fscanf(file, "%le %le %le %le %le\
                                  %le %le %le %le\
                                  %le %le %le %le\
                                  %le %le %le %le\
                                  %le %le %le %le\
                                  %le",
                                 &dum, &dum, &num, &dum, &dum,
                                 &dum, &dum, &dum, &dum,
                                 &dum, &dum, &dum, &dum,
                                 &dum, &dum, &dum, &dum,
                                 &dum, &dum, &dum, &dum,
                                 &dum));
            }
        }
    }
    
    
    printf("The first line of Values %le %le %le %le %le %le %le %le %le\n\
           %le %le %le\n\
           %le %le %le\n\
           %le %le %le\n\
           %le %le %le\n\
           %le\n",
           atmos.lat[0],
           atmos.lon[0],
           num,
           atmos.alt[0],
           atmos.P_3d[0][0][0],
           atmos.T_3d[0][0][0],
           atmos.vel_ew[0][0][0],
           atmos.vel_ns[0][0][0],
           atmos.vel_ve[0][0][0],
           atmos.aero_sw_tau_1[0][0][0],
           atmos.sw_asym_1[0][0][0],
           atmos.sw_pi0_1[0][0][0],
           atmos.aero_sw_tau_2[0][0][0],
           atmos.sw_asym_2[0][0][0],
           atmos.sw_pi0_2[0][0][0],
           atmos.aero_sw_tau_3[0][0][0],
           atmos.sw_asym_3[0][0][0],
           atmos.sw_pi0_3[0][0][0],
           atmos.aero_sw_tau_4[0][0][0],
           atmos.sw_asym_4[0][0][0],
           atmos.sw_pi0_4[0][0][0],
           atmos.incident_frac[0][0][0]);
    
    
    fclose(file);
    
}

/* ------- end ----------------------- ReadTP.c ------------------ */
