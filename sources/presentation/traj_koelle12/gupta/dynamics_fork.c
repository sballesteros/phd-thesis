/* gcc -Wall dynamics_fork.c func.c get_param.c integrator.c
   integrator_print.c int_poinc.c max_min.c seuil.c test.c failed.c
   -lm -lgsl -lgslcblas */

#include "dynamics.h"
int main (void) {
  /*acquisition des parametre et des conditions initiales*/
  double py[DIM];	
  get_param(py);
  param.m2=0.0;
  
  /*pointeur sur la structure de parametres. c'est ca qu'on passe aux fonction de GSL*/
  par_type *p_param = &param;
  
  /*vecteurs des variables d'états*/
  double y[DIM];
  
  /*on calcul le point d'equilibre de la souche residente */
  
#if FLAG_VERBOSE
  printf("equilibre\n");
#endif
  /*si il y a une erreure*/
  if(integrator(y, p_param, py, 0.0, T_EQ , 1e-9, 1e-9)){
#if FLAG_VERBOSE
    printf("erreur t'integration numerique on essaye avec prec 1e-14\n");
#endif
    if(integrator(y, p_param, py, 0.0, T_EQ , 1e-14, 1e-14)){
#if FLAG_VERBOSE
      printf("erreur t'integration numerique on essaye avec prec 1e-18\n");	
#endif
      if(integrator(y, p_param, py, 0.0, T_EQ , 1e-18, 1e-18)){
#if FLAG_VERBOSE
	printf("erreur t'integration numerique on essaye avec prec 1e-20\n");	
#endif
	if(integrator(y, p_param, py, 0.0, T_EQ , 1e-20, 1e-20)){
#if FLAG_VERBOSE
	  printf("on abandonne\n");
#endif
	  failed();
	  goto fin;
	}
      }
    }
  }
  
  
  int i;	
  /*on conserve les valeur des equilibres du resident comme nouvelle
    condition initiales*/
  for(i=0; i<DIM; i++){
    py[i]=y[i];
  }

  double abs_tol, rel_tol;
  /*on introduit le mutant et on regarde le transiant d'invasion*/
  param.m2=param.m1; 
  py[2]=0.000001;
  py[1]-=py[2];
  

  double traj_I1[NB_POINTS+1], traj_I2[NB_POINTS+1], traj_Itot[NB_POINTS+1];

#if FLAG_VERBOSE
  printf("on transiant d'invasion du mutant dans residant a l'equilibre\n");
#endif
  if(int_poinc(y, p_param,  py, 0.0, T_TRAJ, NB_POINTS, 1e-9, 1e-9 ,traj_I1, traj_I2, traj_Itot)){
    abs_tol=1e-14, rel_tol=1e-14;
#if FLAG_VERBOSE
    printf("erreur t'integration numerique on essaye avec prec 1e-14\n");
#endif
    if(int_poinc(y, p_param,  py, 0.0, T_TRAJ, NB_POINTS, 1e-14, 1e-14 ,traj_I1, traj_I2, traj_Itot)){
      abs_tol=1e-18, rel_tol=1e-18;
#if FLAG_VERBOSE
      printf("erreur t'integration numerique on essaye avec prec 1e-18\n");	
#endif
      if(int_poinc(y, p_param,  py, 0.0, T_TRAJ, NB_POINTS, 1e-18, 1e-18 ,traj_I1, traj_I2, traj_Itot)){
	abs_tol=1e-20, rel_tol=1e-20;
#if FLAG_VERBOSE
	printf("erreur t'integration numerique on essaye avec prec 1e-20\n");	
#endif
	if(int_poinc(y, p_param,  py, 0.0, T_TRAJ, NB_POINTS, 1e-20, 1e-20 ,traj_I1, traj_I2, traj_Itot)){
#if FLAG_VERBOSE
	  printf("on abandonne\n");
#endif
	  failed();
	  goto fin;
	}
      }
    }
  }
 
	
#if FLAG_VERBOSE
  printf("trans OK\n");
#endif

#if FLAG_BIF
  double *first_min, *first_max, * first_tmax;
  /*on garde les nb premiers max et min*/
  int nb=4;
	
  first_min = malloc(nb*(sizeof (double)));
  first_max = malloc(nb*(sizeof (double)));
  first_tmax= malloc(nb*(sizeof (double)));

  FILE *ppeak1; 
  ppeak1=fopen("peak1.dat","w");
	
  max_min_transiant(traj_I1, first_min, first_max, first_tmax, nb);
  fprintf(ppeak1, "%g\t", param.s);
  for(i=0;i<(nb-1);i++)
    {
      fprintf(ppeak1, "%g\t%g\t%g\t", first_min[i], first_max[i], first_tmax[i]);
    }
  fprintf(ppeak1, "%g\t%g\t%g\n", first_min[nb-1], first_max[nb-1], first_tmax[i]);
  fclose(ppeak1);
	
  FILE *ppeak2; 
  ppeak2=fopen("peak2.dat","w");
	
  max_min_transiant(traj_I2, first_min, first_max, first_tmax, nb);
  fprintf(ppeak2, "%g\t", param.s);
  for(i=0;i<(nb-1);i++)
    {
      fprintf(ppeak2, "%g\t%g\t%g\t", first_min[i], first_max[i], first_tmax[i]);
    }
  fprintf(ppeak2, "%g\t%g\t%g\n", first_min[nb-1], first_max[nb-1], first_tmax[i]);
  fclose(ppeak2);
	
  free(first_min);
  free(first_max);
  free(first_tmax);
#endif
	
#if FLAG_VERBOSE
  printf("max_min_transiant OK\n");
#endif

#if FLAG_SEUIL
  if(seuil_deter(y, p_param, py, 0.0, T_TRAJ , 1e-9, 1e-9)){
#if FLAG_VERBOSE
    printf("erreur t'integration numerique on essaye avec prec 1e-14\n");
#endif
    if(seuil_deter(y, p_param, py, 0.0, T_TRAJ , 1e-14, 1e-14)){
#if FLAG_VERBOSE
      printf("erreur t'integration numerique on essaye avec prec 1e-18\n");	
#endif
      if(seuil_deter(y, p_param, py, 0.0, T_TRAJ , 1e-18, 1e-18)){
#if FLAG_VERBOSE
	printf("erreur t'integration numerique on essaye avec prec 1e-20\n");	
#endif
	if(seuil_deter(y, p_param, py, 0.0, T_TRAJ , 1e-20, 1e-20)){
#if FLAG_VERBOSE
	  printf("on abandonne seuil\n");
#endif
	  failed();
	  goto fin;
	}
      }
    }
  }
#endif
	
 fin: 
  {}

  return 0;
}
