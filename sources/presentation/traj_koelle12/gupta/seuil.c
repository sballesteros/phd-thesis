#include "dynamics.h"

/*on doit etre sur que l'integration numerique fonctionne, a utiliser apres avoir trouver les abs_tol et rel_tol qui vont bien*/

int seuil_deter(double *y, par_type *pp, double y0[DIM], double t0, double t_end, double abs_tol, double rel_tol){

  int error=0;

  int i;
  /*allocation en memoire---------------------------------------------------------------------------*/  
  /*methode d'integration numerique*/
  const gsl_odeiv_step_type * T = METHODE;
	
  /*adaptive step size control*/
  gsl_odeiv_control * control = gsl_odeiv_control_y_new (abs_tol, rel_tol); /*abs and rel error (eps_abs et eps_rel) */
	
  gsl_odeiv_step * step = gsl_odeiv_step_alloc (T, DIM);
  /*evolution*/
  gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (DIM);

  /*on defini le systeme d'ODE, son jacobien et ses parametres*/
  gsl_odeiv_system sys = {func, jac, DIM, pp};

  /*borne d'integration*/
  double t = t0, t1 = t_end;
	
  /*pas d'incrementation initiale (modifie apres pour etre optimal selon fonction adaptive size*/
  double h = abs_tol;	

  /*on initialise le systeme non lineaire conditions initiales*/
  for(i=0; i<DIM; i++){
    y[i]=y0[i];
  }

  double seuil_I1_1 = -100.0, seuil_I1_2 = -100.0, seuil_I1_3 = -100.0;
  double seuil_I2_1 = -100.0, seuil_I2_2 = -100.0, seuil_I2_3 = -100.0;

  while (t < t1){
    gsl_odeiv_evolve_apply (evolve, control, step, &sys, &t, t1, &h, y);

    if (test(y)) {
      error=1;
      goto merde;
    }

    if(((y[1]+y[3]+y[7])<1e-6) && (seuil_I1_1 == -100.0)) seuil_I1_1=t;
    if(((y[2]+y[3]+y[6])<1e-6) && (seuil_I2_1 == -100.0)) seuil_I2_1=t;

    if(((y[1]+y[3]+y[7])<1e-7) && (seuil_I1_2 == -100.0)) seuil_I1_2=t;
    if(((y[2]+y[3]+y[6])<1e-7) && (seuil_I2_2 == -100.0)) seuil_I2_2=t;

    if(((y[1]+y[3]+y[7])<1e-8) && (seuil_I1_3 == -100.0)) seuil_I1_3=t;
    if(((y[2]+y[3]+y[6])<1e-8) && (seuil_I2_3 == -100.0)) seuil_I2_3=t;

  } /*fin du while*/



  FILE *pseuil; /*fichier des trajectoires*/
  pseuil=fopen("seuil.dat","w");

  if(seuil_I1_1 != -100.0){
    fprintf(pseuil,"%f\t%g\t",param.s, seuil_I1_1);
  }
  else fprintf(pseuil,"%f\tNA\t",param.s);
  if(seuil_I2_1 != -100.0){
    fprintf(pseuil,"%g\t", seuil_I2_1);
  }
  else fprintf(pseuil,"NA\t");


  if(seuil_I1_2 != -100.0){
    fprintf(pseuil,"%g\t", seuil_I1_2);
  }
  else fprintf(pseuil,"NA\t");
  if(seuil_I2_2 != -100.0){
    fprintf(pseuil,"%g\t", seuil_I2_2);
  }
  else fprintf(pseuil,"NA\t");

  if(seuil_I1_3 != -100.0){
    fprintf(pseuil,"%g\t", seuil_I1_3);
  }
  else fprintf(pseuil,"NA\t");
  if(seuil_I2_3 != -100.0){
    fprintf(pseuil,"%g\n", seuil_I2_3);
  }
  else fprintf(pseuil,"NA\n");	

  fclose(pseuil);

  gsl_odeiv_step_free(step);
  gsl_odeiv_evolve_free(evolve);
  gsl_odeiv_control_free(control);

 merde:
  {
  }
	

  return error;
}


