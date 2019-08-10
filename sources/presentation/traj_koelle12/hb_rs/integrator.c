#include "dynamics.h"

/*integration sans rien faire 
  retourne 1 si erreur d'integration numerique
  on met integrator_step car sinon ca merde*/

int integrator(double *y, par_type *pp, double y0[DIM], double t0, double t_end, double abs_tol, double rel_tol){

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

	  
  /*premiere tentative d'integration numerique avec une precision de ABS_TOL et REL_TOL*/
  while (t < t1){
    gsl_odeiv_evolve_apply (evolve, control, step, &sys, &t, t1, &h, y);
    if (test(y)) {
      error=1;
      break;
    }
  }
	  
  gsl_odeiv_step_free(step);
  gsl_odeiv_evolve_free(evolve);
  gsl_odeiv_control_free(control);
	
  return error;
}
