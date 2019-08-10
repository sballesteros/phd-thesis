#include "dynamics.h"

/*ATTENTION
modifier traj_I1 les trajectoire a stocker pour appliquer max_min
*/


/*borne d'integration, on dessine l'attracteur et on en prend une section de poincarré*/

int int_poinc(double *y, par_type *pp, double y0[DIM], double t0, double t_end, int nb_points, double abs_tol, double rel_tol, double *traj_I1 , double *traj_I2, double *traj_Itot){

	#if FLAG_TRAJ
	FILE *ptraj; /*fichier des trajectoires*/
	ptraj=fopen("traj.dat","w");
	#endif

	#if FLAG_POINC
	FILE *ppoinc; 
	ppoinc=fopen("poinc.dat","w");
	int k;
	#endif

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

	#if FLAG_POINC
	double yo[DIM], yp[DIM];
	#endif

	for (i = 0; i <= nb_points; i++) {
		double ti = i * t1 / nb_points;
		while (t < ti){

			#if FLAG_POINC
			/*On conserve l'ancienne valeur */
			double to=t;
			for (k=0; k<DIM; k++){
		     		yo[k] = y[k] ;
		     	}
			#endif

			/*On integre d'un pas de temps */    		
			gsl_odeiv_evolve_apply (evolve, control, step, &sys, &t, ti, &h, y);

			#if FLAG_POINC
			/*Calcul de la section de poincare*/
			if ((cos(2*M_PI*t) > 0.0) && (cos(2*M_PI*to) < 0.0)){
				double alph = -cos(2*M_PI*to)/(cos(2*M_PI*t)-cos(2*M_PI*to)) ;
				for(k=0; k<DIM; k++){
					yp[k] = yo[k] + alph * (y[k]-yo[k]) ;
				}
				fprintf (ppoinc, PRINT_POINC) ;
			}
			#endif
	
		} /*fin du while*/
			
		#if FLAG_BIF
		/*on stock le signal à analyser pour la determination des extremas*/
		traj_I1[i] = y[4];
		traj_I2[i] = y[5];
		traj_Itot[i] = y[4]+y[5];
		#endif 
		
		#if FLAG_TRAJ
		fprintf (ptraj, PRINT_TRAJ);
		#endif

	} /*fin du for sur i NB_POINTS*/

	#if FLAG_TRAJ
	fclose(ptraj);
	#endif

	#if FLAG_POINC
	fclose(ppoinc);
	#endif

	gsl_odeiv_step_free(step);
	gsl_odeiv_evolve_free(evolve);
	gsl_odeiv_control_free(control);

return 0;
}


