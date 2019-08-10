/*
Analyse systeme dynamique 
last modified mar 20 mai 2008 15:56:12 CEST by Sebastien Ballesteros

-affichage des trajectoires
-section de poincaré
-diagramme de bifurcation (max min phase)
-estimation UPCA
-calcul du spectre d'exposant de lyapunov si FLAG_LYAP different de 0
adapté d'apres "determining lyapunov exponents from a time series"
Alan WOLF et al. physica D 1985 285-317 


resolution de systeme d'ODE avec GSL

These sources and GSL are 100% free software
These sources are released in the public domain.
You are free to edit, re-distribute, modify, delete, send me cash,
print it on T-shirts, sing it, learn the code, port it to another
language, whatever. Enjoy!

Sebastien Ballesteros
*/

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_math.h>


#define FLAG_VERBOSE 0
#define FLAG_SEUIL 0 /*si vaut 0 pas de calcul des exposants de lyapunov*/
#define FLAG_TRAJ 1 /*si vaut 0 pas de stockage des trajectoires*/
#define FLAG_POINC 0 /*si vaut 0 pas de section de poincaré*/
#define FLAG_BIF 0 /*si vaut 0 pas de bifurcation et de phases*/


/*on defini les constantes---------------------------------------------------------------------------*/
	#define DIM 9 /*dimension du systeme non lineaire a integrer*/
	#define METHODE gsl_odeiv_step_rkf45

	/*doit etre inferieur a 1e-9*/
	#define ABS_TOL 1e-6
	#define REL_TOL 1e-6
	#define DT0 1e-6

	#define T_EQ 10000.0 /*temps d'integration de la dynamique transiante*/
	#define T_TRAJ 30.0 /*temps d'integration de la trajectoire qu'on garde*/

	#define NB_POINTS 3000

	#define PRECISION 1e-6
	#define CON 11 /*doit etre un nombre impaire !!!!!*/

/*parametrages des sorties*/
	
//	#define PRINT_TRAJ "%f\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", t, y[0], y[1], y[2], y[3], y[4], y[5]
	#define PRINT_TRAJ "%f\t%g\t%g\n", t, y[1]+y[3]+y[7], y[2]+y[3]+y[6]
	#define PRINT_POINC	"%f\t%g\t%g\n", param.s, yp[1]+yp[3]+yp[7], yp[2]+yp[3]+yp[6]	

 	/*sortie de la fonction max_min_transiant*/
	#define PRINT_TRANS_MIN_EQ	"%f\t%g\tNA\n", param.s, (traj_min+traj_max)/ ((double) 2.0)
	#define PRINT_TRANS_MAX_EQ	"%f\t%g\tNA\n", param.s, (traj_min+traj_max)/ ((double) 2.0)
	#define PRINT_TRANS_MAX	"%f\t%g\t%g\n", param.s, max, t_max
	#define PRINT_TRANS_MIN	"%f\t%g\t%g\n", param.s, min, t_min


	/*on regarde si on est sur un point d'equilibre*/
	#define PRINT_MIN_EQ	"%f\t%f\n", param.s, (traj_min+traj_max)/ ((double) 2.0)
	#define PRINT_MAX_EQ	"%f\t%f\n", param.s, (traj_min+traj_max)/ ((double) 2.0) 
	#define PRINT_PHASE_EQ	"%f\t0.0\n", param.s
	#define PRINT_SUMMARY_EQ	"%f\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\t0\n", param.s

	/*si il y a des extrema*/ 
	#define PRINT_PHASE	"%f\t%g\n", param.s, phase	
	#define PRINT_MAX	"%f\t%g\n", param.s, max
	#define PRINT_MIN	"%f\t%g\n", param.s, min
	#define PRINT_SUMMARY	"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n", param.s, traj_min, mean_min, sd_min, maxmin, minmax, mean_max, sd_max, traj_max, phasemin, mean_T, sd_T, phasemax, categorie(input,nb_max)
											
						
	/*si erreur d'integration numerique*/
	#define PRINT_FAIL_TRAJ	"NaN\tNaN\tNaN\n"

	#define PRINT_FAIL_POINC	"%f\tNaN\tNaN\tNaN\n", param.s

	#define PRINT_FAIL_MIN	"%f\tNaN\n", param.s
	#define PRINT_FAIL_MAX	"%f\tNaN\n", param.s 
	#define PRINT_FAIL_PHASE	"%f\tNaN\n", param.s  
	#define PRINT_FAIL_SUMMARY	"%f\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n", param.s  

/*---------------------------------------------------------------------------on defini les constantes*/


/*variables globales----------------------------------------------------------------------------------------*/

/*ce bloc doit etre avant la declaration des prototypes de fonctions sinon il ne reconait pas les par_type *pp */
/*definition des parametres du modele*/
typedef struct {
  double u;
  double b01;
  double b02;
  double v;
  double s;
  double e;
  double m1;
  double m2;
  //		double to; /*pour orbit d'invasion*/
} par_type;

	par_type param;

/*----------------------------------------------------------------------------------------variables globales*/

/*prototypes de fonctions------------------------------------------------------------------------------------------*/  

#if FLAG_BIF
double b_max(double blocc[CON]);
double b_min(double blocc[CON]);
double get_min(double tab[], int dim);
double get_max(double tab[], int dim);
void max_min(double t_traj[NB_POINTS+1], int indice);
void max_min_transiant(double t_traj[NB_POINTS+1], double *first_min, double *first_max, double *first_tmax, int nb);
#endif

#if FLAG_PLOT
void plot_it(void);
#endif


int func(double t, const double y[], double f[], void *params);
void *jac;

int integrator(double *y,  par_type *pp, double y0[DIM], double t0, double t_end, double abs_tol, double rel_tol);
int integrator_print(double *y, par_type *pp, double y0[DIM], double t0, double t_end, int nb_point, double abs_tol, double rel_tol);
int int_poinc(double *y, par_type *pp, double y0[DIM], double t0, double t_end, int nb_points, double abs_tol, double rel_tol,  double *traj_I1 , double *traj_I2, double *traj_Itot);
	
int seuil_deter(double *y, par_type *pp, double y0[DIM], double t0, double t_end, double abs_tol, double rel_tol);

int test(double *y);
void get_param(double *py);
void failed(void);
	
/*------------------------------------------------------------------------------------------prototypes de fonctions*/  

