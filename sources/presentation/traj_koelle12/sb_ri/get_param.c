#include "dynamics.h"

void get_param(double *py) {

/*acquisition des parametres du systeme non lineaire*/
	scanf ("%lf", &param.u) ;
	scanf ("%lf", &param.b01) ;
	scanf ("%lf", &param.b02) ;
	scanf ("%lf", &param.v) ;
	scanf ("%lf", &param.s) ;
	scanf ("%lf", &param.e) ;
	scanf ("%lf", &param.m1) ;
	scanf ("%lf", &param.m2) ;
	
/*acquisition des conditions initiales du systeme non lineaire*/
	scanf ("%lf", &py[0]) ;
	scanf ("%lf", &py[1]) ;
	scanf ("%lf", &py[2]) ;
	scanf ("%lf", &py[3]) ;
	scanf ("%lf", &py[4]) ;
	scanf ("%lf", &py[5]) ;
}

