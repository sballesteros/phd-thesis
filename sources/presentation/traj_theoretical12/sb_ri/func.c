#include "dynamics.h"	

/*definition du systeme a integrer pour le calcul trajectoire, et section de poincaré*/
int func (double t, const double y[], double f[], void *params) {
	
	par_type *p = (par_type *) params;
	double u = p->u;
	double b01 = p->b01;
	double b02 = p->b02;
	double v = p->v;
	double s = p->s;
	double e = p->e;
	double m1 = p->m1;
	double m2 = p->m2;
	
	/*systeme d'equa diff non lineaire*/

	/*R0*/ f[0] = u -(b01*(1+e*cos(2*M_PI*t)))*y[0]*(y[4]+m1) -(b02*(1+e*cos(2*M_PI*t)))*y[0]*(y[5]+m2) -u*y[0];
	/*R1*/ f[1] = (1.0-s)*(b01*(1+e*cos(2*M_PI*t)))*y[0]*(y[4]+m1) -s*(b01*(1+e*cos(2*M_PI*t)))*y[1]*(y[4]+m1) -(b02*(1+e*cos(2*M_PI*t)))*y[1]*(y[5]+m2) -u*y[1];
	/*R2*/ f[2] = (1.0-s)*(b02*(1+e*cos(2*M_PI*t)))*y[0]*(y[5]+m2) -s*(b02*(1+e*cos(2*M_PI*t)))*y[2]*(y[5]+m2) -(b01*(1+e*cos(2*M_PI*t)))*y[2]*(y[4]+m1) -u*y[2];
	/*R12*/ f[3] =  s*(b01*(1+e*cos(2*M_PI*t)))*y[0]*(y[4]+m1) +s*(b02*(1+e*cos(2*M_PI*t)))*y[0]*(y[5]+m2) +s*(b01*(1+e*cos(2*M_PI*t)))*y[1]*(y[4]+m1) +s*(b02*(1+e*cos(2*M_PI*t)))*y[2]*(y[5]+m2) +(b02*(1+e*cos(2*M_PI*t)))*y[1]*(y[5]+m2) +(b01*(1+e*cos(2*M_PI*t)))*y[2]*(y[4]+m1) -u*y[3];
	/*I1*/ f[4] = (b01*(1+e*cos(2*M_PI*t)))*y[0]*(y[4]+m1) +(b01*(1+e*cos(2*M_PI*t)))*y[2]*(y[4]+m1) -v*y[4] -u*y[4];
	/*I2*/ f[5] = (b02*(1+e*cos(2*M_PI*t)))*y[0]*(y[5]+m2) +(b02*(1+e*cos(2*M_PI*t)))*y[1]*(y[5]+m2) -v*y[5] -u*y[5];

return GSL_SUCCESS;
}

