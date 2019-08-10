#include "dynamics.h"	

/*definition du systeme a integrer pour le calcul trajectoire, et section de poincaré*/
int func (double t, const double y[], double f[], void *params) {
	
	par_type *p = (par_type *) params;
	double u = p->u;
	double b01 = p->b01;
	double b02 = p->b02;
	double v = p->v;
	double s = p->s;
	double x = p->x;
	double e = p->e;
	double m1 = p->m1;
	double m2 = p->m2;

	/*systeme d'equa diff non lineaire*/

/*SS*/f[0] = u -(b01*(1+e*cos(2*M_PI*t)))*y[0]*(y[1]+x*y[3]+x*y[7]+m1) -(b02*(1+e*cos(2*M_PI*t)))*y[0]*(y[2]+x*y[3]+x*y[6]+m2) -u*y[0]; /*ok*/
/*IS*/f[1] = (b01*(1+e*cos(2*M_PI*t)))*y[0]*(y[1]+x*y[3]+x*y[7]+m1) -s*(b02*(1+e*cos(2*M_PI*t)))*y[1]*(y[2]+x*y[3]+x*y[6]+m2) -v*y[1] -u*y[1]; /*ok*/
/*SI*/f[2] = (b02*(1+e*cos(2*M_PI*t)))*y[0]*(y[2]+x*y[3]+x*y[6]+m2) -s*(b01*(1+e*cos(2*M_PI*t)))*y[2]*(y[1]+x*y[3]+x*y[7]+m1) -v*y[2] -u*y[2]; /*ok*/
/*II*/f[3] = s*(b02*(1+e*cos(2*M_PI*t)))*y[1]*(y[2]+x*y[3]+x*y[6]+m2)+ s*(b01*(1+e*cos(2*M_PI*t)))*y[2]*(y[1]+x*y[3]+x*y[7]+m1) -v*y[3]-v*y[3]-u*y[3]; /*ok*/
/*RS*/f[4] = v*y[1]-s*(b02*(1+e*cos(2*M_PI*t)))*y[4]*(y[2]+x*y[3]+x*y[6]+m2)-u*y[4]; /*ok*/
/*SR*/f[5] = v*y[2]-s*(b01*(1+e*cos(2*M_PI*t)))*y[5]*(y[1]+x*y[3]+x*y[7]+m1)-u*y[5]; /*ok*/
/*RI*/f[6] = s*(b02*(1+e*cos(2*M_PI*t)))*y[4]*(y[2]+x*y[3]+x*y[6]+m2) -v*y[6] +v*y[3] -u*y[6]; /*ok*/
/*IR*/f[7] = s*(b01*(1+e*cos(2*M_PI*t)))*y[5]*(y[1]+x*y[3]+x*y[7]+m1) -v*y[7] +v*y[3] -u*y[7]; /*ok*/
/*RR*/f[8] = v*y[6]+v*y[7]-u*y[8]; /*ok*/ 	

return GSL_SUCCESS;
}

