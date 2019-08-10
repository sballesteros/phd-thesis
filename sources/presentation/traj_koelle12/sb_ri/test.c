#include "dynamics.h"


int test(double *y){
	int stop =0;
	int i;
	for (i=0; i<DIM; i++){
		if (isnan(y[i]) || y[i]<0) stop=1;
		//if (y[i]<0) stop=1;
	}

return stop;
}

