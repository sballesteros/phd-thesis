#include "dynamics.h"

#if FLAG_BIF
void max_min(double t_traj[NB_POINTS+1], int indice)
{
  /*indice sert dans le cas ou on a plusieurs serie a analyse (typiquement I1, I2 et Itot*/
  char bifmin[20];
  sprintf(bifmin,"bifmin%d.dat",indice);
  char bifmax[20];
  sprintf(bifmax,"bifmax%d.dat",indice);

  FILE *fmin; 
  fmin=fopen(bifmin,"w");
  FILE *fmax; 
  fmax=fopen(bifmax,"w");

  int i, k, stop1, stop2;
  double traj_min, traj_max, t_max;
  double max, min;
  /*on initialise un bloc à analyser par la fonction max_ou_min*/
  double bloc[CON];

  /*on regarde si on est sur un point d'equilibre*/
  traj_min=get_min(t_traj, NB_POINTS);
  traj_max=get_max(t_traj, NB_POINTS);

  if(fabs(traj_max - traj_min) < PRECISION){
    fprintf(fmin, "%g\t%g\t%g\t%g\n", param.e, param.b01, param.v, (traj_min+traj_max)/ ((double) 2.0));
    fprintf(fmax, "%g\t%g\t%g\t%g\tNaN\n", param.e, param.b01, param.v, (traj_min+traj_max)/ ((double) 2.0));  
    goto fin; /*on zap la determination des min et max*/
  }
  else 
    { /*si il y a un extrama*/
      stop1=0, stop2=1;
      for(i=0;  i <= (NB_POINTS-CON) ; i++)
	{
	  /*on initialise un bloc à analyser par la fonction max_ou_min. à optimiser !!! on recopie trop!!!!!*/
	  for (k=0 ; k<CON ; k++)
	    {
	      bloc[k]=t_traj[i+k];	
	    }
	  
	  max=b_max(bloc);
	  if(max != -100.0 && stop1==0)
	    { /*si il y a un max*/ 
	      t_max=(i+(CON-1)/2)*T_TRAJ/((double) NB_POINTS);
	      fprintf(fmax, "%g\t%g\t%g\t%g\t%g\n", param.e, param.b01, param.v, max, t_max);
	      stop1=1;
	      stop2=0;
	    }

	  min=b_min(bloc);	
	  if(min != -100.0 && stop2==0)
	    {/*si il y a un min*/ 
	      fprintf(fmin, "%g\t%g\t%g\t%g\n", param.e, param.b01, param.v, min);
	      stop1=0;
	      stop2=1;
	    }

	} /*fin du for*/
    } /*fin du else*/

 fin : 
  {}
	
  fclose(fmin); 
  fclose(fmax); 
}


/*ATTENTION INTEGRE SUT T_TRAJ*/
/*determine les nb premier max et min*/
void max_min_transiant(double t_traj[NB_POINTS+1], double *first_min, double *first_max, double *first_tmax, int nb)
{
  int i, k, l, stop1, stop2;
  double traj_min, traj_max;
  double max, min;
  /*bloc à analyser par la fonction max_ou_min*/
  double bloc[CON];

  /*on regarde si on est sur un point d'equilibre*/
  traj_min=get_min(t_traj, NB_POINTS);
  traj_max=get_max(t_traj, NB_POINTS);
	
  if(fabs(traj_max - traj_min) < PRECISION)
    {
      for(i=0;i<(nb-1);i++)
	{
	  first_min[i] = (traj_min+traj_max)/ ((double) 2.0);
	  first_max[i] = (traj_min+traj_max)/ ((double) 2.0);
	  first_tmax[i] = 0.0;
	}
      goto fin; /*on zap la determination des min et max*/
    }
  else
    { /*si il y a un extrema*/
      stop1=0, stop2=1;
      l=0;
      for(i=0;  i <= (NB_POINTS-CON) ; i++)
	{
	  /*on initialise un bloc à analyser par la fonction max_ou_min. à optimiser !!! on recopie trop!!!!!*/
	  for (k=0 ; k<CON ; k++)
	    {
	      bloc[k]=t_traj[i+k];	
	    }
	  
	  max=b_max(bloc);
	  if(max != -100.0 && stop1==0)
	    { /*si il y a un max*/ 
	      first_tmax[l] =((double)(i+(CON-1)/2.0)*T_TRAJ)/((double) NB_POINTS);
	      first_max[l] = max;
	      stop1=1;
	      stop2=0;
	    }

	  min=b_min(bloc);	
	  if(min != -100.0 && stop2==0)
	    {/*si il y a un min*/ 
	      first_min[l] = min;
	      l+=1;
	      if (l==nb) goto fin;
	      stop1=0;
	      stop2=1;
	    }

	} /*fin du for*/
    } /*fin du else*/
  
 fin : 
  {}
}

/*fonction qui regarde si l'element centrale du tableau est le max et
  qui le retourne si c'est le cas et retourne -100 sinon on rajoute
  une condition pour prendre que certain max ici ceux plus grand que
  le seuil endemique
*/
double b_max(double blocc[CON])
{
  int i;
  double max;
	
  for (i=0 ; i<(CON-1)/2 ; i++) {
    if (blocc[i]>blocc[i+1]){
      max=-100.0;
      goto fin;	
    }
  }
	
  for (i=(CON-1)/2 ; i<(CON-1);  i++) {
    if (blocc[i]<blocc[i+1]) {
      max=-100.0;
      goto fin;	
    }
  }	

  max = blocc[(CON-1)/2];
  double eq_end= 1.2*(param.u/(param.u+param.v)- param.u/param.b01);
  if (max<eq_end)
    {
      max=-100.0;
    }
 fin : 
  {}
  return max;
}

/*fonction qui regarde si l'element centrale du tableau est le min et qui le retourne*/
double b_min(double blocc[CON])
{
  int i;
  double min;		
  for (i=0 ; i<(CON-1)/2 ; i++) {
    if (blocc[i]<blocc[i+1]) {
      min=-100.0;
      goto fin;	
    }
  }
	
  for (i=(CON-1)/2 ; i<(CON-1);  i++) {
    if (blocc[i]>blocc[i+1])  {
      min=-100.0;
      goto fin;	
    }	
  }	

  min = blocc[(CON-1)/2];

 fin : 
  {}
  return min;
}

double get_min(double tab[], int dim)
{
  int i;
  double min=tab[0];
  for (i=0;i<=dim ;i++){
    if (tab[i]<min) {
      min=tab[i];
    }
  }
  return (min);
}

double get_max(double tab[], int dim){
  int i;
  double max=tab[0];
  for (i=0;i<=dim ;i++){
    if (tab[i]>max){
      max=tab[i];
    }
  }
  return (max);
}
#endif
