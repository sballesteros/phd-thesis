//gcc -Wall sirs_migr.c -lm -lgsl -lgslcblas
/*Algorithme de Gillespie
  les clusters sont modelise par des modele SIR
  on suppose cluster 1 a l'eq, le cluster 2 apparait
  on regarde la dynamique d'invasion
  code source valide*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#define RAND gsl_rng_uniform_pos(randgsl)

#define FLAG_TRAJ 1 /*affichage des trajectoire ON 1 OFF 0*/
#define FLAG_TEXT 0 
#define FLAG_VERBOSE 1 

#define NB_REACTION 4
#define NB_ETAT 3
#define TMAX 200
#define AFFICHE_MAX 100000
#define NB_REALISATION 1

/*declaration des parametres constants CSS de 6M verifier par simulation 1000 simuls 0 ext sur 100 ans*/
double s;

gsl_rng *randgsl; /* global random number generator */
   
/*--------------------------------------------------------------------------------------*/
/*reactions :
  S->S+1, R->R-1 gR
  S+I->I+I bSI/N
  S+n->I bSn/N
  I->R vI
*/
 
/*matrice des reactions (N P colone, reaction lignes)*/
const double reaction[NB_REACTION][NB_ETAT] = {
  /*S  I  R*/
  {1.0,0.0,-1.0},
  {-1.0,1.0,0.0},
  {-1.0,1.0,0.0},
  {0.0,-1.0,1.0},
};

/*--------------------------------------------------------------------------------------*/

int main (void) {

  /*declaration des parametres constants*/
  const double N = 100000.0 ; 
  const double v = 365.0/7.0;
  const double b0 =5.0*v;
  const double g = 0.1;
  const double e = 0.0;
  const double n = 0.01;
  
  /*declaration et initialisation du generateur de nombres aleatoires*/
  const gsl_rng_type *Type; /*type de generateur aleatoire*/
  unsigned long int seed; /*SEED*/
  Type = gsl_rng_mt19937; /*MT19937 generator of Makoto Matsumoto and Takuji Nishimura*/
  seed = (unsigned) time(NULL); /*seed avec l'heure locale*/
  randgsl = gsl_rng_alloc(Type); /*allocation en memoire*/
  gsl_rng_set(randgsl, seed); /*on seed le rng*/
	
  /*fichiers de sorties*/
#if FLAG_TRAJ
  FILE *ptrajectoire; /*fichier des trajectoires*/
  ptrajectoire=fopen("traj.dat","w");
#endif
#if FLAG_TEXT
  FILE *ptemps_extinction; /*fichier des trajectoires*/
  ptemps_extinction=fopen("temps_extinction.dat","w");
#endif
	
  /*variable de calcul*/
  double tau;
  double temps;
  double calc_a0, calc_b0, k;
  int indice_reaction, realisation;
  int j;
	
#if FLAG_TRAJ
  int i;
  double temps_affiche[AFFICHE_MAX];
  /* on calcul les temps que l'on veut plotter*/
  temps_affiche[0]=0.0;
  for (i=0 ; i< (AFFICHE_MAX-1) ; i++){
    temps_affiche[i+1] = temps_affiche[i] + ((double) TMAX)/ ((double) AFFICHE_MAX);
    /*on va de 0 à  (AFFICHE_MAX-1)*(TMAX/ AFFICHE_MAX)*/
  }	
#endif
	
  for (realisation=1 ; realisation <= NB_REALISATION ; realisation++){
		
#if FLAG_VERBOSE		
    printf("realisation %d sur %d\n",realisation, NB_REALISATION);
#endif
		
    /*variable de calcul a reinitialiser a chaque realisations*/
    tau = 0.0;
    temps = 0.0;
		
#if FLAG_TRAJ
    i=0; /*control du temps a ploter*/
#endif
		
    double text1=0.0; /*temps extiction*/
    int stop1=0;	/*stop1 sert pour avoir temps d'extinction*/	
		
    /*----------------------------------------------------------*/
    /*conditions initiales des variables d'etat */
//    double S = N-10.0;
//    double I = 10.0;
//    double R = N-S-I;

    double S = (double) (long) ((v/b0)*N);
    double I = (double) (long) (((g-g*v/b0)/(v+g))*N);
    double R = N-S-I;
    printf("condition initiale : %g\t%g\t%g\n", S, I, R);
    /*----------------------------------------------------------*/
    /* algorithme de GILLESPIE (1976) exact */
    while (temps <= TMAX) 
      {
			
#if FLAG_TRAJ
      if (temps >= (temps_affiche[i])){	
	fprintf(ptrajectoire,"%f\t%g\t%g\t%g\n", temps, S,I,R);
	i += 1;
      }
#endif
	
      /*on cree c et h (a=c.*h) et on defini calc_a0 */
      double a[NB_REACTION];
      
      /*--------------------------------------------------------------------------------------*/
      a[0] = R*g; 
      a[1] = S*(I+n)*b0*(1.0+e*cos(2.0*M_PI*temps))/(S+I+R);
      a[2] = S*n*b0*(1.0+e*cos(2.0*M_PI*temps))/(S+I+R);
      a[3] = I*v;
      /*--------------------------------------------------------------------------------------*/

      /*on fabrique calc_a0*/
      calc_a0=0;
      for (j = 0 ; j < NB_REACTION ; j++)
	calc_a0 += a[j];
			
      tau = -(log(RAND))/calc_a0; /*temps de la prochaine reactions*/
				
      indice_reaction=0;
      calc_b0=a[0]; /*quel est la prochaine reaction ?*/
      k=RAND*calc_a0;
      while (calc_b0 < k) 
	calc_b0 += a[++indice_reaction];
					
      /*--------------------------------------------------------------------------------------*/
      S += reaction[indice_reaction][0];
      I += reaction[indice_reaction][1];
      R += reaction[indice_reaction][2];
      /*--------------------------------------------------------------------------------------*/
			
      temps += tau;	

      if ((I == 0) && stop1<1)
	{text1 = temps;
	  stop1=1;}
				
    } /*fin du while*/
		
#if FLAG_TRAJ
    fprintf(ptrajectoire,"\n\n");	
#endif
		
#if FLAG_TEXT
    fprintf(ptemps_extinction,"%f\t%f\n", temps, text1);
#endif
  } /*fin du for sur realisation*/
					
#if FLAG_TRAJ
  fclose(ptrajectoire);
#endif

#if FLAG_TEXT	
  fclose(ptemps_extinction);	
#endif

  gsl_rng_free(randgsl);

  return 0;
}
