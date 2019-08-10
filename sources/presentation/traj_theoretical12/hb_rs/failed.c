#include "dynamics.h"


void failed(void) {
	#if FLAG_TRAJ
	FILE *ptraj; /*fichier des trajectoires*/
	ptraj=fopen("traj.dat","w");
	#endif

	#if FLAG_POINC
	FILE *ppoinc; 
	ppoinc=fopen("poinc.dat","w");
	#endif
	
	#if FLAG_BIF
	FILE *fbifmin1; 
	fbifmin1=fopen("bifmin1.dat","w");
	FILE *fbifmax1; 
	fbifmax1=fopen("bifmax1.dat","w");

	FILE *fphase1; 
	fphase1=fopen("phase1.dat","w");

	FILE *fsummary1; 
	fsummary1=fopen("summary1.dat","w");
	#endif

        #if FLAG_TRAJ
	fprintf (ptraj, PRINT_FAIL_TRAJ);
	fclose(ptraj); 
	#endif

	#if FLAG_POINC
	fprintf (ppoinc, PRINT_FAIL_POINC) ;
	fclose(ppoinc); 
	#endif
	
	#if FLAG_BIF				
	fprintf(fbifmin1, PRINT_FAIL_MIN);
	fprintf(fbifmax1, PRINT_FAIL_MAX);  
	fprintf(fphase1, PRINT_FAIL_PHASE);  

	fprintf(fsummary1, PRINT_FAIL_SUMMARY); 
	
	fclose(fbifmin1); 
	fclose(fbifmax1); 
	fclose(fphase1);	
	fclose(fsummary1);
	#endif

}


