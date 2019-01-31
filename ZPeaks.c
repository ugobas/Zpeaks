#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Peaks_aux.h"
#include "HMM_aux.h"
#include "HMM.h"
#include "allocate.h"
#include "cluster_score.h"
#include "optimization.h"

int CONTR=0;
int RANDOM=1;
int OMIT_MITO=1; // Omit mitochondrion and chloroplast?

int NORM_PEAK_SCORE=1; // Normalize the score of each peak dividing by length?

#define CODE "ZPeaks"
#define NCHRMAX 80   // Maximum number of chromosomes
#define NPROF 100 // Maximum number of profiles to read
#define VERBOSE 1

char NAMEOUT[80]="Peaks";
#define EXT ".bed"

// Control variables
int STAT_COMP=1;    // Make statistics of comparison?
int PRSCORE=1;      // Print scores in file?

// Global parameters
char *file_e[NCHRMAX], *file_c[NCHRMAX];
int DCLUST=200, DTOL=50;
int SIZE_INI=200; // For removing small fragments
int SIZE_END=300;
int SIZE_STEP=50;
int SIZE_MIN=100;
float THR=2.5;
int WIN=50;

// HMM parameters
int ncl=2;
int *numclus, *numdom;
char namelog[200];
FILE *filelog;

// data
int Nchr;
long nn, *nnch; long *x_coord;
float *y_contr, *y_exper;
short *chr; char **chr_name;
int *Omit_chr=NULL;

// Comparison
struct peak *peak_ref=NULL;
char *file_ref=NULL;
int N_peak_ref=0;

// Profiles
int N_prof;
int nchr_prof[NPROF];
char  **file_prof[NPROF], *name_prof[NPROF];


int Get_peaks(int *cluster, float *y_scr, int *WIN,
	      float Damp, float Thr, float l_EPS, int S_MIN);
struct peak *Clusters2Frags(int *N_peak, int *cluster, float *y_scr, long nn);
struct peak *Center_peaks(int N_peak, struct peak *peaks,
			  float *y_scr, long *x_coord, long nn);
float Smooth(float *y_box,
	     float *y, long N, float DAMP, int *WIN, float l_EPS);
void Smooth_chromosomes(float *y_scr, long nn,
			float Damp, int *WIN, float l_EPS);

float Cluster_likelihood(int *cluster, float *y, int nn, int ncl);
int Get_clusters(int *cluster, float Thr, float *y_scr, int SIZE_MIN);
float Get_likelihood(float Thr, int *cluster, float *y, long nn, int SIZE_MIN);
static void Cluster_size(int *numclus, int *numdom, int *cluster, short *chr,
			 int N,int ncl);
static int Set_control(float *yc, long nn, float Mean);
static int Get_input(int *CONTR, char **file_c, char **file_e,
		     // control and experiment
		     char **file_ref, // reference peaks
		     char ***file_prof, int *N_prof, // profiles
		     char **name_prof, int *nchr_prof,
		     float *THR, int *SIZE_MIN, int *DCLUST,
		     int *PEAKSIZE, int *DTOL, // parameters
		     char *NAMEOUT, char *MODEL, 
		     int argc, char **argv);
static void help();

static void Copy_prof(long *nprof, long **xprof, float **yprof,
		      long *x, float *y, long *nnch, int Nchr);
void Comp_Z_score(float *Z_scr, long nn, float *ye_w, float *yc_w, float sd);
static void Rescale_counts(float *yy, long nn, float scale);
static float Normalize_counts(float *y, long n, float AVE);
static float Mean_counts(float *sd, float *yy, long nn);
static void Split_chr(long *nnch, short *chr, char **chr_name,
		      long nn, int Nchr);
static void Print_score(long *x, float *y, short *chr, long nn,
			char **chr_name, char *NAMEOUT);
static int Read_contr_exp(int CONTR, long **x_coord, float **y_contr,
			  float **y_exper, short **chr, long *nn,
			  char ***chr_name, char *file_c, char *file_e,
			  char *outp);
void Examine_peaks(struct peak *peaks, int N_peak,
		   struct peak *peak_ref, int N_peak_ref, int Nchr,
		   double *Prof_peak, double **Prof_score,
		   double *Discriminant_score,
		   long ***xprof, float ***yprof, long **nprof,
		   int *nchr_prof, char **name_prof, int N_prof,
		   char *NAMEOUT, int WINDOW, int DCLUST,
		   int SIZE_MIN, float Thr, char *what);
int Close_peak(struct peak **peak, long x_coord, float y);
void Compute_scores(int *cluster, short *chr, float **x_VarSam,
		    int Nsam, int Nvar, int ncl, float lik,
		    struct Para Par, FILE *filelog);

//float *Mean_control, *Mean_exp, *sd_control, *sd_exp;
char MODEL[80];
int BOX1; // Size of the wig files

int main(int argc, char **argv)
{
  ZSCORE=1; // Print Z score of properties (1) or raw properties (0)?
  //strcpy(MODEL,"E"); // E=exponential G=Gaussian
  int p, k, j; long i;

  // Read input
  for(p=0; p<NPROF; p++){
    file_prof[p]=malloc(NCHRMAX*sizeof(char *));
    for(k=0; k<NCHRMAX; k++)file_prof[p][k]=NULL;
    nchr_prof[p]=0;
  }
  Get_input(&CONTR, file_c, file_e, &file_ref,
	    file_prof, &N_prof, name_prof, nchr_prof, &THR, &SIZE_MIN,
	    &DCLUST, &PEAKSIZE, &DTOL, NAMEOUT, MODEL, argc, argv);

  // Read experiment and control
  Nchr=Read_contr_exp(CONTR, &x_coord, &y_contr, &y_exper, &chr,
		      &nn, &chr_name, file_c[0], file_e[0], "Statistics.dat");
  BOX1=x_coord[1]-x_coord[0];
  printf("First bin: %d %d\n", x_coord[0], x_coord[1]-1);
  printf("Box size= %d nucleotides %d chromosomes\n", BOX1, Nchr);

  nnch=malloc(Nchr*sizeof(long));
  Split_chr(nnch, chr, chr_name, nn, Nchr);
  
  if(OMIT_MITO){
    printf("Eliminating mitochondria and chloroplast, if any\n");
    Omit_chr=malloc(Nchr*sizeof(int));
    int no=0, imax=0; long nk=0;
    for(i=0; i<Nchr; i++){
      if((strncmp(chr_name[i], "mito", 4)==0)||
	 (strncmp(chr_name[i], "chloro", 6)==0)){
	Omit_chr[i]=1; no++;
      }else{
	Omit_chr[i]=0;
	imax=i+1; nk+=nnch[i];
      }
    }
    printf("%d chromosomes eliminated, %d remaining\n", no, Nchr-no);
    if((imax+no)==Nchr){Nchr=imax; nn=nk;}
    else{
      printf("WARNING, the eliminated chromosomes are not the last ones\n");
    }
  }

  // Read profiles, if any
  N_prof++;
  long  **nprof=malloc(N_prof*sizeof(long *));
  long  ***xprof=malloc(N_prof*sizeof(long **));
  float ***yprof=malloc(N_prof*sizeof(float **));
  p=N_prof-1;
  nchr_prof[p]=Nchr;
  name_prof[p]=malloc(100*sizeof(char));
  strcpy(name_prof[p], "Peakscore");
  for(p=0; p<N_prof; p++){
    int nchr=nchr_prof[p], Mchr=nchr; if(Nchr>nchr)Mchr=nchr;
    printf("Reading profile %s %d files\n", name_prof[p], nchr_prof[p]);
    nprof[p]=malloc(Mchr*sizeof(long));
    xprof[p]=malloc(Mchr*sizeof(long *));
    yprof[p]=malloc(Mchr*sizeof(float *));
    if(p<(N_prof-1)){
      for(k=0; k<nchr; k++){
	nprof[p][k]=Read_file(&xprof[p][k], &yprof[p][k], file_prof[p][k]);
	printf("Reading %d lines %s\n", nprof[p][k], file_prof[p][k]);
      }
      for(k=nchr; k<Mchr; k++)nprof[p][k]=0;
    }
  }
  printf("%d profiles stored\n", N_prof);

  // Read reference peaks, if any
  char namefound[200];
  char **lab_p=NULL;
  if(file_ref!=NULL){
    printf("Reading %s\n", file_ref);
    peak_ref=Read_peaks(&N_peak_ref, file_ref, -1); // Chr=chr-1
    printf("%d reference peaks\n", N_peak_ref);
    sprintf(namefound, "%s_found%s", file_ref, EXT);
  }

  // Normalize to one read per nucleotide
  float AVE=1.0;
  // Normalize both control and experiment per each chromosome
  if(CONTR==0){
    y_contr=malloc(nn*sizeof(float));
    Set_control(y_contr, nn, AVE);
    printf("Using uniform control\n");
  }else{
    printf("Using control\n");
  }
  printf("%d chromosomes\n", Nchr);
  long n0=0;
  for(i=0; i<Nchr; i++){
    long ni=nnch[i], nl=0;
    float sd=Normalize_counts(y_exper+n0, ni, AVE);
    if(CONTR){
      sd=Normalize_counts(y_contr+n0, ni, AVE);
      // Increase controls lower than the average
      for(j=n0; j<n0+ni; j++){
	if(y_contr[j]<AVE){y_contr[j]=0.5*(y_contr[j]+AVE); nl++;}
      }
    }
    n0+=ni;
    printf("%ld control bins over %ld (%.3f) smoothed\n",
	   nl, n0, (float)nl/n0);
  }

  sprintf(namelog,"%s_ZP_log.dat", NAMEOUT);
  filelog=fopen(namelog, "w");

  int *cluster=malloc(nn*sizeof(int));
  float *y_scr=malloc(nn*sizeof(float));
  numdom=malloc(2*sizeof(int));
  numclus=malloc(2*sizeof(int));

  float DAMP_INI=0.1, DAMP_STEP=40, Damp_max=DAMP_INI; int N_damp=10;
  float EPS, EPS_INI=0.02, EPS_END=0.101, EPS_STEP=0.04, EPS_max=EPS_INI;
  float S_max=0, Thr=THR;
  for(EPS=EPS_INI; EPS<=EPS_END; EPS+=EPS_STEP){
    float l_EPS=-log(EPS);
    
    float D1=DAMP_INI,
      S1=Get_peaks(cluster, y_scr, &WIN, D1, Thr, l_EPS, SIZE_MIN);
    if(S1>S_max){S_max=S1; Damp_max=D1; EPS_max=EPS;}
    
    float D2=D1+DAMP_STEP,
      S2=Get_peaks(cluster, y_scr, &WIN, D2, Thr, l_EPS, SIZE_MIN);
    if(S2>S_max){S_max=S2; Damp_max=D2; EPS_max=EPS;}
      
    float D3=D2+DAMP_STEP,
      S3=Get_peaks(cluster, y_scr, &WIN, D3, Thr, l_EPS, SIZE_MIN);
    if(S3>S_max){S_max=S3; Damp_max=D3; EPS_max=EPS;}
      
    float S0_old=0; int kd;
    for(kd=0; kd<N_damp; kd++){
      float D0=Find_max_quad(D1, D2, D3, S1, S2, S3, 1, 1000);
      if(isnan(D0))break;
      float S0=Get_peaks(cluster, y_scr, &WIN, D0, Thr, l_EPS, SIZE_MIN);
      if(S0>S_max){S_max=S0; Damp_max=D0; EPS_max=EPS;}
      if((S0<S0_old)||((S0-S0_old)< 10))break;
      S0_old=S0;
      Rearrange_points(&D1, &D2, &D3, &S1, &S2, &S3, D0, S0);
    } // End loop of Damping
  } // End loop of EPS
  printf("Best: Damp=%.0f EPS=%.2g Np= %.0f (%.3f of all)\n",
	 Damp_max, EPS_max, S_max, S_max/nn);

  float l_EPS=-log(EPS_max);
  int N_frag=Get_peaks(cluster, y_scr, &WIN, Damp_max, THR, l_EPS, SIZE_MIN);
  printf("Best: WIN= %d Damp=%.0f Thr=%.2f Np= %d\n",
	 WIN, Damp_max, THR, N_frag);
  
  /*
  // Find optimal threshold: it is too small, it is not a good strategy
  float L_max=-100, T_opt=THR, T_STEP=0.2;
  float T1=THR-T_STEP, L1=Get_likelihood(T1, cluster, y_scr, nn, SIZE_MIN);
  if(L1>L_max){L_max=L1; T_opt=T1;}
  float T2=T1+T_STEP, L2=Get_likelihood(T2, cluster, y_scr, nn, SIZE_MIN);
  if(L2>L_max){L_max=L2; T_opt=T2;}
  float T3=T2+T_STEP, L3=Get_likelihood(T3, cluster, y_scr, nn, SIZE_MIN);
  if(L3>L_max){L_max=L3; T_opt=T3;}
  float T0_old=0, L0_old=L_max; int Nt=40, it;
  for(it=0; it<Nt; it++){
    float T0=Find_max_quad(T1, T2, T3, L1, L2, L3, 0.5, 5.0);
    if(isnan(T0))break;
    float L0=Get_likelihood(T0, cluster, y_scr, nn, SIZE_MIN);
    if(L0>L_max){L_max=L0; T_opt=T0;}
    if((L0<L0_old)||((L0-L0_old)< 0.00001))break;
    L0_old=L0;
    Rearrange_points(&T1, &T2, &T3, &L1, &L2, &L3, T0, L0);
  } // End loop of Thr
  Thr=T_opt;
  N_frag=Get_peaks(cluster, y_scr, &WIN, Damp_max, Thr, l_EPS, SIZE_MIN);
  printf("After optimizing threshold:\n");
  printf("Best: WIN= %d Damp=%.0f Thr=%.2f Np= %d\n",
	 WIN, Damp_max, Thr, N_frag);
  */

  // Record optimal score
  if(PRSCORE)Print_score(x_coord, y_scr, chr, nn, chr_name, NAMEOUT);
  p=N_prof-1;
  for(k=0; k<Nchr; k++){
    nprof[p][k]=nnch[k];
    xprof[p][k]=malloc(nnch[k]*sizeof(long));
    yprof[p][k]=malloc(nnch[k]*sizeof(float));
  }
  Copy_prof(nprof[p], xprof[p], yprof[p], x_coord, y_scr, nnch, Nchr);

  // Transform into fragments
  int N_peak;
  struct peak *peaks=Clusters2Frags(&N_peak, cluster, y_scr, nn);

  // Set xo, x1, x2, size (next is not set)
  struct peak *peak=peaks;
  while(peak != NULL){Set_peak(peak, 0); peak=peak->next;}

  // Join clusters closer than DCLUST: It worsens performances
  //N_peak=Make_cluster_check(peaks, N_peak, DCLUST, x_coord, y_scr, nn, 0, 0);
  //printf("After joining domains with d<%d: %d\n", DCLUST, N_peak);

  //float *y1=malloc(nn*sizeof(int)); int WW; float DAMP=1;
  //Smooth_chromosomes(y1, nn, DAMP, &WW, l_EPS);


  int Size_Min=SIZE_MIN;
  /* for(Size_Min=SIZE_INI; Size_Min<=SIZE_END; Size_Min+=SIZE_STEP){

    // Eliminate small fragments
    printf("\n\nS=%d\n",  Size_Min);
    if(Size_Min >= BOX1){
      N_peak=Remove_fragments(peaks, N_peak, Size_Min);
      printf("After removing small domains: %d\n", N_peak);
    }
    if(N_peak<2*Nchr){
      printf("ERROR, too few peaks found: %d\nExiting\n", N_peak); 
      break;
      }*/

  // Allocate variables for statistical analysis
  int N_peak_max=N_peak;
  if(N_peak_ref>N_peak_max){N_peak_max=N_peak_ref;}
  double *Discriminant_score=malloc(N_prof*sizeof(double));
  double **Prof_score=malloc(N_prof*sizeof(double *));
  double *Prof_peak=malloc(N_prof*sizeof(double));
  for(p=0; p<N_prof; p++)Prof_score[p]=malloc(N_peak_max*sizeof(double));

  char string[1000];
  sprintf(string, "# S_min Npeaks size dist ");
  for(p=0; p<N_prof; p++)sprintf(string, "%s %s", string, name_prof[p]);
  fprintf(filelog, "%s\n", string);
  if(N_peak_ref)
    sprintf(string,
	    "%s d_12 d_21 over_1 over_2 prod o1/ran o2/ran prod/ran", string);

  if(N_peak==0){
    printf("WARNING, no peaks have been found\n");
    return(0);
  }

  int WINDOW=(2*WIN+1)*BOX1;
  // Peaks are printed here!
  // Center peaks in such a way that the center coincides with the maximum
  // of the not-smoothed score and the peak is included in the previous one
  struct peak *cpeaks=Center_peaks(N_peak, peaks, y_scr, x_coord, nn); //y_scr
  printf("Examine centered peaks\n");
  Examine_peaks(cpeaks, N_peak, peak_ref, N_peak_ref, Nchr, //centered
		Prof_peak, Prof_score, Discriminant_score,
		xprof, yprof, nprof, nchr_prof, name_prof, N_prof,
		NAMEOUT, WINDOW, DCLUST, Size_Min, Thr, "");

  printf("\nExamine not centered peaks\n");
  Examine_peaks(peaks, N_peak, peak_ref, N_peak_ref, Nchr, //not centered
		Prof_peak, Prof_score, Discriminant_score,
		xprof, yprof, nprof, nchr_prof, name_prof, N_prof,
		NAMEOUT, WINDOW, DCLUST, Size_Min, Thr, "NotCentered_");



  return(0);
}

void Examine_peaks(struct peak *peaks, int N_peak,
		   struct peak *peak_ref, int N_peak_ref, int Nchr,
		   double *Prof_peak, double **Prof_score,
		   double *Discriminant_score,
		   long ***xprof, float ***yprof, long **nprof,
		   int *nchr_prof, char **name_prof, int N_prof,
		   char *NAMEOUT, int WINDOW, int DCLUST,
		   int SIZE_MIN, float Thr, char *what)

{
  // Statistics of peaks
  char header[200]; int p;
  sprintf(header, "\n# size_min num size dist");
  for(p=0; p<N_prof; p++)sprintf(header, "%s %s", header, name_prof[p]);
  printf("%s\n", header);

  char string[1000]; int PRINT=1;
  Profile_score(Prof_peak, Prof_score, Discriminant_score,
		xprof, yprof, nprof, nchr_prof, N_prof, peaks, N_peak);
  char nameout[200]; sprintf(nameout, "Properties_%s.dat", NAMEOUT);
  sprintf(string, "%d ", SIZE_MIN);
  Sizediff_stat(peaks, Prof_peak, Prof_score, Discriminant_score,
		N_prof, name_prof, nameout, PRINT, string);
  fprintf(filelog, "%s\n", string);


  if(RANDOM){
    // Random peaks
    struct peak *peakran=Extract_peaks(peaks, N_peak, Nchr); //NULL
    Profile_score(Prof_peak,Prof_score,Discriminant_score,
		  xprof,yprof,nprof,nchr_prof,N_prof,peakran,N_peak);
    sprintf(string, "%d ", SIZE_MIN);
    Sizediff_stat(peakran, Prof_peak, Prof_score, Discriminant_score,
		  N_prof, name_prof, "Properties_Random.dat", 1, string);
    free(peakran);
  }
  
  // Compare with references, if any
  if(N_peak_ref){
    float over_1, over_2, over_1r, over_2r;
    int n_match=0, n_match_ref=0;
    /*float dist=Dist_peaks(peaks,N_peak,peak_ref,N_peak_ref);
    printf(" %.0f", dist);
    dist=Dist_peaks(peak_ref,N_peak_ref,peaks,N_peak);
    printf(" %.0f", dist); */

    // Overlap with random peaks
    struct peak *peak2=Extract_peaks(peak_ref, N_peak_ref, Nchr);
    Count_matches(&n_match,&n_match_ref,DTOL,peaks,N_peak,peak2,N_peak_ref);
    over_1r=(float)n_match_ref/N_peak_ref;
    over_2r=(float)n_match/N_peak;
    free(peak2);

    // Overlap with reference
    Count_matches(&n_match,&n_match_ref,DTOL,peaks,N_peak,peak_ref,
		  N_peak_ref);
    over_1=(float)n_match_ref/N_peak_ref;
    over_2=(float)n_match/N_peak;

    printf("  %.3f %.3f %.3f",  over_1, over_2, over_1*over_2);
    printf("  %.1f %.1f %.1f",  over_1/over_1r, over_2/over_2r,
	   over_1*over_2/(over_1r*over_2r));
  }
  printf("\n");

  if(PRINT){

    // Print peaks
    char Pars[100];
    sprintf(Pars, "W%d_T%.2f_J%d_S%d", WINDOW, Thr, DCLUST, SIZE_MIN);
    sprintf(nameout, "%s_ALL_%s%s.bed", NAMEOUT, what, Pars);
    Print_Peaks(peaks, nameout);

    if(N_peak_ref){
      sprintf(nameout, "%s_notfound.bed", file_ref);
      // Print_Peaks_nomatch(peak_ref, N_peak_ref, nameout);

      char nameold[200];
      sprintf(nameold, "%s_OLD_%s.bed", NAMEOUT, Pars);
      sprintf(nameout, "%s_NEW_%s.bed", NAMEOUT, Pars);
      // Print_Peaks_new(peaks, nameout, nameold);

      // Statistics of peaks present / not present in previous set
      int Np2=0;
      int N_peak_max=N_peak; if(N_peak_ref>N_peak_max)N_peak_max=N_peak_ref;
      struct peak *peak2=malloc(N_peak_max*sizeof(struct peak));

      printf("%s\n", header);

      Np2=Select_peaks_match(peak2, peaks, N_peak, 0); // not matched
      Profile_score(Prof_peak, Prof_score, Discriminant_score,
		    xprof, yprof, nprof, nchr_prof, N_prof, peak2, Np2);
      sprintf(string, "%d ", SIZE_MIN);
      Sizediff_stat(peak2, Prof_peak, Prof_score, Discriminant_score,
		    N_prof, name_prof, "Properties_NewPeaks.dat",PRINT,string);


      Np2=Select_peaks_match(peak2, peaks, N_peak, 1); // matched
      Profile_score(Prof_peak, Prof_score, Discriminant_score,
		    xprof, yprof, nprof, nchr_prof, N_prof, peak2, Np2);
      sprintf(nameout, "Properties_Common%s.dat", NAMEOUT);
      sprintf(string, "%d ", SIZE_MIN);
      Sizediff_stat(peak2, Prof_peak, Prof_score, Discriminant_score,
		    N_prof, name_prof, nameout, PRINT, string);

      Np2=Select_peaks_match(peak2, peak_ref, N_peak_ref, 0); // not matched
      Profile_score(Prof_peak, Prof_score, Discriminant_score,
		    xprof, yprof, nprof, nchr_prof, N_prof, peak2, Np2);
      sprintf(string, "%d ", SIZE_MIN);
      Sizediff_stat(peak2, Prof_peak, Prof_score, Discriminant_score, N_prof,
		    name_prof, "Properties_UnconfirmedPeaks.dat",PRINT,string);

      Profile_score(Prof_peak, Prof_score, Discriminant_score,
		    xprof, yprof,nprof,nchr_prof, N_prof,peak_ref, N_peak_ref);
      sprintf(string, "%d ", SIZE_MIN);
      Sizediff_stat(peak_ref, Prof_peak, Prof_score,
		    Discriminant_score, N_prof, name_prof,
		    "Properties_ReferencePeaks.dat", 1, string);
      free(peak2);
    }
    printf("\n");


    /************* Metaplots ******************/
    sprintf(nameout, "Metaplots_%s_S%d.dat", NAMEOUT, SIZE_MIN);
    printf("\n");
    for(p=0; p<N_prof; p++){
      Plot_profile(peaks, p, xprof[p], yprof[p], nprof[p], nchr_prof[p],
		   name_prof[p], nameout);
      printf("Metaplot of %s\n", name_prof[p]);
    }
    if(N_peak_ref){
      for(p=0; p<N_prof; p++){
	Plot_profile(peak_ref, p, xprof[p], yprof[p], nprof[p], nchr_prof[p],
		     name_prof[p], "Metaplots_reference.dat");
      }
    }

  }
  
}

int Get_input(int *CONTR, char **file_c, char **file_e, // control and exper
	      char **file_ref, // reference peaks
	      char ***file_prof, int *N_prof, // profiles
	      char **name_prof, int *nfile_prof,
	      float *THR, int *SIZE_MIN, int *DCLUST,
	      int *PEAKSIZE, int *DTOL, // parameters
	      char *NAMEOUT, char *MODEL, 
	      int argc, char **argv)      // input
{
  int i;
  if(argc<2)help();
  for(i=1; i<argc; i++){
    if(strncmp(argv[i], "-h", 2)==0)help();
    if(i>1)printf("WARNING, unrecognized option %s\n", argv[i]);
  }

  // Default
  *PEAKSIZE=PEAKSIZE_DEF;

  // Reading parameter file
  char file[NCHAR], dumm[80];
  strcpy(file, argv[1]);
  printf("Reading %s\n", file);
  FILE *file_in=fopen(file, "r");
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  char string[1000];
  int ne=0, nc=0, np=0, n1=0, npf=0;
  int re=0, rc=0, rp=0, r1=0;
  char PATH[NCHAR]; PATH[0]='\0';
  double norm, NSAM; int  wd; float x;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'){
      continue; // Comment line
    }else if(strncmp(string, "END", 3)==0){
      rp=0; r1=0; re=0; rc=0; PATH[0]='\0';
      if(n1){nfile_prof[npf]=n1; npf++; n1=0;}
    }else if(strncmp(string, "DIR", 3)==0){
      char *s=string+4; int k=0;
      while(*s!='\n'){PATH[k]=*s; k++; s++;}
      PATH[k]='\0';
    }else if(strncmp(string, "EXPER", 5)==0){
      re=1;
    }else if(re){
      Read_name(&ne, file_e, PATH, string, NCHRMAX);
    }else if(strncmp(string, "CONTROL:", 8)==0){
      rc=1; 
    }else if(rc){
      Read_name(&nc, file_c, PATH, string, NCHRMAX);
    }else if(strncmp(string, "PREDICTION", 8)==0){
      rp=1;
    }else if(rp){
      if(np){
	printf("ERROR, only one reference file allowed\n"); exit(8);
      }
      Read_name(&np, file_ref, PATH, string, NCHRMAX);
    }else if(strncmp(string, "PROF", 4)==0){
      r1=1;
      name_prof[npf]=malloc(NCHAR*sizeof(char));
      sscanf(string+5, "%s", name_prof[npf]);
    }else if(r1){
      Read_name(&n1, file_prof[npf], PATH, string, NCHRMAX);
    }else if(strncmp(string, "NAME", 4)==0){
      sscanf(string+5, "%s", NAMEOUT);
    }else if(strncmp(string, "PEAKSIZE=", 8)==0){
      sscanf(string+8, "%d", PEAKSIZE);
      printf("Peak size set to %d\n", *PEAKSIZE);
    }else if(strncmp(string, "SIZE_MIN", 8)==0){
      sscanf(string+9, "%d", SIZE_MIN);
      printf("Minimum size required for calling a peak= %d\n", *SIZE_MIN);
      /*}else if(strncmp(string, "SIZE_END", 8)==0){
      sscanf(string+9, "%d", SIZE_END);
      printf("Minimum size required (largest)= %d\n", *SIZE_END);
    }else if(strncmp(string, "SIZE_STEP", 9)==0){
      sscanf(string+10, "%d", SIZE_STEP);
      printf("Minimum size required (step)= %d\n", *SIZE_STEP);*/
    }else if(strncmp(string, "DCLUST=", 7)==0){
      sscanf(string+7, "%d", DCLUST);
    }else if(strncmp(string, "DTOL=", 5)==0){
      sscanf(string+5, "%d", DTOL);
      printf("DTOL= %d\n", *DTOL);
    }else if(strncmp(string, "THR=", 4)==0){
      sscanf(string+4, "%f", THR);
      printf("THR= %.2f\n", *THR);
    }else if(strncmp(string, "MODEL", 5)==0){
      sscanf(string+6, "%s", dumm);
      if((strcmp(dumm, "E")!=0)&&(strcmp(dumm, "G")!=0)){
	printf("WARNING, model %s not implemented\n", dumm);
	printf("ALLowed models are E (exponential) and G (Gaussian)\n");
	printf("Default is %s\n", MODEL);
      }else{
	printf("Model changed from %s to %s\n", MODEL, dumm);
	strcpy(MODEL, dumm);
      }
    }else if(string[0]!='\n'){
      printf("WARNING, unrecognized line: %s\n", string);
    }
  }
  *N_prof=npf;
  if((ne==0)){
    printf("ERROR, no input files specified\n"); help();
  }
  if((ne>1)||(nc>1)){
    printf("ERROR, only one wig file allowed for experiment and control,");
    printf(" found %d and %d\n", ne, nc); exit(8);
  }
  if((nc==0)){
    printf("WARNING, no control files specified\n");
    printf("Using mean experiment box as control\n");
    *CONTR=0;
  }else{
    *CONTR=1;
  }

  printf("Executing %s\n", CODE);
  //printf("PEAKSIZE= %d\n", *PEAKSIZE);
  return(ne);  
}

void help(){
  printf("\nPROGRAM %s\n", CODE);
  printf("Author: Ugo Bastolla,\n");
  printf("Centro de Biologia Molecular Severo Ochoa\n");
  printf("(CSIC-UAM), Madrid Spain\n");
  printf("<ubastolla@cbm.csic.es>\n\n");
  printf("Mandatory argument: parameter file\n");
  printf("FORMAT:\n");
  printf("EXPER:\n");
  printf("DIR=<path of exper files>\n");
  printf("<exper file> (one for each chromosome, MANDATORY)\n");
  printf("END\n");
  printf("CONTROL:\n");
  printf("DIR=<path of control files>\n");
  printf("<control file> (one for each chromosome, optional)\n");
  printf("END\n");
  printf("PREDICTION:\n");
  printf("DIR=<path of prediction file>\n");
  printf("<prediction file> (only one, optional)\n");
  printf("END\n");
  printf("PROF <profile name>:\n");
  printf("DIR=<path of prof results>\n");
  printf("<prof file> (one for each chromosome, optional)\n");
  printf("END\n");
  printf("(Same for all experimental profiles)\n");
  //printf("WINDOW=<#size of window for averaging experiment and control>\n");
  printf("THR=<Threshold for positives> (default %.2f)\n", THR);
  printf("SIZE_MIN=<Minimum size for calling a peak>\n");
  printf("DCLUST=<Distance threshold for joining fragments>\n");
  /*printf("SIZE_INI=<Minimum of minimum size for calling a peak>\n");
  printf("SIZE_END=<Maximum of minimum size for calling a peak>\n");
  printf("SIZE_INI=<Step of minimum size for calling a peak>\n");*/
  printf("DTOL=<Tolerance for comparison>\n");
  printf("MODEL=E ! Score distribution used in the HMM\n");
  printf("# Allowed: E (exponential) and G (Gaussian)\n");
  printf("\n"); exit(8);
}

void Copy_prof(long *nprof, long **xprof, float **yprof,
	       long *x_scr, float *y_scr, long *nnch, int Nchr)
{
  int k; int i;
  long *xs=x_scr; float *ys=y_scr;
  for(k=0; k<Nchr; k++){
    long *xp=xprof[k]; float *yp=yprof[k];
    for(i=0; i<nnch[k]; i++){
      *xp=*xs; xp++; xs++;
      *yp=*ys; yp++; ys++;
    }
  }
}

void Split_chr(long *nnch, short *chr, char **chr_name, long nn, int Nchr){
  long i, n=0; short k, ichr=-1;
  for(i=0; i<nn; i++){
    if(chr[i]!=ichr){
      if(ichr>=0)nnch[ichr]=n;
      n=0; ichr++;
    }
    n++;
  }
  nnch[ichr]=n;  ichr++;
  if(ichr != Nchr){
    printf("ERROR in Split_chr, number of chromosomes= %d expected= %d\n",
	   ichr, Nchr); exit(8);
  }
  for(k=0; k<ichr; k++){
    printf("Chromosome %s %ld fragments\n", chr_name[k], nnch[k]);
  }
}

void Comp_Z_score(float *Z_scr, long nn, float *ye_box, float *yc_box, float sd)
{
  double Z1=0, Z2=0; long i;
  float *ycb=yc_box, *yeb=ye_box, *Z=Z_scr, yc;
  for(i=0; i<nn; i++){
    float ZZ=(*yeb-*ycb);
    //float yc=0.5*(*ycb+1.0), ZZ=(*yeb-yc);
    //float ZZ=(*yeb-*ycb)/(*ycb);
    // WARNING, the last normalization presents a bias towards AT rich
    *Z=ZZ; Z1+=ZZ; Z2+=ZZ*ZZ;
    Z++; yeb++; ycb++;
  }
  Z1/=nn; Z2=sqrt(Z2/nn-Z1*Z1);
  if(Z2>0){
    for(i=0; i<nn; i++)Z_scr[i]=(Z_scr[i]-Z1)/Z2;
  }
}

int Set_control(float *y, long nn, float Mean){
  long i; float *yi=y;
  for(i=0; i<nn; i++){*yi=Mean; yi++;}
  return(0);
}

float Normalize_counts(float *y, long n, float AVE){
  float sd, Mean=Mean_counts(&sd, y, n);
  Rescale_counts(y, n, AVE/Mean);
  return(sd*AVE/Mean);
}

float Mean_counts(float *sd, float *yy, long nn){
  double Y1=0, Y2=0; float *y=yy; long i;
  for(i=0; i<nn; i++){
    Y1+=*y; Y2+=(*y)*(*y); y++;
  }
  Y1/=nn; Y2=sqrt((Y2-nn*Y1*Y1)/(nn-1));
  *sd=Y2;
  return(Y1);
}

void Rescale_counts(float *yy, long nn, float scale){
  long i; float *y=yy;
  for(i=0; i<nn; i++){*yy *= scale; yy++;}
}

void Print_score(long *x, float *y, short *chr, long nn,
		 char **chr_name, char *NAMEOUT)
{
  char nameout[200];
  sprintf(nameout, "%s_score.wig", NAMEOUT);
  printf("Writing %s\n", nameout);
  FILE *file_out=fopen(nameout, "w");
  fprintf(file_out, "track type=wiggle_0\n");
  int step=x[1]-x[0];
  int k=-1; long i;
  for(i=0; i<nn; i++){
    if(chr[i]!=k){
      fprintf(file_out, "fixedStep chrom=%s start=1 step=%d span=%d\n",
	      chr_name[chr[i]], step, step);
      k=chr[i];
    }
    fprintf(file_out, "%.4f\n", y[i]);
  }
  fclose(file_out);
}

int Read_contr_exp(int CONTR, long **x_coord, float **y_contr,
		   float **y_exper, short **chr, long *nn, char ***chr_name,
		   char *file_c, char *file_e, char *outp)
{
  int k; long i;
  FILE *file_out=fopen(outp, "w"); 
  int nch1[NCHRMAX], nch2[NCHRMAX], nch[NCHRMAX];
  short *chr2=NULL; long *x2=NULL; *y_contr=NULL;
  int step1, step2;


  printf("Reading %s\n", file_e);
  int Nchr, Nch1, Nch2;
  Nch1=Read_chroms_wig(nch1, chr_name, &step1, file_e, NCHRMAX);
  for(k=0; k<Nch1; k++)nch[k]=nch1[k];
  printf("file %s %d chromosomes read\n",file_e, Nch1);
  Nchr=Nch1;
  if(CONTR){
    Nch2=Read_chroms_wig(nch2, chr_name, &step2, file_c, NCHRMAX);
    if(step2!=step1){
      printf("ERROR, different step\n");
      printf("Experiment: %s %d chromosomes step=%d\n",file_e,Nch1,step1);
      printf("Control:    %s %d chromosomes step=%d\n",file_c,Nch2,step2);
      exit(8);
    }
    if(Nch2 != Nch1){
      printf("WARNING, different number of chromosomes\n");
      printf("Experiment: %s %d chromosomes step=%d\n",file_e,Nchr,step1);
      printf("Control:    %s %d chromosomes step=%d\n",file_c,Nch2,step2);
      if(Nch2>Nch1){Nchr=Nch2; for(k=Nch1; k<Nch2; k++)nch[k]=0;}
      printf("Setting number of chromosomes to %d\n", Nchr);
      long diff=0;
      for(k=0; k<Nchr; k++)diff+=abs(nch2[k]-nch[k]);
      printf("Difference of number of lines: %ld\n", diff);
      if(diff > abs(Nch1-Nch2)*10000){
	printf("Too many differences, exiting\n"); exit(8);
      }
    }
    for(k=0; k<Nchr; k++)if(nch2[k]>nch[k])nch[k]=nch2[k];
  }
  (*nn)=0; for(k=0; k<Nchr; k++)(*nn)+=nch[k];

  Read_wig_2(x_coord, y_exper, chr, file_e, nch, Nch1);
  if(CONTR){
    Read_wig_2(&x2, y_contr, chr, file_c, nch, Nch2);
    float slope, offset,
      r=Corr_coeff(&slope,&offset,(*y_contr),(*y_exper),(*nn));
    fprintf(file_out, "Corr(exper,control)= %.3f\n", r);
    free(x2); free(chr2);
  }
  fclose(file_out);
  return(Nchr);
}

void Cluster_size(int *numclus, int *numdom, int *cluster, short *chr,
		  int N, int ncl)
{
  int i, k, k0=-1, ch=-1;
  for(k=0; k<ncl; k++){numclus[k]=0; numdom[k]=0;}
  for(i=0; i<N; i++){
    k=cluster[i];
    if((k<0)||(k>=ncl)){
      printf("WARNING, wrong cluster identifier %d (%d clusters expected)\n",
	     k, ncl); continue;
    }
    if(chr[i]!=ch){ch=chr[i]; k0=-1;}
    if(k!=k0){numdom[k]++; k0=k;}
    numclus[k]++;
  }
}

int Close_peak(struct peak **peak, long x_coord, float y)
{
  (*peak)->end=x_coord-1;
  (*peak)->y=y;
  (*peak)->next=(*peak)+1;
  (*peak)++;
  return(1);
}


void Compute_scores(int *cluster, short *chr, float **x_VarSam,
		    int Nsam, int Nvar, int ncl, float lik,
		    struct Para Par, FILE *filelog)
{
  int numclus[ncl], numdom[ncl], k;
  Cluster_size(numclus, numdom, cluster, chr, Nsam, ncl);
  printf("Cluster sizes (peaks/no peaks): ");
  for(k=0; k<ncl; k++)printf(" %d", numdom[k]); printf(" domains ");
  for(k=0; k<ncl; k++)printf(" %d", numclus[k]); printf(" elements\n");

  // Compute likelihood and other scores
 float cscore=Cluster_score(cluster, ncl, x_VarSam, Nsam, Nvar);
 // Effective number of variables
 float lcorr=Get_lcorr(x_VarSam, Nsam, 0);
 float Nsam_eff=Nsam;
 if((lcorr > 0)&&(lcorr < Nsam))Nsam_eff/=lcorr;
 float norm_lik1=Nsam_eff/Nsam, norm_lik2=Nsam_eff;
 lik*=norm_lik1;
 // Compute number of parameters
 int N_para = Nvar + 1; // mu, tau (per cluster)
 if(Par.sig){N_para+= Nvar*(Nvar+1)/2;}
 else if(Par.scale_pos){
   N_para+=Nvar;
   if(Rel_diff(Par.scale_neg[0][0],Par.scale_pos[0][0])>0.02)N_para++;
   if(Rel_diff(Par.scale_neg[1][0],Par.scale_pos[1][0])>0.02)N_para++;
 }
 if(Par.trans)N_para+=(ncl-1);
 N_para=ncl*N_para-1;

 float aic=AIC(lik, N_para, Nsam_eff);
 float bic=BIC(lik, N_para, Nsam_eff);
 lik/=norm_lik2; aic/=norm_lik2; bic/=norm_lik2;
 // Separation score
 double chi=-1;
 if(Par.sig){
   chi=Chi2(Par.mu[0], Par.mu[1], Par.sig[0], Par.sig[1], Nvar);
 }else if(Par.scale_pos){
   float sig0=max(Par.scale_pos[0][0], Par.scale_neg[0][0]);
   float sig1=max(Par.scale_pos[1][0], Par.scale_neg[1][0]);
   chi=Chi2(Par.mu[0], Par.mu[1], &sig0, &sig1, 1);
 }

 // Print
 char txt[400];
 sprintf(txt,
	 "# %d clusters: lik= %.4f AIC=%.4f BIC= %.4f (/NPC) score=%.3f\n",
	 ncl, lik, aic, bic, cscore);
 sprintf(txt,"%s# Discriminative power chi2= %.3f\n", txt, chi);
 printf("%s\n", txt);
 fprintf(filelog, "%s", txt);
 Print_parameters_f(numclus, &Par, ncl, Nvar, filelog);
}

int Get_peaks(int *cluster, float *y_scr, int *WIN,
	      float Damp, float Thr, float l_EPS, int SIZE_MIN)
{
  // Smooth over window, compute Z score and assign clusters
  Smooth_chromosomes(y_scr, nn, Damp, WIN, l_EPS);

  // Remove small domains
  int Np=Get_clusters(cluster, Thr, y_scr, SIZE_MIN);
  printf("Thr=%.2f W=%d D=%.0f N=%d (S>%d)\n",
	 Thr, *WIN, Damp, Np, SIZE_MIN);
  return(Np);
}

void Smooth_chromosomes(float *y_scr, long nn,
			float Damp, int *WIN, float l_EPS)
{
  float *ye_w=malloc(nn*sizeof(float)), *yc_w;
  if(CONTR){yc_w=malloc(nn*sizeof(float));}
  else{yc_w=y_contr;}
  long n0=0; int i;
  for(i=0; i<Nchr; i++){
    long ni=nnch[i];
    if((Omit_chr)&&(Omit_chr[i])){n0+=ni; continue;}
    float sd=Smooth(ye_w+n0, y_exper+n0, ni, Damp, WIN, l_EPS);
    if(CONTR)sd=Smooth(yc_w+n0, y_contr+n0, ni, Damp, WIN, l_EPS);
    Comp_Z_score(y_scr+n0, ni, ye_w+n0, yc_w+n0, sd);
    n0+=ni;
  }
  free(ye_w); if(CONTR)free(yc_w);
}

float Smooth(float *y_box,
	     float *y, long N, float DAMP, int *WIN, float l_EPS)
{
  /* Weighted average of the number of reads in neighboring boxes
   */
  double y1=0, y2=0;

  float DampFact=exp(-BOX1/DAMP);
  *WIN=l_EPS*DAMP/BOX1;

  int k, ic; long i, j;
  float *Damp=malloc((*WIN+1)*sizeof(float));
  Damp[0]=1.0;
  for(k=1; k<= *WIN; k++){
    Damp[k]=Damp[k-1]*DampFact;
  }

  float *y_new=y_box, *y_old=y;
  for(i=0; i<N; i++){
    float ww=1, ysum=(*y_old), *w=Damp+1, *y_ptr=y_old+1;
    j=i+1;
    for(k=1; k<=*WIN; k++){
      if(j>=N)break;
      ysum+=(*y_ptr)*(*w); ww+=(*w); w++; j++; y_ptr++;
    }
    w=Damp+1; y_ptr=y_old-1; j=i-1; 
    for(k=1; k<=*WIN; k++){
      if(j<0)break;
      ysum+=(*y_ptr)*(*w); ww+=(*w); w++; j--; y_ptr--;
    }
    *y_new=(ysum/ww);
    y1+=(*y_new); y2+=(*y_new)*(*y_new);
    y_new++; y_old++;
  }
  free(Damp);
  y1/=N; y2=(y2-N*y1*y1)/(N-1);
  return(sqrt(y2));
}

float Get_likelihood(float Thr, int *cluster, float *y, long nn, int SIZE_MIN)
{
  int Np=Get_clusters(cluster, Thr, y, SIZE_MIN);
  float lik=Cluster_likelihood(cluster, y, nn, 2)/nn;
  printf("Thr=%.2f likelihood/N: %.4f\n", Thr, lik);
  return(lik);
}

int Get_clusters(int *cluster, float Thr, float *y_scr, int SIZE_MIN)
{
  int Np=0, i; 
  int S=SIZE_MIN/BOX1, m; long n0=0, j;
  for(i=0; i<Nchr; i++){
    long ni=nnch[i];
    if((Omit_chr)&&(Omit_chr[i])){
      for(j=n0; j<(n0+ni); j++)cluster[j]=0;
      n0+=ni; continue;
    }
    for(j=n0; j<(n0+ni); j++){
      if(y_scr[j]>Thr){Np++; cluster[j]=1;}
      else{cluster[j]=0;}
    }
    int k0=-1, k, j0=n0;
    for(j=n0; j<n0+ni; j++){
      k=cluster[j];
      if(k!=k0){
	if(k==0){ // Peak ends, check whether size >= S
	  if((j-j0)<S)
	    for(m=j0; m<j; m++){cluster[m]=0; Np--;}
	}else{ // Peak starts
	  j0=j;
	}
	k0=k;
      }
    }
    if((k==1)&&((j-j0)<S))
      for(m=j0; m<j; m++){cluster[m]=0; Np--;}
    n0+=ni;
  }
  return(Np);
}

float Cluster_likelihood(int *cluster, float *y, int nn, int ncl)
{
  printf("Computing likelihood\n");
  double *y1=malloc(ncl*sizeof(double));
  double *y2=malloc(ncl*sizeof(double));
  double *ll=malloc(ncl*sizeof(double));
  int *nc=malloc(ncl*sizeof(int)), i, k;
  for(k=0; k<ncl; k++){y1[k]=0; y2[k]=0; nc[k]=0;}
  for(i=0; i<nn; i++){
    k=cluster[i];
    if((k<0)||(k>=ncl))printf("ERROR, cluster %d out of bound\n", k);
    y1[k]+=y[i]; y2[k]+=y[i]*y[i]; nc[k]++;
  }
  for(k=0; k<ncl; k++){
    y1[k]/=nc[k]; y2[k]=y2[k]/nc[k]-y1[k]*y1[k]; ll[k]=0;
  }
  for(i=0; i<nn; i++){
    k=cluster[i]; float z=y[i]-y1[k];
    ll[k]+=z*z;
  }
  double lik=0;
  for(k=0; k<ncl; k++){
    lik+=ll[k]/y2[k]+nc[k]*log(y2[k]);
  }
  lik *= (-0.5);
  lik -= 0.5*nn*log(6.283);
  free(y1); free(y2); free(ll); free(nc);
  return(lik);
}

struct peak *Center_peaks(int N_peak, struct peak *peaks,
			  float *y_scr, long *x_coord, long nn)
{
  struct peak *cpeaks=malloc(N_peak*sizeof(struct peak));
  struct peak *peak=peaks, *cpeak=cpeaks;
  int i, B=BOX1/2;
  for(i=0; i<N_peak; i++){
    // Find maximum
    int jmax=peak->ifrag, j=jmax+1; float y_max=y_scr[jmax];
    while(1){
      if(y_scr[j]>y_max){y_max=y_scr[j]; jmax=j;} j++;
      if((j==nn)||(x_coord[j]>peak->end)||(chr[j]>peak->chr))break;
    }
    long x0=x_coord[jmax]+B; // Half the way
    int d=x0-peak->ini, d2=peak->end-x0; if(d2 < d)d=d2;
    if(d<B)d=B;
    /*
    // Find maximum distance at which score is > DAMP*max in both dir.
    float y_thr=y_max*DAMP;
    int j1; for(j1=jmax-1; j1>=0; j1--)if(y_scr[j1]<y_thr)break; j1++;
    int j2; for(j2=jmax+1; j2<nn; j2++)if(y_scr[j2]<y_thr)break; j2--;
    int jd=jmax-j1, d2=j2-jmax; if(d2 < jd)jd=d2;*/

    // Set peak
    //cpeak->ini=x_coord[jmax-jd];
    //cpeak->end=x_coord[jmax+jd]+BOX1-1;
    cpeak->chr=peak->chr;
    cpeak->ini=x0-d;
    cpeak->end=x0+d;
    cpeak->size=cpeak->end-cpeak->ini+1;
    cpeak->xo=(cpeak->ini+cpeak->end)/2;
    cpeak->x1=cpeak->xo-PEAKSIZE;
    cpeak->x2=cpeak->xo+PEAKSIZE;
    //cpeak->x1=cpeak->ini;
    //cpeak->x2=cpeak->end;
    cpeak->y=y_max;
    cpeak->next=cpeak+1;
    if(0){
      printf("peak: %d %ld %ld %d %.1f  ",
	     peak->chr, peak->ini, peak->end, peak->size, peak->y); 
      printf("cpeak: %d %ld %ld %d %.1f",
	     cpeak->chr, cpeak->ini, cpeak->end, cpeak->size, cpeak->y);
      printf("   %d %d\n", peak->size-cpeak->size, peak->xo-cpeak->xo);
    } 
    peak=peak->next;
    cpeak=cpeak->next;
  }
  (cpeak-1)->next=NULL;
  return(cpeaks);
}

struct peak *Clusters2Frags(int *N_peak, int *cluster, float *y_scr, long nn)
{
  Cluster_size(numclus, numdom, cluster, chr, nn, 2);
  *N_peak=numdom[1];
  int N_frag=numclus[1];
  printf("%d peaks\n", *N_peak);
  struct peak *peaks=malloc(*N_peak*sizeof(struct peak));
  struct peak *peak=peaks;
  int i, m=0, k, k1=-1, n=0, ch1=-1, nch=-1; double y;
  int NCH=1000, numch[NCH]; for(i=0; i<NCH; i++)numch[i]=0;
  for(i=0; i<nn; i++){
    if(chr[i]!=ch1){ // New chromosome starts
      ch1=chr[i]; nch++;
      if(k1==1){
	n+=Close_peak(&peak, x_coord[i-1]+BOX1, y/m); k1=-1;
      }
    }
    k=cluster[i];
    if(k==1){
      if(k1!=1){ // k=1, k1=0 : start peak
	peak->ifrag=i;
	peak->ini=x_coord[i];
	peak->chr=chr[i];
	y=0; m=0;
	numch[nch]++;
      }
      y+=y_scr[i]; m++;
    }else if(k1==1){ // k=0, k1=1: close peak
      n+=Close_peak(&peak, x_coord[i], y/m);
    }
    k1=k;
  }
  if(k1==1)n+=Close_peak(&peak, x_coord[i-1]+BOX1, y/m);
  (peak-1)->next=NULL;
  printf("%d peaks\n", n);
  printf("%d chromosomes in Clusters2Frag\n", nch+1);
  for(i=0; i<=nch; i++)printf("%d peaks in chr %d\n", numch[i], i+1);
  return(peaks);
}
