/***************************************************************
 *Computer Physics Communications 135 (2001) 40-57             *
 *Extension to the weighted histogram analysis method:         *
 *combining umbrella sampling with free energy calculations    *
 *Marc Souaille, Benot Roux                                    *
 *-------------------------------------------------------------*
 * CURRENT SUPPORTED FUNCTIONALITIES:
 *
 *1. 1D or 2D pmf along the restraint variables
 *2. 1D or 2D reweighted PMF along other variables
 *3. UMS, SITS-UMS 
 *4. WHAM for H-REX with local biasing (added on July-12-2014)
 *   metafile: file_i 0. lambda_i (stored in spring: 1D)
 *   metafile: file_i 0. 0. lambda_i 0. (stored in spring: 2D)
 *   file_i: step PMF-variable energy (1D)
 *   file_i: step PMF-variable1 PMF-variable2 energy (2D)
 *
 * *************************************************************/
 /*[rstr=H/C/R]:
 *H--Harmonic
 *C--Cosine
 *H--H-REX with local biasing
 */

//to compile: gcc -fopenmp -lm -O2 -o execfile (gcc version: 12.2.0,intel compiler is still not used)

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

//#include "wham.h"

#define MAXNITER 10000
#define RADIANS   6.28318530717959
#define PI 3.141592653589793
#define KBOLTZ 0.0019872064689 
#define LINESIZE 1024
#define TOL 1.e-4
#define EQU(a,b) (strncasecmp(a,b,strlen(b))==0)
#define COMMAND_LINE "wham metafile freefile [temp=k] [type=n] [nkm=n] [nkn=n] [ITS=0/1] [uDim=n] [rDim=m] [rstr=H/C/R] [histx=minx maxx num_binsx] [histy=miny maxy num_binsy] [histz=minz maxz num_binsz] [pmfx=minx maxx num_binsx] [pmfy=miny maxy num_binsy] [pmfz=minz maxz num_binsz] [PxNo|Px=xx] [PyNo|Py=xx] [PzNo|Pz=xx] [nproc=int]\n"
#define DETAILS "metafile: filename1 locx locy locz kfx kfy kfz [temp]\nfreefile: outputfile\ntemp: temperature of simulation\ntype: 0-normal UMS (default); 1-UMS-SITS; 2-reweight UMS; 3-reweight UMS-SITS; 10-REMD; 20-ITS-US for different nk at each window; 100-H-REX \nuDim: dim. of UMS simulation\nrDim: dim. of UMS reweight\nrstr: H-harmonic (0); C-Cosine (1); R-REX(2)\n"
#define DEFAULT "temp=300; type=0; uDim=1; rDim=1; rstr=0 (H);\n histx: -PI, PI, 60; histy: -PI, PI, 60; Px=2*PI; Py=2*PI\npmfx: 0., 0., 1.e200; pmfy: 0., 0., 1.e200\n"
                 
/*1. ITS-US with different nk values at each window
-----------------
[type=n] [nkm=n] [nkn=n] [ITS=0/1]:
type=20: invoke this application
nkm=m: m in SITS
nkn=n: n in SITS
ITS=0/1: 0 for ITS, 1 for SITS
-----------------format of files:
metafile: datafile nkfile loc kf
datafile: time-index RCs Ucomp0 Ucomp1 Ucomp2
nkfile: log(nk-values) betaK-values
-----------------related variables
wham->nk_num[]: number of nk used for each window
wham->betaK[][]: betaK values for each window
wham->nkValues[][]: nk values for each window
wham->Ucomp0[][]: energy component0 (in SITS)
wham->Ucomp1[][]: energy component1 (SITS-nonSITS)
wham->Ucomp2[][]: energy component2 (nonSITS)
wham->nk: can be replaced by calcType=20
wham->ITS: 0 for ITS; 1 for SITS
wham->m:
wham->n:
a0, b0, a1, b1
gftt, gft[], nkValues, e, betaK
wham->nkFileName[]: files storing nk values
wham->Usits-normal[][][]: store the calculated Usits - Unormal
*/

struct wham_info
{
    int nDim;
    int rstrType;
    int periodic[3];
    double period[3];
    double hist_min[3], hist_max[3], bin_width[3];
    int num_bins[3];
    int num_windows;
    int calcType;
    double T;
    /* number of separate biased trajectory windows*/
    double (*bias_locations)[3];
    /*array of locations of the bias for each window*/
    double (*spring_constants)[3];
    /*array of spring constants for the biases for each window*/
    double *F;
    /* array of free energy perturbations due to restraint*/
    double *kT;
    /* array of sampling temperatures, in kcal/mol */
    double *partition;
    double **U0, **Usits, ***X, ***rewX, **U;
    /*nt: num. of samples at each window*/
    int *nt;
    /*ebf: exp(beta*F(k)); ebf2: exp(-beta*F(k)); fact: nk*Exp(beta*F(k))*/
    /*ebw_ilk: store Exp(-beta*Wk(R_i, l)), i: window index, l: sample index 
     * in ith window, k: window index*/
    double *ebf, *ebfo, *ebf2, *fact, ***ebw;
    char metaFileName[LINESIZE], fepFileName[LINESIZE], (*dataFileName)[LINESIZE], (*nkFileName)[LINESIZE];
    FILE *metaFile, *fepFile;

    double ***prob, ***pmf, *probRed1, *probRed2, *probRed3, *pmfRed1, *pmfRed2, *pmfRed3;
    /*for calcType == 20*/
    double **nkValues, **betaK, m, n;
    int *nk_num;
    int nk; /*0: same set of nk values at each window of ITS-US
              1: different set of nk values for each window of ITS-US*/
    int ITS; /*0: ITS; 1: SITS*/
    double **Ucomp0, **Ucomp1, **Ucomp2, ***UsitsNormal;
    
};

struct pmf_rew_data
{
    int num_bins[3];
    double hist_min[3], hist_max[3], bin_width[3];
    int nDim;
    char metaFileName[LINESIZE], whamFileName[LINESIZE];
    FILE *metaFile, *whamFile;
    int num_windows;
    double (*bias_pos)[3], (*force_const)[3];
    double *part;
    double *pmf1, *prob1, **histogram1, **histogramb1;
    double **pmf2, **prob2, ***histogram2, ***histogramb2;
    double ***pmf3, ***prob3, ****histogram3, ****histogramb3;//20270723
};

double calc_bias_pot(struct wham_info *h, int index, double *coor, int BIAS);
int get_numwindows(FILE *file);
int is_metadata(char *line);
int CalcReweightPMF (struct wham_info *wham, struct pmf_rew_data *pmf);

int main(int argc, char *argv[])
{
    int i, j, k, kj, kjl, l, m, n, m1, m2, m3, n1, n2, n3, i1, i2, i3, j1, j2 , j3, k1, k2, k3, nthread;//20240806
    int current_window, vals, nline, conv;
    char *c, *line, filename[LINESIZE], filename2[LINESIZE];
    double loc, loc1, loc2, loc3, spring, spring1, spring2, spring3, correl_time, temp;//20240723 loc3 and spring3
    double beta, e_rstr, coor[3], step, e_sits;//20240723 coor[3]
    double ebfk, bottom, summ, delta;
    double ***histogram, denominator, fmin, fmin1, fmin2, fmin3, fk, nk;
    double a, b, nkValues, betaK, e, x, y, z;
    double gft[2000], gftt, a0, b0, a1, b1;
    int index, index1, index2, index3;
    struct wham_info *wham;
    struct pmf_rew_data *pmf;
    FILE *file, *file2;

    if(argc < 3) {
        printf("-------------------------********************************COMMAND_LINE\n");
        printf(COMMAND_LINE);
        printf("********************************-------------------------DETAILS\n");
        printf(DETAILS);
        printf("-------------------------********************************DEFAULTS\n");
        printf(DEFAULT);
        printf("#############################END@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        exit(0);
    }
    fflush(stdout);

    wham = (struct wham_info*) calloc(1, sizeof(struct wham_info));
    pmf = (struct pmf_rew_data*) calloc(1, sizeof(struct pmf_rew_data));

    /*set default values*/
    wham->T = 300.;
    wham->calcType = 0;
    wham->nDim = 1;
    wham->rstrType = 0;
    /*different nk values*/
    wham->nk = 0;
    wham->m = 2.;
    wham->n = 1.;
    wham->ITS = 0;

    for(i=0; i<3; i++) {//20240812
        wham->hist_min[i] = -PI;
        wham->hist_max[i] = PI;
        wham->num_bins[i] = 60;
        wham->periodic[i] = 1;
        wham->period[i] = 2*PI;

        pmf->hist_min[i] = -PI;
        pmf->hist_max[i] = PI;
        pmf->num_bins[i] = 60;
    }
    pmf->nDim = 1;

    /*print out command line*/
    printf("#");
    for (i=0; i<argc; i++) {
        printf(" %s", argv[i]);
    }
    printf("\n");

    /*read into info.*/
    strcpy (wham->metaFileName, argv[1]);
    strcpy (wham->fepFileName, argv[2]);
    argc = argc - 2;
    argv = argv + 2;

    if(argc >= 2) {
        for(i=1; i<argc; i=i+0) {
            if(EQU(argv[i], "temp")) {
                c = &(argv[i][5]);
                wham->T = atof(c);
                printf("temp: %lf\n", wham->T);
                argc--;
                argv++;
            } else if(EQU(argv[i], "type")) {
                c = &(argv[i][5]);
                wham->calcType = atoi(c);
                printf("calcType: %d\n", wham->calcType);
                argc--;
                argv++;
            } else if(EQU(argv[i], "nkm")) {
                c = &(argv[i][4]);
                wham->m = atof(c);
                printf("nkM: %f\n", wham->m);
                argc--;
                argv++;
            } else if(EQU(argv[i], "nkn")) {
                c = &(argv[i][4]);
                wham->n = atof(c);
                printf("nkN: %f\n", wham->n);
                argc--;
                argv++;
            } else if(EQU(argv[i], "ITS")) {
                c = &(argv[i][4]);
                wham->ITS = atoi(c);
                printf("ITS: %d\n", wham->ITS);
                argc--;
                argv++;
            } else if(EQU(argv[i], "uDim")) {
                c = &(argv[i][5]);
                wham->nDim = atoi(c);
                printf("uDim: %d\n", wham->nDim);
                argc--;
                argv++;
            } else if(EQU(argv[i], "rDim")) {
                c = &(argv[i][5]);
                pmf->nDim = atoi(c);
                printf("rDim: %d\n", pmf->nDim);
                argc--;
                argv++;
            } else if(EQU(argv[i], "rstr")) {
                c = &(argv[i][5]);
                wham->rstrType = atoi(c);
                printf("rstrType: %d\n", wham->rstrType);
                argc--;
                argv++;
            } else if(EQU(argv[i], "px")) {
                if(EQU(argv[i], "pxno")) {
                    wham->periodic[0] = 0;
                } else {
                    c = &(argv[i][3]);
                    wham->period[0] = atof(c);
                }
                printf("periodic: %d %lf\n", wham->periodic[0], wham->period[0]);
                argc--;
                argv++;
            } else if(EQU(argv[i], "py")) {
                if(EQU(argv[i], "pyno")) {
                    wham->periodic[1] = 0;
                } else {
                    c = &(argv[i][3]);
                    wham->period[1] = atof(c);
                }
                printf("periodic: %d %lf\n", wham->periodic[1], wham->period[1]);
                argc--;
                argv++;
            } else if (EQU(argv[i], "pz")) {//20240812
                if(EQU(argv[i], "pzno")) {
                    wham->periodic[2] = 0;
                } else {
                    c = &(argv[i][3]);
                    wham->period[2] = atof(c);
                }
                printf("periodic: %d %lf\n", wham->periodic[2], wham->period[2]);
                argc--;
                argv++;
            } else if(EQU(argv[i], "hist")) {
                if(EQU(argv[i], "histx")) {
                    c = &(argv[i][6]);
                    wham->hist_min[0] = atof(c);
                    wham->hist_max[0] = atof(argv[i+1]);
                    wham->num_bins[0] = atoi(argv[i+2]);
                    printf("histx: %lf %lf %d\n", wham->hist_min[0], wham->hist_max[0],
                            wham->num_bins[0]);
                } else if(EQU(argv[i], "histy")) {
                    c = &(argv[i][6]);
                    wham->hist_min[1] = atof(c);
                    wham->hist_max[1] = atof(argv[i+1]);
                    wham->num_bins[1] = atoi(argv[i+2]);
                    printf("histy: %lf %lf %d\n", wham->hist_min[1], wham->hist_max[1],
                            wham->num_bins[1]);
                } else if(EQU(argv[i], "histz")) {//20240723
                    c = &(argv[i][6]);
                    wham->hist_min[2] = atof(c);
                    wham->hist_max[2] = atof(argv[i+1]);
                    wham->num_bins[2] = atoi(argv[i+2]);
                    printf("histz: %lf %lf %d\n", wham->hist_min[2], wham->hist_max[2],
                            wham->num_bins[2]);
                } else {
                    c = &(argv[i][5]); /*hist=min max num_bin*/
                    wham->hist_min[0] = atof(c);
                    wham->hist_max[0] = atof(argv[i+1]);
                    wham->num_bins[0] = atoi(argv[i+2]);
                    printf("histx: %lf %lf %d\n", wham->hist_min[0], wham->hist_max[0],
                            wham->num_bins[0]);
                }

                argc = argc - 3;
                argv = argv + 3;
            } else if(EQU(argv[i], "pmf")) {
                if(EQU(argv[i], "pmfx")) {
                    c = &(argv[i][5]);
                    pmf->hist_min[0] = atof(c);
                    pmf->hist_max[0] = atof(argv[i+1]);
                    pmf->num_bins[0] = atoi(argv[i+2]);
                    printf("pmfhistx: %lf %lf %d\n", pmf->hist_min[0], pmf->hist_max[0],
                            pmf->num_bins[0]);
                } else if(EQU(argv[i], "pmfy")) {
                    c = &(argv[i][5]);
                    pmf->hist_min[1] = atof(c);
                    pmf->hist_max[1] = atof(argv[i+1]);
                    pmf->num_bins[1] = atoi(argv[i+2]);
                    printf("pmfhisty: %lf %lf %d\n", pmf->hist_min[1], pmf->hist_max[1],
                            pmf->num_bins[1]);
                } else if(EQU(argv[i], "pmfz")) {//20240723
                    c = &(argv[i][5]);
                    pmf->hist_min[2] = atof(c);
                    pmf->hist_max[2] = atof(argv[i+1]);
                    pmf->num_bins[2] = atoi(argv[i+2]);
                    printf("pmfhistz: %lf %lf %d\n", pmf->hist_min[2], pmf->hist_max[2],
                            pmf->num_bins[2]);
                } else {
                    c = &(argv[i][4]); /*pmf=minx maxx num_binsx*/
                    pmf->hist_min[0] = atof(c);
                    pmf->hist_max[0] = atof(argv[i+1]);
                    pmf->num_bins[0] = atoi(argv[i+2]);
                    printf("pmfhistx: %lf %lf %d\n", pmf->hist_min[0], pmf->hist_max[0],
                            pmf->num_bins[0]);
                }
                argc = argc - 3;
                argv = argv + 3;
            } else if (EQU(argv[i],"nproc")){
                c = &(argv[i][6]);
                nthread = atoi(c);
                printf("threads used: %d\n",nthread);
                argc--;
                argv++; 
            }else {
                printf("unknown input keyword %s\n", argv[i]);
                exit(0);
            }
        }
    }
    //printf("%d",nthread);
    /*calculate bin_width*/
    for(i=0; i<wham->nDim; i++) {
        wham->bin_width[i] = (wham->hist_max[i] - wham->hist_min[i])/(double)wham->num_bins[i];
    }

    for(i=0; i<pmf->nDim; i++) {
        pmf->bin_width[i] = (pmf->hist_max[i] - pmf->hist_min[i])/(double)pmf->num_bins[i];
    }

    /*parse the read info*/
    if(wham->nDim == 1) {
        wham->hist_min[1] = 0.;
        wham->hist_max[1] = 0.;
        wham->bin_width[1] = 1.e200;
        wham->num_bins[1] = 1;
        wham->periodic[1] = 0;
    }

    for(i=0; i<wham->nDim; i++) {
        if(wham->periodic[i] == 0) wham->period[i] = 1.e200;
    }

    if(pmf->nDim == 1) {
        pmf->hist_min[1] = 0.;
        pmf->hist_max[1] = 0.;
        pmf->bin_width[1] = 1.e200;
        pmf->num_bins[1] = 1;
    }

    if(wham->calcType < 2 || wham->calcType == 10 || wham->calcType == 20 ||
       wham->calcType == 100) { /*no reweight is required*/
        pmf->nDim = 0;
        pmf->hist_min[0] = 0.;
        pmf->hist_max[0] = 0.;
        pmf->bin_width[0] = 1.e200;
        pmf->num_bins[0] = 1;
    }

    /*for ITS/SITS-US*/
    if(wham->calcType == 20) { /*different set of nk values in ITS-US*/
        wham->nk = 1;
        if (wham->ITS == 0) { /* for ITS*/
            a0 = 0.; b0 = 0.;
            a1 = 1.; b1 = 1.;
        } else { /*for SITS*/
            a0 = 1./(wham->m); b0 = 1./(wham->n);
            a1 = (wham->m - 1.)/(wham->m); b1 = (wham->n - 1.)/(wham->n);
        }
    }

    /*print out the info*/
    if(wham->calcType == 0) {
        printf("WHAM for normal UMS simulations\n");
    } else if(wham->calcType == 1) {
        printf("WHAM for UMS-SITS simulations\n");
    } else if(wham->calcType == 2) {
        printf("WHAM for reweighting based on normal UMS simulations\n");
    } else if(wham->calcType == 3) {
        printf("WHAM for reweighting based on UMS-SITS simulations\n");
    } else if(wham->calcType == 10) {
        printf("WHAM for REMD simulations\n");
    } else if(wham->calcType == 20) {
        printf("WHAM for ITS-US simulations with diff. nk values\n");
        printf("ITS = 0 or ITS = 1 for SITS: %d\n", wham->ITS);
        printf("m = %lf n = %lf\n", wham->m, wham->n);
        printf("a0 = %lf b0 = %lf a1 = %lf b1 = %lf\n", a0, b0, a1, b1);
    } else if(wham->calcType == 100) {
	printf("WHAM for H-REX\n");
    } else {
        printf("WHAM: unknown type of task----%d----exit\n", wham->calcType);
        exit(0);
    }
    printf("metafile %20.20s \nfepfile %20.20s\n", wham->metaFileName, wham->fepFileName);
    printf("T: %lf WHAM-nDim: %d PMF-nDim: %d\n", wham->T, wham->nDim, pmf->nDim);
    if(wham->rstrType == 0) printf("Harmonic restraint used\n");
    else if(wham->rstrType == 1) printf("Cosine restraint used\n");
    else if(wham->rstrType == 2) printf("H-REX with local biasing\n");
    else {
        printf("Unknonw restraint type %d - exit\n", wham->rstrType);
        exit(0);
    }

    for(i=0; i<wham->nDim; i++) {
        printf("WHAM: %d hist_min %lf hist_max %lf bin_width %lf num_bins %d\n", i+1,
                wham->hist_min[i], wham->hist_max[i], wham->bin_width[i], wham->num_bins[i]);
    }

    for(i=0; i<pmf->nDim; i++) {
        printf("PMF: %d hist_min %lf hist_max %lf bin_width %lf num_bins %d\n", i+1,
                pmf->hist_min[i], pmf->hist_max[i], pmf->bin_width[i], pmf->num_bins[i]);
    }
    fflush(stdout);

    /*read the data from metafile*/
    wham->metaFile = fopen(wham->metaFileName, "r");
    if(wham->metaFile == (FILE*)NULL) {
        printf("Can not open the data file %30.20s - exit\n", wham->metaFileName);
        exit(0);
    }

    /*count number of windows in the data file*/
    wham->num_windows = get_numwindows(wham->metaFile);

    printf("WHAM: num. of windows used: %d \n", wham->num_windows);
    fflush(stdout);

    /*allocate space*/
    wham->bias_locations = calloc(wham->num_windows, sizeof(double[3]));//20240723
    wham->spring_constants = calloc(wham->num_windows, sizeof(double[3]));//20240723
    wham->F = calloc(wham->num_windows, sizeof(double));
    wham->kT = calloc(wham->num_windows, sizeof(double));
    wham->partition = calloc(wham->num_windows, sizeof(double));
    wham->nt = calloc(wham->num_windows, sizeof(double));
    wham->ebf = calloc(wham->num_windows, sizeof(double));
    wham->ebfo = calloc(wham->num_windows, sizeof(double));
    wham->ebf2 = calloc(wham->num_windows, sizeof(double));
    wham->fact = calloc(wham->num_windows, sizeof(double));
    wham->dataFileName = calloc(wham->num_windows, sizeof(char[LINESIZE]));
    wham->nkFileName = calloc(wham->num_windows, sizeof(char[LINESIZE]));
    wham->nk_num = calloc(wham->num_windows, sizeof(int));
    wham->prob = (double***)calloc(wham->num_bins[0], sizeof(double**));
    wham->pmf = (double***)calloc(wham->num_bins[0], sizeof(double**));
    for(i=0; i<wham->num_bins[0]; i++) {
        wham->pmf[i] = (double**)calloc(wham->num_bins[1], sizeof(double*));
        wham->prob[i] = (double**)calloc(wham->num_bins[1], sizeof(double*));
        for(j=0;j<wham->num_bins[1];j++) {
           wham->pmf[i][j] = (double*)calloc(wham->num_bins[2], sizeof(double));
           wham->prob[i][j] = (double*)calloc(wham->num_bins[2], sizeof(double));
        }
    }
        

    if(wham->nDim == 1) {
        for(i=0; i<wham->num_windows; i++) {
            wham->spring_constants[i][1] = 0.;
            wham->bias_locations[i][1] = 0.;
        }
    }

    /*read into info in metafile and count the number of samples in each window*/
    rewind(wham->metaFile);
    line = (char *) calloc(LINESIZE, sizeof(char));
    current_window = 0;
    while (fgets(line, LINESIZE, wham->metaFile)) {
        if (is_metadata(line)) {
            if(wham->nDim == 1 && wham->nk == 0) { /*filename, loc, kf*/
                vals = sscanf(line, "%s %lf %lf %lf %lf", filename, &loc, &spring,
                        &correl_time, &temp);
                if (vals >= 3) {
                    strcpy (wham->dataFileName[current_window], filename);
                    wham->bias_locations[current_window][0] = loc;
                    wham->spring_constants[current_window][0] = spring;
                    if (vals == 5) {
                        wham->kT[current_window] = temp;
                    }
                }
            } else if(wham->nDim == 1 && wham->nk == 1) { /*filename, nkfile, loc, kf*/
                vals = sscanf(line, "%s %s %lf %lf %lf %lf", 
                        filename, filename2, &loc, &spring,
                        &correl_time, &temp);
                if (vals >= 4) {
                    strcpy (wham->dataFileName[current_window], filename);
                    strcpy (wham->nkFileName[current_window], filename2);
                    wham->bias_locations[current_window][0] = loc;
                    wham->spring_constants[current_window][0] = spring;
                    if (vals == 6) {
                        wham->kT[current_window] = temp;
                    }
                }
            } else if(wham->nDim == 2 && wham->nk == 0) { 
                /*filename, locx, locy, kfx, kfy*/
                vals = sscanf(line, "%s %lf %lf %lf %lf %lf %lf", filename, 
                        &loc1, &loc2, &spring1, &spring2,
                        &correl_time, &temp);
                if (vals >= 5) {
                    strcpy (wham->dataFileName[current_window], filename);
                    wham->bias_locations[current_window][0] = loc1;
                    wham->bias_locations[current_window][1] = loc2;
                    wham->spring_constants[current_window][0] = spring1;
                    wham->spring_constants[current_window][1] = spring2;
                    if (vals == 7) {
                        wham->kT[current_window] = temp;
                    }
                }
            } else if(wham->nDim == 2 && wham->nk == 1) { 
                /*filename, nkfile, locx, locy, kfx, kfy*/
                vals = sscanf(line, "%s %s %lf %lf %lf %lf %lf %lf", 
                        filename, filename2,
                        &loc1, &loc2, &spring1, &spring2,
                        &correl_time, &temp);
                if (vals >= 6) {
                    strcpy (wham->dataFileName[current_window], filename);
                    strcpy (wham->nkFileName[current_window], filename2);
                    wham->bias_locations[current_window][0] = loc1;
                    wham->bias_locations[current_window][1] = loc2;
                    wham->spring_constants[current_window][0] = spring1;
                    wham->spring_constants[current_window][1] = spring2;
                    if (vals == 8) {
                        wham->kT[current_window] = temp;
                    }
                }    
            } else if(wham->nDim == 3 && wham->nk == 0) { //20240723
                /*filename, locx, locy, locz, kfx, kfy, kfz*/
                vals = sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf %lf", filename, 
                        &loc1, &loc2, &loc3, &spring1, &spring2, &spring3,
                        &correl_time, &temp);
                if (vals >= 7) {
                    strcpy (wham->dataFileName[current_window], filename);
                    wham->bias_locations[current_window][0] = loc1;
                    wham->bias_locations[current_window][1] = loc2;
                    wham->bias_locations[current_window][2] = loc3;
                    wham->spring_constants[current_window][0] = spring1;
                    wham->spring_constants[current_window][1] = spring2;
                    wham->spring_constants[current_window][2] = spring3;
                    if (vals == 9) {
                        wham->kT[current_window] = temp;
                    }
                }
            } else if(wham->nDim == 3 && wham->nk == 1) { //20240723
                /*filename, nkfile, locx, locy, locz, kfx, kfy, kfz*/
                vals = sscanf(line, "%s %s %lf %lf %lf %lf %lf %lf %lf %lf", 
                        filename, filename2,
                        &loc1, &loc2, &loc3, &spring1, &spring2, &spring3,
                        &correl_time, &temp);
                if (vals >= 8) {
                    strcpy (wham->dataFileName[current_window], filename);
                    strcpy (wham->nkFileName[current_window], filename2);
                    wham->bias_locations[current_window][0] = loc1;
                    wham->bias_locations[current_window][1] = loc2;
                    wham->bias_locations[current_window][2] = loc3;
                    wham->spring_constants[current_window][0] = spring1;
                    wham->spring_constants[current_window][1] = spring2;
                    wham->spring_constants[current_window][2] = spring3;
                    if (vals == 10) {
                        wham->kT[current_window] = temp;
                    }
                }
            }
            current_window++;
        }
    }

    if(current_window != wham->num_windows) {
        printf("Inconsistent no. of windows: %d %d\n", current_window, wham->num_windows);
        exit(0);
    }

    /*count samples in each window*/
    for(i=0; i<wham->num_windows; i++) {
        file = fopen(wham->dataFileName[i], "r");
        if(file == (FILE*)NULL) {
            printf("Can not open file: %30.25s - exit\n", wham->dataFileName[i]);
            exit(0);
        }

        nline = 0;
        while (fgets(line, LINESIZE, file)) {
            if(line[0] != '#') nline++;
        }

        wham->nt[i] = nline;
        wham->partition[i] = nline;
        fclose(file);
    }

    /*allocate space for ebw_ilk*/
    wham->ebw = (double ***) calloc(wham->num_windows, sizeof(double**));
    for(i=0; i<wham->num_windows; i++) {
        wham->ebw[i] = (double **) calloc(wham->nt[i], sizeof(double*));
        for(j=0; j<wham->nt[i]; j++) {
            wham->ebw[i][j] = (double *) calloc(wham->num_windows, sizeof(double));
        }
    }

    if (wham->calcType == 20 && wham->nk == 1) {
        wham->UsitsNormal = (double ***) calloc(wham->num_windows, sizeof(double**));
        for(i=0; i<wham->num_windows; i++) {
            wham->UsitsNormal[i] = (double **) calloc(wham->nt[i], sizeof(double*));
            for(j=0; j<wham->nt[i]; j++) {
                wham->UsitsNormal[i][j] = (double *) calloc(wham->num_windows, sizeof(double));
            }
        }
    }

    wham->X = (double ***) calloc(wham->num_windows, sizeof(double**));
    for(i=0; i<wham->num_windows; i++) {
        wham->X[i] = (double **) calloc(wham->nt[i], sizeof(double*));
        for(j=0; j<wham->nt[i]; j++) {
            wham->X[i][j] = (double *) calloc(wham->nDim, sizeof(double));
        }
    }
    if(pmf->nDim > 0) {
        wham->rewX = (double ***) calloc(wham->num_windows, sizeof(double**));
        for(i=0; i<wham->num_windows; i++) {
            wham->rewX[i] = (double **) calloc(wham->nt[i], sizeof(double*));
            for(j=0; j<wham->nt[i]; j++) {
                wham->rewX[i][j] = (double *) calloc(pmf->nDim, sizeof(double));
            }
        }
    }

    if(wham->calcType == 1 || wham->calcType == 3) {
        wham->U0 = (double **) calloc(wham->num_windows, sizeof(double*));
        wham->Usits = (double **) calloc(wham->num_windows, sizeof(double*));
        for(i=0; i<wham->num_windows; i++) {
            wham->U0[i] = (double *) calloc(wham->nt[i], sizeof(double));
            wham->Usits[i] = (double *) calloc(wham->nt[i], sizeof(double));
        }
    } else if (wham->calcType == 20) {
        wham->Ucomp0 = (double **) calloc(wham->num_windows, sizeof(double*));
        wham->Ucomp1 = (double **) calloc(wham->num_windows, sizeof(double*));
        wham->Ucomp2 = (double **) calloc(wham->num_windows, sizeof(double*));
        wham->nkValues = (double **) calloc(wham->num_windows, sizeof(double*));
        wham->betaK = (double **) calloc(wham->num_windows, sizeof(double*));
        for(i=0; i<wham->num_windows; i++) {
            wham->nkValues[i] = (double *) calloc(2000, sizeof(double));
            wham->betaK[i] = (double *) calloc(2000, sizeof(double));
            wham->Ucomp0[i] = (double *) calloc(wham->nt[i], sizeof(double));
            wham->Ucomp1[i] = (double *) calloc(wham->nt[i], sizeof(double));
            wham->Ucomp2[i] = (double *) calloc(wham->nt[i], sizeof(double));
        }
    } else if (wham->calcType == 10 || wham->calcType == 100) {
        wham->U = (double **) calloc(wham->num_windows, sizeof(double*));
        for(i=0; i<wham->num_windows; i++) {
            wham->U[i] = (double *) calloc(wham->nt[i], sizeof(double));
        }
    }
    
    /*read into nk values*/
    if(wham->nk == 1) {
        for(i=0; i<wham->num_windows; i++) {
            file = fopen(wham->nkFileName[i], "r");
            if(file == (FILE*)NULL) {
                printf("Can not open file: %30.25s - exit\n", wham->dataFileName[i]); exit(0);
            }

            nline = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                   vals = sscanf(line, "%lf %lf", &nkValues, &betaK);
                   wham->nkValues[i][nline] = nkValues;
                   wham->betaK[i][nline] = betaK;
                   nline++;
                }
            }
            wham->nk_num[i] = nline;
            fclose(file);
        } 
        /*code debug 1*/
        /*
        for(i=0; i<wham->num_windows; i++) {
            for (j=0; j<wham->nk_num[i]; j++) {
                   printf("%d nk_num %d nk_value %lf beta_value %lf\n", i, wham->nk_num[i], 
                   wham->nkValues[i][j], wham->betaK[i][j]);
            }
        }
        */
        /*dend debug 1*/
    }
    

    /*compute ebw*/
    beta = 1./(wham->T*KBOLTZ);
    if(wham->kT[0] > 1.e-3) {
        for(i=0; i<wham->num_windows; i++) {
            wham->kT[i] = 1./(wham->kT[i]*KBOLTZ);
            printf("REMD temp: %d %lf\n", i, wham->kT[i]);
        }
    }
    fflush(stdout);

    for(i=0; i<wham->num_windows; i++) {
        file = fopen(wham->dataFileName[i], "r");
        if(file == (FILE*)NULL) {
            printf("Can not open file: %30.25s - exit\n", wham->dataFileName[i]);
            exit(0);
        }

        if(wham->nDim == 1 && wham->calcType == 0 && wham->nk == 0) { /*step data*/

            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf", &step, coor+0);
                    wham->X[i][l][0] = coor[0];
                    for(k=0; k<wham->num_windows; k++) {
                        e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                        wham->ebw[i][l][k] = exp(-beta*e_rstr);
                    }
                    l++;
                }
            }

        } else if(wham->nDim == 1 && wham->calcType == 1 && wham->nk == 0) { /*step data U0 Usits*/

            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf %lf", &step, coor+0, 
                            wham->U0[i]+l, wham->Usits[i]+l);
                    wham->X[i][l][0] = coor[0];
                    e_sits = wham->Usits[i][l] - wham->U0[i][l];
                    for(k=0; k<wham->num_windows; k++) {
                        e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                        wham->ebw[i][l][k] = exp(-beta*(e_rstr + e_sits));
                    }
                    l++;
                }
            }

        } else if(wham->nDim == 1 && wham->calcType == 20 && wham->nk == 1) { /*step data U0 U1 U2*/

            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf %lf %lf", &step, coor+0, 
                            wham->Ucomp0[i]+l, wham->Ucomp1[i]+l, wham->Ucomp2[i]+l);
                    //wham->U0[i][l] = wham->Ucomp0[i][l] +  wham->Ucomp1[i][l] +  wham->Ucomp2[i][l];
                    wham->X[i][l][0] = coor[0];
                    //e_sits = wham->Usits[i][l] - wham->U0[i][l];
                    for(k=0; k<wham->num_windows; k++) {
                        /*calculate e_sits*/
                        /*e = e0 + a1*e1 + b1*e2;*/
                        e = wham->Ucomp0[i][l] + a1*wham->Ucomp1[i][l] + b1*wham->Ucomp2[i][l];
                        for(k1=0; k1<wham->nk_num[k]; k1++) {
                           gft[k1] = -wham->betaK[k][k1]*e + (wham->nkValues[k][k1]);
                        }
                        gftt = gft[0];
                        for(k1=1; k1<wham->nk_num[k]; k1++) {
                            if(gftt > gft[k1]) gftt = gftt + log(1. + exp(gft[k1] - gftt));
                            else gftt = gft[k1] + log(1. + exp(gftt - gft[k1]));
                        }
                        /*e_sits = a0*wham->Ucomp1[i][l] + b0*wham->Ucomp2[i][l];*/
                        e_sits = (-1./beta*gftt);
                        /*code debug 2*/
                        /*
                        printf("ilk: %d %d %d U_sits %f U_norm %f diff %f\n", i,
                        l, k, e_sits, e, e_sits-e);
                        */
                        /*end debug 2*/
                        e_sits -= e;

                        e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                        wham->UsitsNormal[i][l][k] = e_sits;
                        /*code debug 3*/
                         //printf("ilk: %d %d %d diff %lf e_rstr %lf \n", i, l, k, e_sits, e_rstr);
                        /*end debug 3*/

                        wham->ebw[i][l][k] = exp(-beta*(e_rstr + e_sits));
                    }
                    l++;
                }
            }

        } else if(wham->nDim == 1 && wham->calcType == 10 && wham->nk == 0) { /*step data U*/

            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf", &step, coor+0, wham->U[i]+l);
                    wham->X[i][l][0] = coor[0];
                    for(k=0; k<wham->num_windows; k++) {
                        a = wham->kT[k];
                        b = wham->U[i][l];
                        wham->ebw[i][l][k] = exp(-a*b);
                    }
                    l++;
                }
            }

        } else if(wham->nDim == 2 && wham->calcType == 10 && wham->nk == 0) { /*step x1 x2 U*/

            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf %lf", &step, coor+0, coor+1, wham->U[i]+l);
                    wham->X[i][l][0] = coor[0];
                    wham->X[i][l][1] = coor[1];
                    for(k=0; k<wham->num_windows; k++) {
                        a = wham->kT[k];
                        b = wham->U[i][l];
                        wham->ebw[i][l][k] = exp(-a*b);
                    }
                    l++;
                }
            }
        } else if(wham->nDim == 3 && wham->calcType == 10 && wham->nk == 0) { /*step x1 x2 x3 U*/

            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf %lf %lf", &step, coor+0, coor+1, coor+2, wham->U[i]+l);
                    wham->X[i][l][0] = coor[0];
                    wham->X[i][l][1] = coor[1];
                    wham->X[i][l][2] = coor[2];//20240723
                    for(k=0; k<wham->num_windows; k++) {
                        a = wham->kT[k];
                        b = wham->U[i][l];
                        wham->ebw[i][l][k] = exp(-a*b);
                    }
                    l++;
                }
            }

        } else if(wham->nDim == 1 && wham->calcType == 1) { /*step data U0 Usits*/
        } else if(wham->nDim == 1 && wham->calcType == 2) { /*step data U0 Usits*/
            printf("undone\n"); exit(0);
        } else if(wham->nDim == 1 && wham->calcType == 3) { /*step data U0 Usits*/
            printf("undone\n"); exit(0);
        } else if(wham->nDim == 1 && wham->calcType == 100 && wham->nk == 0) { 
		/*step data U0
		 *for H-REX
		 */

            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf", &step, coor+0, 
                            wham->U[i]+l);
                    wham->X[i][l][0] = coor[0];
		    coor[0] = wham->U[i][l];
                    for(k=0; k<wham->num_windows; k++) {
                        e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                        wham->ebw[i][l][k] = exp(-beta*e_rstr);
                    }
                    l++;
                }
            }
        } else if(wham->nDim == 2 && wham->calcType == 100 && wham->nk == 0) { 
		/*step data1 data2 U0
		 *for HREX
		 */

            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf %lf", &step, coor+0, coor+1, 
                            wham->U[i]+l);
                    wham->X[i][l][0] = coor[0];
                    wham->X[i][l][1] = coor[1];
		    coor[0] = wham->U[i][l];
                    for(k=0; k<wham->num_windows; k++) {
                        e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                        wham->ebw[i][l][k] = exp(-beta*e_rstr);
                    }
                    l++;
                }
            }
        } else if(wham->nDim == 3 && wham->calcType == 100 && wham->nk == 0) { //20240723
		/*step data1 data2 data3 U0
		 *for HREX
		 */

            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf %lf %lf", &step, coor+0, coor+1, coor+2,
                            wham->U[i]+l);
                    wham->X[i][l][0] = coor[0];
                    wham->X[i][l][1] = coor[1];
                    wham->X[i][l][2] = coor[2];
		    coor[0] = wham->U[i][l];
                    for(k=0; k<wham->num_windows; k++) {
                        e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                        wham->ebw[i][l][k] = exp(-beta*e_rstr);
                    }
                    l++;
                }
            }


        } else if(wham->nDim == 2 && wham->calcType == 0) { /*step data U0 Usits*/
            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf\n", &step, coor+0, coor+1);
                    wham->X[i][l][0] = coor[0];
                    wham->X[i][l][1] = coor[1];
                    for(k=0; k<wham->num_windows; k++) {
                        e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                        wham->ebw[i][l][k] = exp(-beta*e_rstr);
                    }
                    l++;
                }
            }
            //printf("undone\n"); exit(0);
        } else if(wham->nDim == 2 && wham->calcType == 1) { /*step data U0 Usits*/
            printf("undone\n"); exit(0);
        } else if(wham->nDim == 2 && wham->calcType == 2) { /*step data U0 Usits*/
            printf("undone\n"); exit(0);
        } else if(wham->nDim == 2 && wham->calcType == 3) { /*step data U0 Usits*/
            printf("undone\n"); exit(0);
        } else if(wham->nDim == 2 && wham->calcType == 20 && wham->nk == 1) { 
            /*step data1 data2 U0 U1 U2*/

            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf %lf %lf %lf", 
                            &step, coor+0, coor+1,
                            wham->Ucomp0[i]+l, wham->Ucomp1[i]+l, wham->Ucomp2[i]+l);
                    //wham->U0[i][l] = wham->Ucomp0[i][l] +  wham->Ucomp1[i][l] +  wham->Ucomp2[i][l];
                    wham->X[i][l][0] = coor[0];
                    wham->X[i][l][1] = coor[1];
                    //e_sits = wham->Usits[i][l] - wham->U0[i][l];
                    for(k=0; k<wham->num_windows; k++) {
                        /*calculate e_sits*/
                        /*e = e0 + a1*e1 + b1*e2;*/
                        e = wham->Ucomp0[i][l] + a1*wham->Ucomp1[i][l] + b1*wham->Ucomp2[i][l];
                        for(k1=0; k1<wham->nk_num[k]; k1++) {
                           gft[k1] = -wham->betaK[k][k1]*e + (wham->nkValues[k][k1]);
                        }
                        gftt = gft[0];
                        for(k1=1; k1<wham->nk_num[k]; k1++) {
                            if(gftt > gft[k1]) gftt = gftt + log(1. + exp(gft[k1] - gftt));
                            else gftt = gft[k1] + log(1. + exp(gftt - gft[k1]));
                        }
                        /*e_sits = a0*wham->Ucomp1[i][l] + b0*wham->Ucomp2[i][l];*/
                        e_sits = (-1./beta*gftt);
                        /*code debug 2*/
                        /*
                        printf("ilk: %d %d %d U_sits %f U_norm %f diff %f\n", i,
                        l, k, e_sits, e, e_sits-e);
                        */
                        /*end debug 2*/
                        e_sits -= e;

                        e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                        wham->UsitsNormal[i][l][k] = e_sits;
                        /*code debug 3*/
                         //printf("ilk: %d %d %d diff %lf e_rstr %lf \n", i, l, k, e_sits, e_rstr);
                        /*end debug 3*/

                        wham->ebw[i][l][k] = exp(-beta*(e_rstr + e_sits));
                    }
                    l++;
                }
            }

        } else if(wham->nDim == 3 && wham->calcType == 0) { /*step data U0 Usits*/
            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf %lf\n", &step, coor+0, coor+1, coor+2);
                    wham->X[i][l][0] = coor[0];
                    wham->X[i][l][1] = coor[1];
                    wham->X[i][l][2] = coor[2];
                    for(k=0; k<wham->num_windows; k++) {
                        e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                        wham->ebw[i][l][k] = exp(-beta*e_rstr);
                    }
                    l++;
                }
            }
            //printf("undone\n"); exit(0);
        } else if(wham->nDim == 3 && wham->calcType == 1) { /*step data U0 Usits*/
            printf("undone\n"); exit(0);
        } else if(wham->nDim == 3 && wham->calcType == 2) { /*step data U0 Usits*/
            printf("undone\n"); exit(0);
        } else if(wham->nDim == 3 && wham->calcType == 3) { /*step data U0 Usits*/
            printf("undone\n"); exit(0);

        
        } else if(wham->nDim == 3 && wham->calcType == 20 && wham->nk == 1) { //20240723
            /*step data1 data2 data3 U0 U1 U2*/

            l = 0;
            while (fgets(line, LINESIZE, file)) {
                if(line[0] != '#') {
                    vals = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", 
                            &step, coor+0, coor+1, coor+2,
                            wham->Ucomp0[i]+l, wham->Ucomp1[i]+l, wham->Ucomp2[i]+l);
                    //wham->U0[i][l] = wham->Ucomp0[i][l] +  wham->Ucomp1[i][l] +  wham->Ucomp2[i][l];
                    wham->X[i][l][0] = coor[0];
                    wham->X[i][l][1] = coor[1];
                    wham->X[i][l][2] = coor[2];//202407230
                    //e_sits = wham->Usits[i][l] - wham->U0[i][l];
                    for(k=0; k<wham->num_windows; k++) {
                        /*calculate e_sits*/
                        /*e = e0 + a1*e1 + b1*e2;*/
                        e = wham->Ucomp0[i][l] + a1*wham->Ucomp1[i][l] + b1*wham->Ucomp2[i][l];
                        for(k1=0; k1<wham->nk_num[k]; k1++) {
                           gft[k1] = -wham->betaK[k][k1]*e + (wham->nkValues[k][k1]);
                        }
                        gftt = gft[0];
                        for(k1=1; k1<wham->nk_num[k]; k1++) {
                            if(gftt > gft[k1]) gftt = gftt + log(1. + exp(gft[k1] - gftt));
                            else gftt = gft[k1] + log(1. + exp(gftt - gft[k1]));
                        }
                        /*e_sits = a0*wham->Ucomp1[i][l] + b0*wham->Ucomp2[i][l];*/
                        e_sits = (-1./beta*gftt);
                        /*code debug 2*/
                        /*
                        printf("ilk: %d %d %d U_sits %f U_norm %f diff %f\n", i,
                        l, k, e_sits, e, e_sits-e);
                        */
                        /*end debug 2*/
                        e_sits -= e;

                        e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                        wham->UsitsNormal[i][l][k] = e_sits;
                        /*code debug 3*/
                         //printf("ilk: %d %d %d diff %lf e_rstr %lf \n", i, l, k, e_sits, e_rstr);
                        /*end debug 3*/

                        wham->ebw[i][l][k] = exp(-beta*(e_rstr + e_sits));
                    }
                    l++;
                }
            }

        }
        if(l != wham->nt[i]) {
            printf("Inconsistent no. of lines in data file: %30.25s %d %d - exit\n", 
                    wham->dataFileName[i], l, wham->nt[i]);
        }
        fclose(file);
    }

    /*iteration for wham*/
    /*initialization of fk */
    for(i=0; i<wham->num_windows; i++) {
        wham->F[i] = 0.;
        wham->ebf[i] = exp(beta*wham->F[i]);
        wham->ebfo[i] = wham->ebf[i];
        wham->fact[i] = wham->nt[i]*wham->ebf[i];
    }

    for(n=0; n<=MAXNITER; n++) { /*start iteration*/

        for(k=0; k<wham->num_windows; k++) { /*loop over windows: exp(-beta*fk)*/ 
            ebfk = 0.;

            # pragma omp parallel for num_threads(nthread) \
                reduction(+: ebfk) schedule(guided)
            for(i=0; i<wham->num_windows; i++) { /*loop over windows: sum over windows*/
                for(l=0; l<wham->nt[i]; l++) { /*loop over time series*/
                    bottom = 0.;
                    for(j=0; j<wham->num_windows; j++) { /*most inner loop*/
                        bottom += wham->ebw[i][l][j]*wham->fact[j];
                    }
                    ebfk += wham->ebw[i][l][k]/bottom;
                }
            }
            wham->ebf2[k] = ebfk; /*exp(-beta*fk_new)*/
            wham->ebf[k] = 1./(wham->ebf[0]*ebfk); /*exp(beta*(fk_new - f0_new))*/
            wham->fact[k] = wham->nt[k]*wham->ebf[k]; /*nk*exp(beta*fk_new)*/
        }

        /*check the convergence*/
        conv = 1;
        summ = 0.;
        for(k=0; k<wham->num_windows; k++) {
            wham->ebf[k] = wham->ebf2[0]/wham->ebf2[k]; /*shifted, new exp(beta*fk)*/

            delta = fabs(1./beta*log(wham->ebf[k]/wham->ebfo[k]));

            if(delta >= TOL) conv = 0;
            wham->ebfo[k] = wham->ebf[k];
            summ += delta;
        }

        printf("### %d %lf\n", n, summ);
        if(n%100 == 0) {
            for(i=0; i<wham->num_windows; i++) {
                printf("$$$ %lf\n", 1./beta*log(wham->ebf[i]));
            }
        }
        if(conv > 0) break;
        fflush(stdout);
    }

    printf("#########converged F[k] values: ###########\n");
    for(i=0; i<wham->num_windows; i++) {
        printf("%lf %lf %d\n", wham->F[i], 1./beta*log(wham->ebf[i]), wham->nt[i]);
        wham->F[i] = 1./beta*log(wham->ebf[i]);
    }
    fflush(stdout);


    /*calculate pmf*/

    histogram = (double***) calloc(wham->num_bins[0], sizeof(double**));//20240806
    for (i=0; i<wham->num_bins[0]; i++) {
        histogram[i]=(double**) calloc(wham->num_bins[1],sizeof(double*));
        for (int j=0;j<wham->num_bins[1];j++){
            histogram[i][j]=(double*) calloc(wham->num_bins[2],sizeof(double));
        }
        //histogram[i] = (double*) calloc(wham->num_bins[1], sizeof(double));
    }

    for (i=0; i<wham->num_windows; i++) {//20240806
        for(k=0; k<wham->num_bins[0]; k++)
           for(kj=0; kj<wham->num_bins[1]; kj++)
                for (int kjl=0;kjl<wham->num_bins[2];kjl++)
                    histogram[k][kj][kjl] = 0.;

        /*debug*/
        //sprintf(filename2, "%s%d.dat", "t", i);
        //file2 = fopen(filename2, "w");
        /*end debug*/


        n1 = 0;
        for(j=0; j<wham->nt[i]; j++) {
            k1 = 0; k2 = 0; k3 = 0;
            if(wham->X[i][j][0] >= wham->hist_min[0] && 
               wham->X[i][j][0] < wham->hist_max[0]) {
                index1 = 
                (int) ((wham->X[i][j][0] - wham->hist_min[0])/wham->bin_width[0]);
                k1 = 1;
                /*debug*/
                //fprintf(file2, "%d %lf\n", index, wham->X[i][j][0]);
                /*end debug*/
                n1++;
            }

            if(wham->nDim == 1) {
                index2 = 0;
                index3 = 0;//20240806
                k2 = 1;
                k3 = 1;//20240806
            } else if(wham->X[i][j][1] >= wham->hist_min[1] && 
                      wham->X[i][j][1] < wham->hist_max[1] && wham->nDim == 2) {
                index2 = 
                (int) ((wham->X[i][j][1] - wham->hist_min[1])/wham->bin_width[1]);
                k2 = 1;
                index3 = 0;
                k3 = 1;
            } else if(wham->X[i][j][1] >= wham->hist_min[1] && 
                      wham->X[i][j][1] < wham->hist_max[1] && wham->X[i][j][2] >= wham->hist_min[2] && 
                      wham->X[i][j][2] < wham->hist_max[2] && wham->nDim == 3) {
                        index2 = (int) ((wham->X[i][j][1] - wham->hist_min[1])/wham->bin_width[1]);
                        index3 = 
                        (int) ((wham->X[i][j][2] - wham->hist_min[2])/wham->bin_width[2]);
                        k3 = 1;
                        k2 = 1;
            }

            denominator = 0.;
            coor[0] = wham->X[i][j][0];
            coor[1] = wham->X[i][j][1];
            coor[2] = wham->X[i][j][2];//20240806
            if(wham->calcType == 0) {
                for (k=0; k<wham->num_windows; k++) {
                    e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                    fk = wham->F[k];
                    nk = wham->nt[k];
                    denominator += nk*exp(-beta*(e_rstr - fk));
                }
            } else if (wham->calcType == 20 && wham->nk == 1) {
                for (k=0; k<wham->num_windows; k++) {
                    e_sits = wham->UsitsNormal[i][j][k];
                    e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                    fk = wham->F[k];
                    nk = wham->nt[k];
                    denominator += nk*exp(-beta*(e_rstr- fk + e_sits));
                }
            } else if (wham->calcType == 100 && wham->nk == 0) {
		coor[0] = wham->U[i][j];
                for (k=0; k<wham->num_windows; k++) {
                    e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                    fk = wham->F[k];
                    nk = wham->nt[k];
                    denominator += nk*exp(-beta*(e_rstr- fk));
                }
            } else {
                for (k=0; k<wham->num_windows; k++) {
                    e_sits = wham->Usits[i][j] - wham->U0[i][j];
                    e_rstr = calc_bias_pot(wham, k, coor, wham->rstrType);
                    fk = wham->F[k];
                    nk = wham->nt[k];
                    denominator += nk*exp(-beta*(e_rstr- fk + e_sits));
                }
            }

            if(k1 == 1 && k2 == 1 && k3 == 1) {//20240806
                histogram[index1][index2][index3] += 1./denominator;
            }
        }
        for(k=0; k<wham->num_bins[0]; k++) {
            for(kj=0; kj<wham->num_bins[1]; kj++) {//20240806
                for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                    wham->prob[k][kj][kjl] += histogram[k][kj][kjl];//20240806
                }
            }
        }
        /*debug*/
        //fclose(file2);
        /*end debug*/
    }
    /*normalization*/
    summ = 0.;
    for(k=0; k<wham->num_bins[0]; k++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                summ += wham->prob[k][kj][kjl];
            }//20240806
        }
    }

    for(k=0; k<wham->num_bins[0]; k++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                wham->prob[k][kj][kjl] /= summ;
            }//20240806
        }
    }

    for(k=0; k<wham->num_bins[0]; k++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                printf("prob: %26.18f \n", wham->prob[k][kj][kjl]);
            }//20240806
        }
    }

    fmin = 1.e200;
    for(i=0; i<wham->num_bins[0]; i++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                if(wham->prob[i][kj][kjl] > 1.e-200) {
                    wham->pmf[i][kj][kjl] = -1./beta*log(wham->prob[i][kj][kjl]);//20240806
                } else {
                    wham->pmf[i][kj][kjl] = 99999999.;
                }
                if(fmin > wham->pmf[i][kj][kjl]) fmin = wham->pmf[i][kj][kjl];
            }
        }
    }

    for(i=0; i<wham->num_bins[0]; i++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){//20240806
                wham->pmf[i][kj][kjl] -= fmin;
            }
        }
    }

    for(i=0; i<wham->num_bins[0]; i++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                x = (i + 0.5)*wham->bin_width[0] + wham->hist_min[0];
                y = (kj + 0.5)*wham->bin_width[1] + wham->hist_min[1];
                z = (kjl+ 0.5)*wham->bin_width[2] + wham->hist_min[2];//20240806
                if(wham->nDim == 1) y = 0.;
                if(wham->nDim == 1 || wham->nDim == 2){z=0.;}//20240806
                printf("PMFF: %lf %lf %lf %26.18f %26.18f\n", x, y, z, wham->pmf[i][kj][kjl], wham->prob[i][kj][kjl]);//20240806
            }
        }
    }

    /*reduced 1d pmf*/
    /*if(wham->nDim == 2) {
        wham->probRed1 = (double*) calloc(wham->num_bins[0], sizeof(double));
        wham->probRed2 = (double*) calloc(wham->num_bins[1], sizeof(double));
        wham->pmfRed1 = (double*) calloc(wham->num_bins[0], sizeof(double));
        wham->pmfRed2 = (double*) calloc(wham->num_bins[1], sizeof(double));
        
        fmin1 = 1.e200; fmin2 = 1.e200;
        for(i=0; i<wham->num_bins[0]; i++) {
            for(kj=0; kj<wham->num_bins[1]; kj++) {
                wham->probRed1[i] += wham->prob[i][kj];
            }

            if(wham->probRed1[i] > 1.e-200) {
                wham->pmfRed1[i] = -1./beta*log(wham->probRed1[i]);
            } else {
                wham->pmfRed1[i] = 99999999.;
            }
            if(fmin1 > wham->pmfRed1[i]) fmin1 = wham->pmfRed1[i];
        }

        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for(i=0; i<wham->num_bins[0]; i++) {
                wham->probRed2[kj] += wham->prob[i][kj];
            }

            if(wham->probRed2[kj] > 1.e-200) {
                wham->pmfRed2[kj] = -1./beta*log(wham->probRed2[kj]);
            } else {
                wham->pmfRed2[kj] = 99999999.;
            }
            if(fmin2 > wham->pmfRed2[kj]) fmin2 = wham->pmfRed2[kj];
        }*/
        //20240806
        if(wham->nDim == 3) {
        wham->probRed1 = (double*) calloc(wham->num_bins[0], sizeof(double));
        wham->probRed2 = (double*) calloc(wham->num_bins[1], sizeof(double));
        wham->probRed3 = (double*) calloc(wham->num_bins[2], sizeof(double));
        wham->pmfRed1 = (double*) calloc(wham->num_bins[0], sizeof(double));
        wham->pmfRed2 = (double*) calloc(wham->num_bins[1], sizeof(double));
        wham->pmfRed3 = (double*) calloc(wham->num_bins[2], sizeof(double));
        
        fmin1 = 1.e200; fmin2 = 1.e200, fmin3=1.e200;
        for(i=0; i<wham->num_bins[0]; i++) {
            for(kj=0; kj<wham->num_bins[1]; kj++) {
                for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                    wham->probRed1[i] += wham->prob[i][kj][kjl];//20240806
                }
            }

            if(wham->probRed1[i] > 1.e-200) {
                wham->pmfRed1[i] = -1./beta*log(wham->probRed1[i]);
            } else {
                wham->pmfRed1[i] = 99999999.;
            }
            if(fmin1 > wham->pmfRed1[i]) fmin1 = wham->pmfRed1[i];
        }

        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for(i=0; i<wham->num_bins[0]; i++) {
                for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                    wham->probRed2[kj] += wham->prob[i][kj][kjl];//20240806
                }
            }

            if(wham->probRed2[kj] > 1.e-200) {
                wham->pmfRed2[kj] = -1./beta*log(wham->probRed2[kj]);
            } else {
                wham->pmfRed2[kj] = 99999999.;
            }
            if(fmin2 > wham->pmfRed2[kj]) fmin2 = wham->pmfRed2[kj];
        }

        for(kjl=0; kjl<wham->num_bins[2]; kjl++) {//20240806
            for(kj=0; kj<wham->num_bins[1]; kj++) {
                for (i=0; i<wham->num_bins[0]; i++){
                    wham->probRed2[kjl] += wham->prob[i][kj][kjl];//20240806
                }
            }

            if(wham->probRed3[kjl] > 1.e-200) {
                wham->pmfRed3[kjl] = -1./beta*log(wham->probRed3[kjl]);
            } else {
                wham->pmfRed3[kjl] = 99999999.;
            }
            if(fmin3 > wham->pmfRed3[kjl]) fmin2 = wham->pmfRed3[kjl];
        }

        for(i=0; i<wham->num_bins[0]; i++) {
            printf("Reduced1: %lf  %lf\n", 
                   wham->pmfRed1[i] - fmin1, wham->probRed1[i]);
        }
        for(i=0; i<wham->num_bins[1]; i++) {
            printf("Reduced2: %lf  %lf\n",
            wham->pmfRed2[i] - fmin2, wham->probRed2[i]);
        }
        for (i=0; i<wham->num_bins[2]; i++){
            printf("Reduced3: %lf  %lf\n",
            wham->pmfRed3[i] - fmin3, wham->probRed3[i]);//20240806 un
        }
    }

    /*print out the results for reweighting use*/
    printf("#############################PMF#############################\n");

    /*1st line*/
    printf("dimension: %d\n", wham->nDim);

    /*2nd line*/
    printf("restraint: %d\n", wham->rstrType);
    /*3rd line*/
    if(wham->nDim == 1) {
        printf("periodic: %d %f\n", wham->periodic[0], wham->period[0]);
    } else if(wham->nDim == 2){
        printf("periodic: %d %d %f %f\n", 
        wham->periodic[0], wham->periodic[1], wham->period[0], wham->period[1]);
    } else {//20240806
        printf("perodic: %d %d %d %f %f %f\n",
        wham->periodic[0],wham->periodic[1], wham->periodic[2], wham->period[0], wham->period[1], wham->period[2]);
    }
    /*4th line*/
    if(wham->nDim == 1) {
        printf("bins: %f %f %d\n", wham->hist_min[0], wham->hist_max[0], wham->num_bins[0]);
    } else if(wham->nDim == 2) {
        printf("bins: %f %f %d %f %f %d\n", wham->hist_min[0], wham->hist_max[0], 
        wham->num_bins[0], wham->hist_min[1], wham->hist_max[1], wham->num_bins[1]);
    } else {
        printf("bins: %f %f %d %f %f %d %f %f %d\n", wham->hist_min[0], wham->hist_max[0], 
        wham->num_bins[0], wham->hist_min[1], wham->hist_max[1], 
        wham->num_bins[1], wham->hist_min[2], wham->hist_max[2], wham->num_bins[2]);//20240806
    }
    /*5th*/
    printf("windows: %d\n", wham->num_windows);
    /*6th line*/
    if(wham->nDim == 1) {
        printf("WIN_i, bias_loc, spring_const, F, partition, temp\n");
    } else if(wham->nDim == 2){
        printf("WIN_i, bias_locx, bias_locy, spring_constx, spring_consty, F, partition, temp\n");
    } else {//20240806
        printf("WIN_i, bias_locx, bias_locy, bias_locz, spring_constx, spring_consty, spring_constz, F, partition, temp\n");
    }
    
    if(wham->nDim == 1) {
        for (i=0; i<wham->num_windows; i++) {
            printf("%d %f %f %26.18f %f %f\n", i, wham->bias_locations[i][0], 
            wham->spring_constants[i][0], wham->F[i], wham->partition[i], wham->kT[i]);
        }
    } else if(wham->nDim == 2){
        for (i=0; i<wham->num_windows; i++) {
            printf("%d %f %f %f %f %26.18f %f %f\n", i, wham->bias_locations[i][0],
            wham->bias_locations[i][1], wham->spring_constants[i][0],
            wham->spring_constants[i][1], wham->F[i], wham->partition[i], wham->kT[i]);
        }
    } else {//20240806
        for (i=0; i<wham->num_windows; i++) {
            printf("%d %f %f %f %f %f %f %26.18f %f %f\n", i, wham->bias_locations[i][0],
            wham->bias_locations[i][1], wham->bias_locations[i][2], wham->spring_constants[i][0],
            wham->spring_constants[i][1], wham->spring_constants[i][2], wham->F[i], wham->partition[i], wham->kT[i]);
        }
    }

    return(1);
}


/****************************************************/
int get_numwindows(FILE *file)
{
    char *line;
    int  num_windows;

    num_windows = 0;
    line = malloc(sizeof(char)*LINESIZE);
    if (!line)
    {
        printf("couldn't allocate space for line\n");
        exit(-1);
    }

    /* make sure we're at the beginning of the file*/
    rewind(file);
    while (fgets(line,LINESIZE,file))
    {
        if (is_metadata(line)) num_windows++;
    }
    free(line);
    return num_windows;
}
 
  
int is_metadata(char *line)
{
    int i,length;
    int not_space;

    /* check if it's a comment*/
    if (line[0] == '#') return 0;
    /* verify it's not a blank line*/
    length = strlen(line);
    not_space = 0;
    i = 0;
    while (!not_space && i < length)
    {
        if (!isspace(line[i])) not_space = 1;
        i++;
    }
    return not_space;
}

double calc_bias_pot(struct wham_info *h, int index, double *coor, int BIAS)
{
    double springx, springy, springz;
    double dx,dy,dz;
    
    /*for 1d case, h->spring_constants[index][1] should be zero*/
    springx = h->spring_constants[index][0];
    springy = h->spring_constants[index][1];
    springz = h->spring_constants[index][2];
    dx = coor[0] - h->bias_locations[index][0];
    dy = coor[1] - h->bias_locations[index][1];
    dz = coor[2] - h->bias_locations[index][2];

    /*for non-periodic case, period[i] should be a huge number*/
    if (BIAS == 0) {
        dx = fabs(dx);
        if (dx > h->period[0]/2.)
        {
            dx -= h->period[0];
        }

        dy = fabs(dy);
        if (dy > h->period[1]/2.)
        {
            dy -= h->period[1];
        }

        dz = fabs(dz);//20240806
        if (dz > h->period[2]/2.)
        {
            dz -= h->period[2];
        }
        return 0.5*(springx*(dx*dx) + springy*(dy*dy) + springz*(dz*dz));
    } else if(BIAS == 1) {
        return ((springx*(1. - cos(dx)) + springy*(1. - cos(dy))));
    } else if(BIAS == 2) {
	    return (springx*dx);
    } else {
	    printf("unknown rstr type %d -- exit \n", BIAS);
	    exit (1);
    }

    return (1);
}



