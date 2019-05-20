//
//Code for submission to TPAMI: Image Denoising using the Higher Order Singular Value Decomposition, (version %0.0.1):
//---------------------------------------------------
//Copyright (C) 2012 Ajit Rajwade, Anand Rangarajan and Arunava Banerjee
// Authors: Ajit Rajwade, Anand Rangarajan and Arunava Banerjee
// Date:    June 4th 2012
// 
// Contact Information:
//
//Ajit Rajwade: avr@cise.ufl.edu
//Anand Rangarajan: anand@cise.ufl.edu
//Arunava Banerjee: arunava@cise.ufl.edu
// Terms:         
// 
//The source code is provided under the
//terms of the GNU General Public License (version 3).
//

/*    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

#define ODDNUMBER 1
#define EVENNUMBER 2
#define EPSILON 1e-4
#define GAUSSIAN_NOISE 1
#define POISSON_NOISE 2
#define NEGATIVE_EXPONENTIAL_NOISE 3
#define NUM_NOISE_INSTANCES 5000
# define NUMBINS 100
#define FILESTRINGLENGTH 3000
#define VERBOSE 1
#define MY_INFINITY 10000000
#if !defined(_MACROS_)
#define _MACROS_
#define check_ptr(ptr,fn) \
        if((ptr)==NULL) { \
                printf("%10s: NULL pointer\n",(char *)(fn)); \
                perror((char *)(fn)); \
                exit(1); \
        }
#endif
#define CHOSEN_SIGNIFICANCE_LEVEL 0.3

char hosvd_filename[100];
int MAXK;
int SR;
int ps,halfps;
int K;
int NOISEMODEL;
double THRESHOLD, THRESHOLD2;
double SIGMA;
double DISTANCE_THRESHOLD;
double averageKval;

double **T1g, **Sg, **UTg,**S_noisy,**S_clean,**gPatch;
double **UDCT;

struct OP
{
	double ***UVs;
};
typedef struct OP OP;

struct ImageLL
{
	double **patch;
	struct ImageLL *next;
};
typedef struct ImageLL ImageLL;


struct thresholdFunction
{
	double *deltas;
	int *binindices;
	int num_nonzero;
	double minval;
	double maxval;
};
typedef struct thresholdFunction thresholdFunction;

//alloc.c
float **allocate_2d_float(int N,int M,char zero);
char **allocate_2d_char(int N,int M,char zero);
void free_2d_float(float **a,int N);
double **allocate_2d_double(int N,int M,char zero);
void free_2d_double(double **a,int N);
int **allocate_2d_int(int N,int M,char zero);
void free_2d_int(int **a,int N);
float *allocate_1d_float(int N,char zero);
double *allocate_1d_double(int N,char zero);
int *allocate_1d_int(int N,char zero);
unsigned char *allocate_1d_uchar(int N,char zero);
double ***allocate_3d_double(int N1, int N2, int N3, char zero);
int ***allocate_3d_int(int N1, int N2, int N3, char zero);

// pgmfile.c
double ** readPGM (char *fname, int *height, int *width);
void writePGM (char *fname, double **im, int height, int width);

// qualitymetrics.c
double MSE (double **im1, double **im2, int r, int c);
double PSNR (double mse);
double MSSIM (double **im1, double **im2, int r, int c);

// dctfilter2.c 
int getNearestNeighbor (double **patches, double *currpatch, int numpatches, int patchsize);
thresholdFunction* getThresholdFunctions (int numbins, char *tfilename, int numpatches);
//void dctFilter2 (double **im, double **im2, int ps, int H, int W, int numbins, char *patchfilename, char *tfilename);
void dctFilter2 (double **, double **, int, int, int, int, char *, char *, double**);
void getDCTFilter (double **A, int ps);

// hyptest.c
int mycompare (const void *a, const void *b);
double probks(double alam);
void ksGaussian (double *data, int n, double *d, double *prob, double sigma);
void kstwo (double *data1, int n1, double *data2, int n2, double *d, double *prob);

// matrices.c
void matrix_transpose (double **A, double **B, int r, int c);
void matrix_multiply (double **A, double **B, double **C, int r1, int c1, int c2);
void printmatrix (double **A, int r, int c);
void GSL2matrix (gsl_matrix *gF, double **F, int n1, int n2);
void matrix2GSL (gsl_matrix *gF, double **F, int n1, int n2);
void setmatrix_zero (double **A, int n1, int n2);
void kronecker_product (double **B, double **A1, double **A2, int n1, int n2);
void vector_outerproduct (double *a, double *b, double **A, int n);

// qualitymetrics.c
double MSE (double **im1, double **im2, int r, int c);
double PSNR (double mse);
double MSSIM (double **im1, double **im2, int r, int c);
void compute_residuals (double **res, double **imnoisy, double **imdenoised, int H, int W);

// hyptest.c
double KS_similar_patches (double ***residualpatch, int refindex, int numpatches);
double corrcoeff_Noiseness (double **res, int H, int W, int ps);
double corrcoeff_Noiseness_multiscale (double **res, int H, int W, int lowest_scale, int biggest_scale);

ImageLL **allocate_2d_ImageLL(int N,int M,char zero);
ImageLL ***allocate_3d_ImageLL(int N1, int N2, int N3, char zero);

// matrix3d
void getFolding1 (double ***A, double **B, int m1, int m2, int m3);
void getReverseFolding1 (double ***A, double **B, int m1, int m2, int m3);
void getFolding2 (double ***A, double **B, int m1, int m2, int m3);
void getFolding3 (double ***A, double **B, int m1, int m2, int m3);

