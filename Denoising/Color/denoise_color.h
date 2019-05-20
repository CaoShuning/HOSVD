//
//Code for submission to TPAMI: Image Denoising using the Higher Order Singular Value Decomposition, (version %0.0.1):
//---------------------------------------------------
//Copyright (C) 2012 Ajit Rajwade, Anand Rangarajan and Arunava Banerjee
// Authors: Ajit Rajwade, Anand Rangarajan and Arunava Banerjee
// Date:    June 4th 2012
// 
// Contact Information:
//
//Ajit Rajwade:	avr@cise.ufl.edu
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

int SR;
int ps;
int K;
double THRESHOLDS[3];
double SIGMA;
double DISTANCE_THRESHOLD;

double **KLT_V, *KLT_mu;

//alloc.c
float **allocate_2d_float(int N,int M,char zero);
void free_2d_float(float **a,int N);
char **allocate_2d_char(int N,int M,char zero);
double **allocate_2d_double(int N,int M,char zero);
double ***allocate_3d_double(int N1, int N2, int N3, char zero);
int ***allocate_3d_int(int N1, int N2, int N3, char zero);
void free_2d_double(double **a,int N);
void free_3d_double (double ***a, int N1, int N2);
void free_3d_int (int ***a, int N1, int N2);
int **allocate_2d_int(int N,int M,char zero);
void free_2d_int(int **a,int N);
float *allocate_1d_float(int N,char zero);
double *allocate_1d_double(int N,char zero);
int *allocate_1d_int(int N,char zero);
unsigned char *allocate_1d_uchar(int N,char zero);

// pgmfile.c
void readPPM (char *fname, int *height, int *width,double **im[3]);
void writePPM (char *fname, double **im[3], int height, int width);

// qualitymetrics.c
double MSE (double **im1[3], double **im2[3], int r, int c);
double PSNR (double mse);


