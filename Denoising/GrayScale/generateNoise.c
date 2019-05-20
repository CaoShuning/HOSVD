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


#include "denoise.h"

int main (int argc, char **argv)
{
	FILE *fp;
	char opfname[1000];
	double **imt,x;
	double SIGMA;
	int i,j,H,W,seed;
	// time variables
	time_t start, end;
	double dtime;
	struct timeval tv;
	struct timezone tz;
	// GSL random number routines
       const gsl_rng_type * T;
       gsl_rng * r;	

	if (argc < 3)
	{
		printf ("\nUSAGE: ./generateNoise.o <clean file name> <sigma> <noise model>");
		printf ("\nFor noisemodel: (1) Gaussian, (2) Poisson, (3) Negative exponential. SIGMA is not relevant for Poisson model");
		exit (0);
	}

	imt = readPGM (argv[1],&H,&W);
	SIGMA = atof(argv[2]);
	NOISEMODEL = atoi(argv[3]);
	printf ("\nNoise model = %d",NOISEMODEL);

       /* create a generator chosen by the 
          environment variable GSL_RNG_TYPE */
       gsl_rng_env_setup();
       T = gsl_rng_default;
       r = gsl_rng_alloc (T);
	 gettimeofday(&tv,&tz);
	 printf ("\nTime  = %ld",tv.tv_sec); fflush(stdout);
	 gsl_rng_set (r, tv.tv_sec);

//	printf ("\nEnter seed: "); scanf ("%d",&seed);
//	gsl_rng_set (r, seed);

	// generating a noise instance [Gaussian Noise]
	for (i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			if (NOISEMODEL == GAUSSIAN_NOISE) imt[i][j] = imt[i][j] + gsl_ran_gaussian (r,SIGMA);
			else if (NOISEMODEL == POISSON_NOISE) imt[i][j] = (double)gsl_ran_poisson (r,imt[i][j]);
			else if (NOISEMODEL == NEGATIVE_EXPONENTIAL_NOISE) imt[i][j] = imt[i][j] + gsl_ran_exponential (r,SIGMA);
			else { printf ("\nIncorrect noise model value!"); exit (0); }
		}
	}
	printf ("\nWriting noise file"); fflush (stdout);
	if (NOISEMODEL == GAUSSIAN_NOISE) sprintf (opfname,"noisy_Sig%.2lf_%s",SIGMA,argv[1]);
	else if (NOISEMODEL == POISSON_NOISE) sprintf (opfname,"noisyPoisson_%s",argv[1]);
	else if (NOISEMODEL == NEGATIVE_EXPONENTIAL_NOISE) sprintf (opfname,"noisyExp_Mu%.2lf_%s",SIGMA,argv[1]);
	writePGM (opfname, imt, H, W);
}
