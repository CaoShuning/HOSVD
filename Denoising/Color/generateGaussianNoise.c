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

#include "denoise_color.h"

int main (int argc, char **argv)
{
	double noisy_mse, noisy_psnr;
	FILE *fp;
	char opfname[1000];
	double **imt[3],**imt2[3];
	double SIGMA;
	int i,j,H,W,k;
	// time variables
	time_t start, end;
	double dtime;
	struct timeval tv;
	struct timezone tz;
	// GSL random number routines
       const gsl_rng_type * T;
       gsl_rng * r;	

	readPPM (argv[1],&H,&W,imt);
	SIGMA = atoi(argv[2]);

	for(i=0;i<3;i++)
	{
		imt2[i] = allocate_2d_double(H,W,'0');
	}

       /* create a generator chosen by the 
          environment variable GSL_RNG_TYPE */
       gsl_rng_env_setup();
       T = gsl_rng_default;
       r = gsl_rng_alloc (T);
	 gettimeofday(&tv,&tz);
	 printf ("\nTime  = %ld",tv.tv_sec); fflush(stdout);
	 gsl_rng_set (r, tv.tv_sec);

	// generating a noise instance [Gaussian Noise]
	for (i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			for(k=0;k<3;k++)
			{
	 			imt2[k][i][j] = imt[k][i][j] + gsl_ran_gaussian (r,SIGMA);
			}
		}
	}
	printf ("\nWriting noise file"); fflush (stdout);
	sprintf (opfname,"noisy_Sig%.2lf_%s",SIGMA,argv[1]);
	writePPM (opfname, imt2, H, W);

	noisy_mse = MSE(imt,imt2,H,W); noisy_psnr = PSNR(noisy_mse); 
	printf ("\nBefore denoising, MSE = %lf, PSNR = %lf",noisy_mse,noisy_psnr); fflush (stdout);
}

