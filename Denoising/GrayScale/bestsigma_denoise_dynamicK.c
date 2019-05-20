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

void copyimage (double **dest, double **src, int H, int W)
{
	int i,j;

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			dest[i][j] = src[i][j];
		}
	}

}


void getMinK (double *dist, int **indices, int **minK, int n)
{
	int i,j,minindex;
	double nextmin;

	for (i=0;i<K;i++)
	{
		nextmin = INFINITY;
		for(j=0;j<n;j++)
		{
			if (dist[j] < nextmin) {nextmin = dist[j]; minindex = j;}
		}

		minK[i][0] = indices[minindex][0];
		minK[i][1] = indices[minindex][1];
		dist[minindex] = INFINITY+10;
	}
}

void getNearestNeighbors (double *dist,int **indices, int **nearestNeigh, int n)
{
	int i,j=0;

	for (i=0;i<n;i++)
	{
		if (dist[i] < DISTANCE_THRESHOLD) 
		{
			nearestNeigh[j][0] = indices[i][0];
			nearestNeigh[j][1] = indices[i][1];
			j++;
		}
	}
}

void getNearestNeighborsProb (double *prob, double *dist,int **indices, int **nearestNeigh, int n)
{
	int i,j=0;

	for (i=0;i<n;i++)
	{
		if ((dist[i] < DISTANCE_THRESHOLD && prob[i] > EPSILON)) // > CHOSEN_SIGNIFICANCE_LEVEL) || dist[i] == 0.0) 
		{
			nearestNeigh[j][0] = indices[i][0];
			nearestNeigh[j][1] = indices[i][1];

			dist[j] = dist[i]; prob[j] = prob[i];
			j++;
		}
	}
}

void getPatch(double **B, double **im, int y, int x)
{
	int i,j;

	for(i=0;i<ps;i++)
	{
		for(j=0;j<ps;j++)
		{
			B[i][j] = im[y+i][x+j];
		}
	}
}

void addPatch(double **B, double **im2, int **numcount,int y, int x)
{
	int i,j;

	for (i=0;i<ps;i++)
	{
		for(j=0;j<ps;j++)
		{
			im2[y+i][x+j] += B[i][j];
			numcount[y+i][x+j] += 1;			
		}
	}
}

void addPatchWeight (double **B, double **im2, double **numcount,int y, int x, double w)
{
	int i,j;

	for (i=0;i<ps;i++)
	{
		for(j=0;j<ps;j++)
		{
			im2[y+i][x+j] += B[i][j]*w;
			numcount[y+i][x+j] += w;			
		}
	}
}

void addPatchUneven (double **B, double **im2, double **numcount,int y, int x)
{
	int i,j;
	double w;

	for (i=0;i<ps;i++)
	{
		for(j=0;j<ps;j++)
		{
			w = sqrt((i-ps/2)*(i-ps/2) + (j-ps/2)*(j-ps/2));
			if (w != 0) w = 1.0/w; else w = 1;

			im2[y+i][x+j] += B[i][j]*w;
			numcount[y+i][x+j] += w;			
		}
	}
}


void print_to_file (double **B, int y, int x, FILE *fp)
{
	int k1,k2,count=0;

	for (k1=0;k1<ps;k1++)
	{
		for(k2=0;k2<ps;k2++)
		{
			fprintf (fp,"\n%d %d %.3lf %d",y+k1,x+k2,B[k1][k2],count);
			count++;
		}
	}
}

void selectbestsigma_denoise_routine (double **im, double **refim, double **opim, int H, int W, 
double smallest_sigma, double biggest_sigma, char *filename, double **imt, int ps)
{
	double sig_step = 2.0,sig;
	int num_sigmas = (int)((biggest_sigma-smallest_sigma)/sig_step);
	double pr;
	double data1[ps*ps],data2[ps*ps],ksstat,prob[(2*SR+1)*(2*SR+1)];
	// useful variables
	int i,j,k,l,k1,k2,num_neighbors=0,num_neighbors2,minindex_y,minindex_x,s;
	int left,right,top,bottom,count,p1,p2,count2;
	// variables for data structures
	double a,w,*distances,mindist;
	int **indices, **nearestNeigh,nonzero;
	double **numcount;
	double **F,**G,**B,**BB1,**BB2,**BT,**V1,**V2,**V1T,**temp1,**S;
	// time variables
	time_t start, end;
	double dtime;
	struct timeval tv;
	struct timezone tz;
	// GSL random number routines
       const gsl_rng_type * T;
       gsl_rng * r;
	// GSL eigenvector routines
	gsl_eigen_symmv_workspace *ws;
	gsl_matrix *gF, *gG, *gV1, *gV2;
	gsl_vector *gEv,*work;
	FILE *fp;
	double mse[num_sigmas+1],psnr[num_sigmas+1],Kval[num_sigmas+1],Pval[num_sigmas+1],cc[num_sigmas+1];
	char opfname[600],fname[500];
	double minval; 
	int minindex;

	sprintf (fname,"scales_noiseness%.2lf_%s.txt",SIGMA,filename);
	fp = fopen (fname,"a"); if (!fp) { fp = fopen (fname,"w"); }
	printf ("\nOpening the file %s",fname);

	printf ("\nNum sigmas = %d, ps = %d",num_sigmas,ps); fflush(stdout);

	for (s = 0; s <= num_sigmas; s++)
	{
		sig = s*sig_step+smallest_sigma;
		printf ("\n****************************************");	
		printf ("\nSCALE = %lf, s = %d",sig,s); 
		fflush(stdout);

		THRESHOLD = sig*sqrt(4*log(ps));

		if (NOISEMODEL == GAUSSIAN_NOISE) DISTANCE_THRESHOLD = 3*sig*sig*ps*ps;
		else if (NOISEMODEL == NEGATIVE_EXPONENTIAL_NOISE) DISTANCE_THRESHOLD = 5*sig*sig*ps*ps;
		printf ("\nDistance threshold = %lf",DISTANCE_THRESHOLD);
		printf ("\nTHRESHOLD = %lf",THRESHOLD);

		//if (VERBOSE) {printf ("\nAllocating data structures\n"); fflush (stdout);}
		F = allocate_2d_double(ps,ps,'0');
		G = allocate_2d_double(ps,ps,'0');
		B = allocate_2d_double(ps,ps,'0');
		BB1 = allocate_2d_double(ps,ps,'0');
		BB2 = allocate_2d_double(ps,ps,'0');
		BT = allocate_2d_double(ps,ps,'0');
		V1 = allocate_2d_double(ps,ps,'0');
		V2 = allocate_2d_double(ps,ps,'0');
		V1T = allocate_2d_double(ps,ps,'0');
		temp1 = allocate_2d_double(ps,ps,'0');
		S = allocate_2d_double(ps,ps,'0');

		nearestNeigh = allocate_2d_int((2*SR+1)*(2*SR+1),2,'0');
		distances = (double*)calloc((2*SR+1)*(2*SR+1),sizeof(double));
		indices = allocate_2d_int((2*SR+1)*(2*SR+1),2,'0');
		numcount = allocate_2d_double(H,W,'0');

		ws = gsl_eigen_symmv_alloc (4*ps);
		gF  = gsl_matrix_calloc (ps,ps);
		gG  = gsl_matrix_calloc (ps,ps);
		gV1  = gsl_matrix_calloc (ps,ps);
		gV2  = gsl_matrix_calloc (ps,ps);
		gEv = gsl_vector_calloc (ps);

		for(i=0;i<H;i++)
		{
			for(j=0;j<W;j++)
			{
				opim[i][j] = 0;
				numcount[i][j] = 0;
			}
		}

	
		for (i=0;i<=H-ps;i=i+1)
		{
			//if (VERBOSE) {printf ("%d ",i); fflush(stdout);}
			for(j=0;j<=W-ps;j=j+1)
			{
				num_neighbors = 0; num_neighbors2 = 0;
				top = i-SR; if (top < 0) top = 0;
				left = j-SR; if (left < 0) left = 0;
				bottom = i+SR; if (bottom > H-ps)  bottom = H-ps;
				right = j+SR; if (right > W-ps) right = W-ps;
						
				count = 0;
				mindist = MY_INFINITY;
				for (k=top;k<=bottom;k++)
				{
					for(l=left;l<=right;l++)
					{
						distances[count] = 0;
						indices[count][0] = k; indices[count][1] = l;
	
						count2 = 0;					
						for (p1=0;p1<ps;p1++)
						{
							for(p2=0;p2<ps;p2++)
							{
								a = im[i+p1][j+p2]-im[k+p1][l+p2];
								distances[count] += a*a;
								if (distances[count] > DISTANCE_THRESHOLD) goto label1;
								data1[count2] = a; 
								count2++;
							}
						}
						if (distances[count] < DISTANCE_THRESHOLD) 
						{
							num_neighbors2++;
							if (distances[count] == 0 || ((i==k) && (j==l))) pr = 1; 
							else if (NOISEMODEL == GAUSSIAN_NOISE) 
								ksGaussian (data1,ps*ps,&ksstat,&pr,sig*sqrt(2.0));
							else if (NOISEMODEL == NEGATIVE_EXPONENTIAL_NOISE) 
								ksLaplacian (data1,ps*ps,&ksstat,&pr,1/sig);	
							if (pr > EPSILON) 
							{
								nearestNeigh[num_neighbors][0] = k;
								nearestNeigh[num_neighbors][1] = l;
								prob[num_neighbors] = pr;
								num_neighbors++;
							}
						}
	
						if (distances[count] < mindist) 
						{ mindist = distances[count]; minindex_y = indices[count][0]; minindex_x = indices[count][1];}

label1:					count++;
					}
				}
				if (num_neighbors == 0)
				{
					num_neighbors = 1;
					nearestNeigh[0][0] = minindex_y;  nearestNeigh[0][1] = minindex_x; prob[0] = 1;
				}
	
				// set the covariance matrices to zero first
				for (k1=0;k1<ps;k1++)
				{
					for(k2=0;k2<ps;k2++)
					{
						F[k1][k2] = 0;
						G[k1][k2] = 0;						
					}
				}
	
				if (num_neighbors > 5)
				{
					// compute the covariance matrices F and G
					w = 0;
					for (k=0;k<num_neighbors;k++)
					{
						getPatch(B,im,nearestNeigh[k][0],nearestNeigh[k][1]);
						matrix_transpose(B,BT,ps,ps);				
						matrix_multiply(B,BT,BB1,ps,ps,ps);
						matrix_multiply(BT,B,BB2,ps,ps,ps);
	
						for (k1=0;k1<ps;k1++)
						{
							for(k2=0;k2<ps;k2++)
							{
								F[k1][k2] =F[k1][k2] + BB1[k1][k2]*prob[k];
								G[k1][k2] =G[k1][k2] + BB2[k1][k2]*prob[k];						
							}
						}
						w += prob[k];
					}
	
					// divide by w
					for (k1=0;k1<ps;k1++)
					{
						for(k2=0;k2<ps;k2++)
						{
							F[k1][k2] /= w;
							G[k1][k2] /= w;						
						}
					}
	
				       for (k1 = 0; k1 < ps; k1++)
					 {
				       	for (k2=0;k2< ps; k2++)
						{		        
					   		gsl_matrix_set (gF, k1, k2, F[k1][k2]);
					   		gsl_matrix_set (gG, k1, k2, G[k1][k2]);
						}
					}
		
					gsl_eigen_symmv (gF,gEv,gV1,ws);
					gsl_eigen_symmv (gG,gEv,gV2,ws);
		
				       for (k1 = 0; k1 < ps; k1++)
					 {
				       	for (k2=0;k2< ps; k2++)
						{		  
							V1[k1][k2] = gsl_matrix_get (gV1, k1, k2);
							V2[k1][k2] = gsl_matrix_get (gV2, k1, k2);
						}
					}
				}
				else
				{
					getDCTFilter (V1,ps);
					matrix_transpose(V1,V2,ps,ps);
				}
	
				k = 0; nearestNeigh[k][0] = minindex_y; nearestNeigh[k][1] = minindex_x;
				{
					// compute V1'*patch*V2 - store in S
					getPatch(B,im,nearestNeigh[k][0],nearestNeigh[k][1]);
					matrix_transpose(V1,V1T,ps,ps); // V1T = V1'
					matrix_multiply(V1T,B,temp1,ps,ps,ps); // temp1 = V1T*B
					matrix_multiply(temp1,V2,S,ps,ps,ps);	// S = temp1*V2
	
					nonzero = ps*ps;
					for(k1=0;k1<ps;k1++)
					{
						for(k2=0;k2<ps;k2++)
						{
							if (fabs(S[k1][k2]) < THRESHOLD) { S[k1][k2] = 0; nonzero = nonzero-1;}
						}
					}	
	
					matrix_multiply(V1,S,temp1,ps,ps,ps); // temp1 = V1*S
					matrix_transpose(V2,V1T,ps,ps); // V1T = V2'
					matrix_multiply(temp1,V1T,B,ps,ps,ps); //B = temp1*V2T				
						
					addPatchWeight(B,opim,numcount,nearestNeigh[k][0],nearestNeigh[k][1],1.0);
				}
			}
		}

		for(i=0;i<H;i++)
		{	
			for(j=0;j<W;j++)
			{
				a = opim[i][j] /(double)numcount[i][j];
	
				if (numcount[i][j] == 0)
				{
					if (VERBOSE) {printf ("\n%lf %lf [%d %d]",opim[i][j],numcount[i][j],i,j); fflush(stdout);}
				}
	
				opim[i][j] = a;		
			}		
		}	

		mse[s] = MSE(opim,imt,H,W); psnr[s] = PSNR(mse[s]); 
		printf ("\nAfter denoising, MSE = %lf, PSNR = %lf",mse[s],psnr[s]); 
		fflush (stdout);

		sprintf (opfname,"smoothed_P%dx%d_Sig%.2lf_D%.2lf_iter%d_%s",ps,ps,sig,DISTANCE_THRESHOLD,i,filename);
		writePGM (opfname, opim, H, W);
		compute_residuals (refim,im,opim,H,W);
		sprintf (opfname,"residual_P%dx%d_Sig%.2lf_D%.2lf_iter%d_%s",ps,ps,sig,DISTANCE_THRESHOLD,i,filename);
		writePGM (opfname, refim, H, W);

		// the above routine to compute the residuals is for normalized residuals only which are useful for visualization purpose.
		// we are interested purely in normalized residuals.
		compute_residuals_nonnormalized (refim,im,opim,H,W);

		KSTest_Noiseness_multiscale (refim,H,W,8,16,&Kval[s],&Pval[s]);
		cc[s] = corrcoeff_Noiseness_multiscale(refim,H,W,8,16);
		printf ("\nK = %lf, P = %lf, cc = %lf",Kval[s],Pval[s],cc[s]);
		fflush(stdout);

		fprintf (fp,"%s %d %lf %lf\n",filename,ps,psnr[s],cc[s]);

		free_2d_double(F,ps);
		free_2d_double(G,ps);
		free_2d_double(B,ps);
		free_2d_double(BB1,ps);
		free_2d_double(BB2,ps);
		free_2d_double(BT,ps);
		free_2d_double(V1,ps);
		free_2d_double(V2,ps);
		free_2d_double(V1T,ps);
		free_2d_double(temp1,ps);
		free_2d_double(S,ps);
		free (distances);
		free_2d_int(indices,2*SR+1);
		free_2d_double(numcount,H);
		free_2d_int(nearestNeigh,2*SR+1);

		gsl_eigen_symmv_free (ws);
		gsl_matrix_free (gF);
		gsl_matrix_free (gG);
		gsl_matrix_free (gV1);
		gsl_matrix_free (gV2);
		gsl_vector_free (gEv);
	}	

	minval = INFINITY;
	count = 0;
	for(s=0;s<=num_sigmas;s++)
	{
		if (mse[s] <= minval) { minval = mse[s]; minindex = s;}
	}
	printf ("\nMSE: %lf %lf ",minval,smallest_sigma + minindex*sig_step); fflush(stdout);
	fprintf (fp,"%s %lf %lf ",filename,minval,smallest_sigma + minindex*sig_step);

	minval = INFINITY;
	for(s=0;s<=num_sigmas;s++)
	{
		if (cc[s] <= minval) { minval = cc[s]; minindex = s;}
	}
	printf ("\nCC: %lf %lf ",minval,smallest_sigma + minindex*sig_step); fflush(stdout);
	fprintf (fp,"%lf %lf ",minval,smallest_sigma + minindex*sig_step);

	minval = INFINITY;
	for(s=0;s<=num_sigmas;s++)
	{
		if (Kval[s] <= minval) { minval = Kval[s]; minindex = s;}
	}
	printf ("\nK: %lf %lf ",minval,smallest_sigma + minindex*sig_step); fflush(stdout);
	fprintf (fp,"%lf %lf ",minval,smallest_sigma + minindex*sig_step);

	minval = INFINITY;
	for(s=0;s<=num_sigmas;s++)
	{
		if (Pval[s] <= minval) { minval = Pval[s]; minindex = s;}
	}
	printf ("\nlog P: %lf %lf ",minval,smallest_sigma + minindex*sig_step); fflush(stdout);
	fprintf (fp,"%lf %lf ",minval,smallest_sigma + minindex*sig_step);

	fclose (fp);
}


int main (int argc, char **argv)
{
	// useful variables
	char opfname[500];
	int i,j,choiceDCT,numIters;
	int H,W;
	FILE *fp;
	// variables for data structures
	double **im,**im2,**imt,**refim;
	// quality metrics
	double mse, psnr,noisy_mse,noisy_psnr,noisy_ssim,ssim;
	// time vars
	time_t start,end;
	double dtime;
     
	im = readPGM (argv[1],&H,&W);
	imt = readPGM(argv[2],&H,&W);
	SIGMA = atof(argv[3]);
	ps = atoi(argv[4]);
	choiceDCT = atoi(argv[5]); // 0 for DCT ; anything else will use noisy values
	NOISEMODEL = atoi(argv[6]);
	printf ("\nNoise std. dev. = %lf, K = %d, patchsize = %d, file %s",SIGMA,K,ps,argv[2]);

	// set critical global variables. SR is the radius of the search window to look for similar patches
	// K is the number of nearest neighbors.
	SR = 20;
    
	im2 = allocate_2d_double(H,W,'0');
	refim = allocate_2d_double(H,W,'0');

	if (choiceDCT == 0)
	{
		if (VERBOSE) {printf ("\nDoing the DCT filtering"); fflush (stdout);}
		dctFilter (im,refim,ps,H,W);
		if (VERBOSE) {printf ("\nWriting DCT file output"); fflush (stdout);}
		sprintf (opfname,"DCT_P%dx%d_Sig%.2lf_K%d_iter%d_%s",ps,ps,SIGMA,K,i,argv[2]);
		writePGM (opfname, refim, H, W);
		mse = MSE(imt,refim,H,W); psnr = PSNR(mse); 
		printf ("\nAfter DCT, MSE = %lf, PSNR = %lf",mse,psnr); fflush (stdout);

		numIters = 1;
	}
	else
	{
		copyimage (refim,im,H,W);
		numIters = 1;
	}

	for (i = 0; i < numIters;i++)
	{
		printf ("\n"); fflush(stdout);
		time (&start);
		selectbestsigma_denoise_routine (im,refim,im2,H,W,2,40,argv[2],imt,ps);
		time (&end);
		dtime = difftime (end,start);
	}

}
