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

void getAngles(double **V1, double **V2, int n, double *angles)
{
	int k,l;
	double sum;

	for(k=0;k<n;k++)
	{
		sum = 0;
		for(l=0;l<n;l++)
		{
			sum += V1[k][l]*V2[k][l];
		}
		sum = fabs(sum); if (sum > 1.0) sum = 1.0;
		angles[k] = acos(fabs(sum))*180/M_PI;
	}
}

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

void addHypotheses(double ***imhyps, double **B, double **numcount, int y, int x)
{
	int k1,k2;

	for(k1=0;k1<ps;k1++)
	{
		for(k2=0;k2<ps;k2++)
		{
			imhyps[y+k1][x+k2][(int)numcount[y+k1][x+k2]] = B[k1][k2]; 
			numcount[y+k1][x+k2] = numcount[y+k1][x+k2] + 1;
		}
	}
}


void denoise_routine (double **im, double **refim, double **opim, int H, int W, int ps)
{
	double angles_u[ps],angles_v[8];
	int oddeven; 
	double pr;
	double data1[ps*ps],data2[ps*ps],ksstat,prob[(2*SR+1)*(2*SR+1)];
	// useful variables
	int i,j,k,l,k1,k2,num_neighbors=0,num_neighbors2,minindex_y,minindex_x;
	int left,right,top,bottom,count,p1,p2,count2;
	// variables for data structures
	double a,w,*distances,mindist,cm;
	int **indices, **nearestNeigh;
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
	double ***imhyps,meanval,stdval;
	double **DCT_V1, **DCT_V2;

	if (VERBOSE) {printf ("\nAllocating data structures\n"); fflush (stdout);}
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
	imhyps = allocate_3d_double(H,W,ps*ps,'0');

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

	DCT_V1 = allocate_2d_double(ps,ps,'0');
	DCT_V2 = allocate_2d_double(ps,ps,'0');
	getDCTFilter (DCT_V1,ps);
	matrix_transpose(DCT_V1,DCT_V2,ps,ps);	

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			opim[i][j] = 0;
			numcount[i][j] = 0;
		}
	}

	if (NOISEMODEL == GAUSSIAN_NOISE) DISTANCE_THRESHOLD = 3*SIGMA*SIGMA*ps*ps;
	else if (NOISEMODEL == NEGATIVE_EXPONENTIAL_NOISE) DISTANCE_THRESHOLD = 5*SIGMA*SIGMA*ps*ps;
	printf ("\nDistance threshold = %lf",DISTANCE_THRESHOLD);

	//fp = fopen ("DCT_distance.txt","w");
	for (i=0;i<=H-ps;i=i+1)
	{
		if (VERBOSE) {printf ("%d ",i); fflush(stdout);}
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
							//data2[count2] = gsl_ran_gaussian (r,sqrt(2)*SIGMA);;
							count2++;
						}
					}
					if (distances[count] < DISTANCE_THRESHOLD) 
					{
						num_neighbors2++;
						if (distances[count] == 0 || ((i==k) && (j==l))) pr = 1; 
						else if (NOISEMODEL == GAUSSIAN_NOISE) 
							ksGaussian (data1,ps*ps,&ksstat,&pr,SIGMA*sqrt(2.0));
						else if (NOISEMODEL == NEGATIVE_EXPONENTIAL_NOISE) 
							ksLaplacian (data1,ps*ps,&ksstat,&pr,1/SIGMA);	
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

label1:				count++;
				}
			}
			if (num_neighbors == 0)
			{
				num_neighbors = 1;
				nearestNeigh[0][0] = minindex_y;  nearestNeigh[0][1] = minindex_x; prob[0] = 1;
			}

			//getNearestNeighborsProb (prob,distances,indices,nearestNeigh,count);
			//if (VERBOSE) {printf ("\n[%d %d]: %d %d",i,j,num_neighbors,num_neighbors2); fflush(stdout);}
			//printf ("\nnum_neighbors  = %d",num_neighbors); fflush(stdout);

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
					//printf ("\n%lf ",prob[k]); fflush(stdout);
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

//			getAngles(V1,DCT_V1,ps,angles_u);
//			getAngles(V2,DCT_V2,ps,angles_v);
//			for(k=0;k<ps;k++)
//			{
//				fprintf (fp,"%.4lf %.4lf ",angles_u[k], angles_v[k]);
//			}
//			fprintf (fp,"\n");

			k = 0; 
			// this change below made on Oct 8th at 22:51 - error discovered from Mohsen's image
			//nearestNeigh[k][0] = minindex_y; nearestNeigh[k][1] = minindex_x;
			nearestNeigh[k][0] = i; nearestNeigh[k][1] = j;

			//if (minindex_y != i || minindex_x != j) 
			//{ printf ("\n** %d --> %d,  %d --> %d",i,minindex_y,j,minindex_x); exit (0);}
			//for(k=0;k<num_neighbors;k++)
			{
				// compute V1'*patch*V2 - store in S
				getPatch(B,im,nearestNeigh[k][0],nearestNeigh[k][1]);
				matrix_transpose(V1,V1T,ps,ps); // V1T = V1'
				matrix_multiply(V1T,B,temp1,ps,ps,ps); // temp1 = V1T*B
				matrix_multiply(temp1,V2,S,ps,ps,ps);	// S = temp1*V2

				for(k1=0;k1<ps;k1++)
				{
					for(k2=0;k2<ps;k2++)
					{
						if (fabs(S[k1][k2]) < THRESHOLD) S[k1][k2] = 0; 
					}
				}	

				matrix_multiply(V1,S,temp1,ps,ps,ps); // temp1 = V1*S
				matrix_transpose(V2,V1T,ps,ps); // V1T = V2'
				matrix_multiply(temp1,V1T,B,ps,ps,ps); //B = temp1*V2T				
					
				addPatchWeight(B,opim,numcount,nearestNeigh[k][0],nearestNeigh[k][1],1.0);
				//print_to_file (B,nearestNeigh[k][0],nearestNeigh[k][1],fp);
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
	THRESHOLD = SIGMA*(sqrt(4*log(ps)));
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
		denoise_routine (im,refim,im2,H,W,ps);
		time (&end);
		dtime = difftime (end,start);
		printf ("\nTime taken = %lf seconds", dtime);
 
		sprintf (opfname,"smoothed_P%dx%d_Sig%.2lf_D%.2lf_iter%d_%s",ps,ps,SIGMA,DISTANCE_THRESHOLD,i,argv[2]);
		writePGM (opfname, im2, H, W);

		noisy_mse = MSE(im,imt,H,W); noisy_psnr = PSNR(noisy_mse); noisy_ssim = MSSIM (im,imt,H,W);
		printf ("\nBefore denoising, MSE = %lf, PSNR = %lf, SSIM = %lf",noisy_mse,noisy_psnr,noisy_ssim); fflush (stdout);
		mse = MSE(im2,imt,H,W); psnr = PSNR(mse); ssim = MSSIM (im2,imt,H,W);
		printf ("\nAfter denoising, MSE = %lf, PSNR = %lf, SSIM = %lf",mse,psnr,ssim); fflush (stdout);

		compute_residuals (refim,im,im2,H,W);
		sprintf (opfname,"residual_P%dx%d_Sig%.2lf_D%.2lf_iter%d_%s",ps,ps,SIGMA,DISTANCE_THRESHOLD,i,argv[2]);
		writePGM (opfname, refim, H, W);

		copyimage (im,im2,H,W);

		fp = fopen ("RESULTS_DYNAMIC_K_NEW.TXT","a"); 
		if (!fp) 
		{
			fp = fopen ("RESULTS_DYNAMIC_K_NEW.TXT","w");
		}
	// <FILENAME> <distance> <ps> <SR> <sigma> <iteration number> <noisy mse> <noisy psnr> <noisy ssim> <filtered mse> <filtered psnr> <filtered ssim>
		fprintf (fp,"\n%s %f %d %d %lf %d %lf %lf %lf %lf %lf %lf",
		argv[2],DISTANCE_THRESHOLD,ps,SR,SIGMA,i,noisy_mse,noisy_psnr,noisy_ssim,mse,psnr,ssim);
		fclose (fp);
	}

}
