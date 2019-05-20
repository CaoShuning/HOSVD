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

void denoise_routine (double **im, double **refim, double **opim, int H, int W, int ps)
{
	// useful variables
	int i,j,k,l,k1,k2;
	int left,right,top,bottom,count,p1,p2;
	// variables for data structures
	double a,*distances;
	int **indices, **minK,**numcount;
	double **F1,**F2,**F3,**B,**U1,**U2,**U3,**U1T,**U23,**U23T,**	temp1,**S,***patchstack,**S23,**F1T,**F2T,**F3T;
	double **FF1,**FF2,**FF3;
	double **meanpatch;
	// GSL eigenvector routines
	gsl_matrix *gFF1, *gFF2, *gFF3,*gV1, *gV2, *gV3;
	gsl_eigen_symmv_workspace  *ws, *wsK;
	gsl_vector *gEv, *gEvK;


	if (VERBOSE) {printf ("\nAllocating data structures\n"); fflush (stdout);}
	FF1 = allocate_2d_double(ps,ps,'0');
	FF2 = allocate_2d_double(ps,ps,'0');
	B = allocate_2d_double(ps,ps,'0');
	U1 = allocate_2d_double(ps,ps,'0');
	U2 = allocate_2d_double(ps,ps,'0');
	U1T = allocate_2d_double(ps,ps,'0');
	meanpatch = allocate_2d_double(ps,ps,'0');

	distances = (double*)calloc((2*SR+1)*(2*SR+1),sizeof(double));
	indices = allocate_2d_int((2*SR+1)*(2*SR+1),2,'0');
	minK = allocate_2d_int((2*SR+1)*(2*SR+1),2,'0');
	numcount = allocate_2d_int(H,W,'0');

	gFF1  = gsl_matrix_calloc (ps,ps);
	gFF2  = gsl_matrix_calloc (ps,ps);
	gV1  = gsl_matrix_calloc (ps,ps);
	gV2  = gsl_matrix_calloc (ps,ps);
	ws = gsl_eigen_symmv_alloc (4*ps);
	gEv = gsl_vector_calloc (ps);

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			opim[i][j] = 0;
			numcount[i][j] = 0;
		}
	}

	DISTANCE_THRESHOLD = 3*ps*ps*SIGMA*SIGMA;
	for (i=0;i<=H-ps;i++)
	{
		if (VERBOSE) {printf ("%d ",i); fflush(stdout);}

		for(j=0;j<=W-ps;j++)
		{
			top = i-SR; if (top < 0) top = 0;
			left = j-SR; if (left < 0) left = 0;
			bottom = i+SR; if (bottom > H-ps)  bottom = H-ps;
			right = j+SR; if (right > W-ps) right = W-ps;
					
			count = 0;
			for (k=top;k<=bottom;k++)
			{
				for(l=left;l<=right;l++)
				{
					distances[count] = 0;
					minK[count][0] = indices[count][0] = k; 
					minK[count][1] = indices[count][1] = l;

					for (p1=0;p1<ps;p1++)
					{
						for(p2=0;p2<ps;p2++)
						{
							a = refim[i+p1][j+p2]-refim[k+p1][l+p2];
							distances[count] += a*a;
						}
					}
					if (distances[count] < DISTANCE_THRESHOLD) count++;
				}
			}
			// minimum K distances
			if (count > 30) { K = 30; getMinK (distances,indices,minK,count);}
			else K = count;

			F1 = allocate_2d_double(ps,ps*K,'0');
			F2 = allocate_2d_double(ps,ps*K,'0');
			F3 = allocate_2d_double(K,ps*ps,'0');
			FF3 = allocate_2d_double(K,K,'0');
			F1T = allocate_2d_double(ps*K,ps,'0');
			F2T = allocate_2d_double(ps*K,ps,'0');
			F3T = allocate_2d_double(ps*ps,K,'0');
			U3 = allocate_2d_double(K,K,'0');
			U23 = allocate_2d_double(K*ps,K*ps,'0');
			U23T = allocate_2d_double(K*ps,K*ps,'0');
			temp1 = allocate_2d_double(ps,K*ps,'0');
			S = allocate_2d_double(ps,K*ps,'0');
			patchstack = allocate_3d_double(ps,ps,K,'0');
			S23 = allocate_2d_double(ps,K*ps,'0');
			gFF3  = gsl_matrix_calloc (K,K);
			gV3  = gsl_matrix_calloc (K,K);
			wsK = gsl_eigen_symmv_alloc (4*K);
			gEvK = gsl_vector_calloc (K);

			// initialize mean patch to 0
			for(k1=0;k1<ps;k1++)
			{
				for(k2=0;k2<ps;k2++)
				{
					meanpatch[k1][k2] = 0;
				}
			}

			/*// compute the mean patch
			for (k=0;k<K;k++)
			{
				getPatch(B,refim,minK[k][0],minK[k][1]);
				for(k1=0;k1<ps;k1++)
				{
					for(k2=0;k2<ps;k2++)
					{
						meanpatch[k1][k2] += B[k1][k2];
					}
				}
			}			

			for(k1=0;k1<ps;k1++)
			{
				for(k2=0;k2<ps;k2++)
				{
					meanpatch[k1][k2] /= K;
				}
			}*/


			for (k=0;k<K;k++)
			{
				getPatch(B,refim,minK[k][0],minK[k][1]);
				for (k1=0;k1<ps;k1++)
				{
					for(k2=0;k2<ps;k2++)
					{
						// patchstack contains the patches stacked up together in an array of size ps x ps x K
						patchstack[k1][k2][k] = B[k1][k2]-meanpatch[k1][k2];		
					}
				}
			}

			getFolding1 (patchstack,F1,ps,ps,K); //F1 has size ps by (ps x K)
			getFolding2 (patchstack,F2,ps,ps,K);//F2 has size ps by (ps x K)
			getFolding3 (patchstack,F3,ps,ps,K);//F3 has size K by (ps x ps)
	
			getDCTFilter(U1,ps);
			getDCTFilter(U2,ps);
			getDCTFilter(U3,K);

			kronecker_product(U23,U2,U3,ps,K);

			matrix_transpose(U1,U1T,ps,ps); // V1T = V1'
			matrix_multiply(U1T,F1,temp1,ps,ps,K*ps); // temp1 = V1T*B
			matrix_multiply(temp1,U23,S23,ps,K*ps,K*ps);	// S = temp1*V23

			THRESHOLD = sqrt(2*log(ps*ps*K))*SIGMA; 
			//THRESHOLD = 3*SIGMA; 
			for(k1=0;k1<ps;k1++)
			{
				for(k2=0;k2<K*ps;k2++)
				{
					if (fabs(S23[k1][k2]) < THRESHOLD) { S23[k1][k2] = 0;}
				}
			}	

			matrix_multiply(U1,S23,temp1,ps,ps,K*ps); // temp1 = V1*S
			matrix_transpose(U23,U23T,K*ps,K*ps); // V1T = V2'
			matrix_multiply(temp1,U23T,F1,ps,K*ps,K*ps); //B = temp1*V2T				
			getReverseFolding1 (patchstack,F1,ps,ps,K);

			for(k=0;k<K;k++)
			{
				for(k1=0;k1<ps;k1++)
				{
					for(k2=0;k2<ps;k2++)
					{
						B[k1][k2] = patchstack[k1][k2][k]+meanpatch[k1][k2];
					}
				}
				addPatch(B,opim,numcount,minK[k][0],minK[k][1]);
			}	

			free_2d_double(F1,ps);
			free_2d_double(F2,ps);
			free_2d_double(F3,K);
			free_2d_double(FF3,K);
			free_2d_double(F1T,ps*K);
			free_2d_double(F2T,ps*K);
			free_2d_double(F3T,ps*ps);
			free_2d_double(U3,K);
			free_2d_double(U23,ps*K);
			free_2d_double(U23T,ps*K);
			free_2d_double(temp1,ps);
			free_2d_double(S,ps);
			free_3d_double(patchstack,ps,ps);
			free_2d_double(S23,ps);

			gsl_matrix_free (gFF3);
			gsl_matrix_free (gV3);
			gsl_eigen_symmv_free (wsK);
			gsl_vector_free (gEvK);
		}
	}

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			if (numcount[i][j] == 0)
			{
				if (VERBOSE) {printf ("\n%lf %d [%d %d]",opim[i][j],numcount[i][j],i,j); fflush(stdout);}
				a = im[i][j];
			}
			else
			{
				a = opim[i][j] /(double)numcount[i][j];
			}

			opim[i][j] = a;		
		}		
	}
	
	free_2d_double(FF1,ps);
	free_2d_double(FF2,ps);
	free_2d_double(B,ps);
	free_2d_double(U1,ps);
	free_2d_double(U2,ps);
	free_2d_double(U1T,ps);
	free (distances);
	free_2d_int(indices,(2*SR+1)*(2*SR+1));
	free_2d_int(minK,(2*SR+1)*(2*SR+1));
	free_2d_int(numcount,H);
	free_2d_double(meanpatch,ps);

	gsl_matrix_free (gFF1);
	gsl_matrix_free (gFF2);
	gsl_matrix_free (gV1);
	gsl_matrix_free (gV2);
	gsl_eigen_symmv_free (ws);
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
     
	if (argc < 5)
	{
		printf ("\nUSAGE: ./denoise_hosvd_dynamicK.o <noisy file> <clean file> <SIGMA> <patchsize> <choiceDCT (1)/(0)>");
		exit (0);
	}

	im = readPGM (argv[1],&H,&W);
	imt = readPGM(argv[2],&H,&W);
	SIGMA = atof(argv[3]);
	ps = atoi(argv[4]);
	choiceDCT = atoi(argv[5]); // 0 for DCT ; anything else will use noisy values
	printf ("\nNoise std. dev. = %lf, K = %d, patchsize = %d, file %s, thresh = %lf",SIGMA,K,ps,argv[2],THRESHOLD);
	fflush(stdout);

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
		printf ("\nTime taken = %lf seconds", dtime); fflush (stdout);
 
		sprintf (opfname,"DCTdynamicK_P%dx%d_Sig%.2lf_D%.2lf_iter%d_%s",ps,ps,SIGMA,DISTANCE_THRESHOLD,i,argv[2]);
		writePGM (opfname, im2, H, W);

		noisy_mse = MSE(im,imt,H,W); noisy_psnr = PSNR(noisy_mse); noisy_ssim = MSSIM (im,imt,H,W);
		printf ("\nBefore denoising, MSE = %lf, PSNR = %lf, SSIM = %lf",noisy_mse,noisy_psnr,noisy_ssim); fflush (stdout);
		mse = MSE(im2,imt,H,W); psnr = PSNR(mse); ssim = MSSIM (im2,imt,H,W);
		printf ("\nAfter denoising, MSE = %lf, PSNR = %lf, SSIM = %lf",mse,psnr,ssim); fflush (stdout);

		compute_residuals (refim,im,im2,H,W);
sprintf (opfname,"residualDCTdynamicK_P%dx%d_Sig%.2lf_D%.2lf_iter%d_%s",ps,ps,SIGMA,DISTANCE_THRESHOLD,i,argv[2]);
		writePGM (opfname, refim, H, W);

		copyimage (refim,im2,H,W);

		fp = fopen ("RESULTS_DCT_DYNAMICK.TXT","a");
	// <FILENAME> <K> <ps> <SR> <sigma> <iteration number> <noisy mse> <noisy psnr> <filtered mse> <filtered psnr>
		fprintf (fp,"\n%s %d %d %d %lf %d %lf %lf %lf %lf",argv[2],K,ps,SR,SIGMA,i,noisy_mse,noisy_psnr,mse,psnr);
		fclose (fp);
	}

}
