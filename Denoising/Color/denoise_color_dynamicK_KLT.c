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

void copyimage (double **dest[3], double **src[3], int H, int W)
{
	int i,j,k;

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			for(k=0;k<3;k++)
			{
				dest[k][i][j] = src[k][i][j];
			}
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


void getPatch(double **B, double **im[3], int y, int x, int channel)
{
	int i,j;

	for(i=0;i<ps;i++)
	{
		for(j=0;j<ps;j++)
		{
			B[i][j] = im[channel][y+i][x+j];
		}
	}
}

void addPatch(double **B, double **im2[3], int **numcount,int y, int x, int channel)
{
	int i,j;

	for (i=0;i<ps;i++)
	{
		for(j=0;j<ps;j++)
		{
			im2[channel][y+i][x+j] += B[i][j];
			numcount[y+i][x+j] += 1;			
		}
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


void denoise_routine (double **im[3], double **refim[3], double **opim[3], int H, int W, int ps)
{
	// useful variables
	int i,j,k,l,k1,k2,c,num_neighbors=0;;
	int left,right,top,bottom,count,p1,p2;
	// variables for data structures
	double b,*distances,a;
	int **indices, **nearestNeigh,***numcount;
	double **F,**G,**B,**BB1,**BB2,**BT,**V1,**V2,**V1T,**temp1,**S;
	// time variables
	time_t start, end;
	double dtime;
	// GSL random number routines
       const gsl_rng_type * T;
       gsl_rng * r;
	// GSL eigenvector routines
	gsl_eigen_symmv_workspace *ws;
	gsl_matrix *gF, *gG, *gV1, *gV2;
	gsl_vector *gEv;

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

	distances = (double*)calloc((2*SR+1)*(2*SR+1),sizeof(double));
	indices = allocate_2d_int((2*SR+1)*(2*SR+1),2,'0');
	nearestNeigh = allocate_2d_int((2*SR+1)*(2*SR+1),2,'0');
	numcount = allocate_3d_int(3,H,W,'0');

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
			opim[0][i][j] = opim[1][i][j]  = opim[2][i][j] = 0;
			numcount[0][i][j] = numcount[1][i][j] = numcount[2][i][j] = 0;
		}
	}

	for (i=0;i<=H-ps;i++)
	{
		if (VERBOSE) {printf ("%d ",i); fflush(stdout);}
		for(j=0;j<=W-ps;j++)
		{
			for (c = 0; c < 3; c++)
			{
				num_neighbors = 0;
				top = i-SR; if (top < 0) top = 0;
				left = j-SR; if (left < 0) left = 0;
				bottom = i+SR; if (bottom > H-ps)  bottom = H-ps;
				right = j+SR; if (right > W-ps) right = W-ps;
					
				if (c == 0) DISTANCE_THRESHOLD = 
				3*(SIGMA*SIGMA*(KLT_V[0][0]*KLT_V[0][0] + KLT_V[1][0]*KLT_V[1][0] + KLT_V[2][0]*KLT_V[2][0]))*ps*ps;
				else if (c == 1) DISTANCE_THRESHOLD = 
				3*(SIGMA*SIGMA*(KLT_V[0][1]*KLT_V[0][1] + KLT_V[1][1]*KLT_V[1][1] + KLT_V[2][1]*KLT_V[2][1]))*ps*ps;
				else if (c == 2) DISTANCE_THRESHOLD = 
				3*(SIGMA*SIGMA*(KLT_V[0][2]*KLT_V[0][2] + KLT_V[1][2]*KLT_V[1][2] + KLT_V[2][2]*KLT_V[2][2]))*ps*ps;
	
				count = 0;
				for (k=top;k<=bottom;k++)
				{
					for(l=left;l<=right;l++)
					{
						distances[count] = 0;
						indices[count][0] = k; indices[count][1] = l;
		
						for (p1=0;p1<ps;p1++)
						{
							for(p2=0;p2<ps;p2++)
							{
								// DISTANCE MEASURED IN 'KLT' SPACE
								b = im[c][i+p1][j+p2]-im[c][k+p1][l+p2];
								distances[count] += b*b;
							}
						}
						if (distances[count] < DISTANCE_THRESHOLD) num_neighbors++;
						count++;
					}
				}
				getNearestNeighbors (distances,indices,nearestNeigh,count);

				// set the covariance matrices to zero first
				for (k1=0;k1<ps;k1++)
				{
					for(k2=0;k2<ps;k2++)
					{
						F[k1][k2] = 0;
						G[k1][k2] = 0;						
					}
				}

				// compute the covariance matrices F and G
				for (k=0;k<num_neighbors;k++)
				{
					getPatch(B,im,nearestNeigh[k][0],nearestNeigh[k][1],c);
					matrix_transpose(B,BT,ps,ps);				
					matrix_multiply(B,BT,BB1,ps,ps,ps);
					matrix_multiply(BT,B,BB2,ps,ps,ps);

					for (k1=0;k1<ps;k1++)
					{
						for(k2=0;k2<ps;k2++)
						{
							F[k1][k2]   = F[k1][k2] + BB1[k1][k2];
							G[k1][k2]  = G[k1][k2] + BB2[k1][k2];						
						}
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

//				gsl_eigen_symmv (gsl_matrix * A, gsl_vector * eval, gsl_matrix * evec, gsl_eigen_symmv_workspace * w);
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
	
				//for(k=0;k<num_neighbors;k++)
				{
					// compute V1'*patch*V2 - store in S
					getPatch(B,im,i,j,c);
					matrix_transpose(V1,V1T,ps,ps); // V1T = V1'
					matrix_multiply(V1T,B,temp1,ps,ps,ps); // temp1 = V1T*B
					matrix_multiply(temp1,V2,S,ps,ps,ps);	// S = temp1*V2
	
					for(k1=0;k1<ps;k1++)
					{
						for(k2=0;k2<ps;k2++)
						{
							if (fabs(S[k1][k2]) < THRESHOLDS[c]) S[k1][k2] = 0;
						}
					}	
	
					matrix_multiply(V1,S,temp1,ps,ps,ps); // temp1 = V1*S
					matrix_transpose(V2,V1T,ps,ps); // V1T = V2'
					matrix_multiply(temp1,V1T,B,ps,ps,ps); //B = temp1*V2T				
						
					addPatch(B,opim,numcount[c],i,j,c);
				} // close k
			} // close c
		} // close j
	} // close i
 
	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			for(k=0;k<3;k++)
			{
				a = opim[k][i][j] /((double)numcount[k][i][j]);
	
				if (numcount[k][i][j] == 0)
				{
					if (VERBOSE) {printf ("\n%lf %d [%d %d]",opim[k][i][j],numcount[k][i][j],i,j); fflush(stdout);}
				}

				opim[k][i][j] = a;		
			}
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
	free_3d_int(numcount,3,H);

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
	double **im[3],**im2[3],**imt[3],** refim[3],**res[3];
	// quality metrics
	double mse, psnr,noisy_mse,noisy_psnr;
	// time vars
	time_t start,end;
	double dtime;
  
	readPPM (argv[1],&H,&W,im); 
	readPPM (argv[2],&H,&W,imt);
	SIGMA = atof(argv[3]);
	ps = atoi(argv[4]);
	choiceDCT = atoi(argv[5]); // 0 for DCT ; anything else will use noisy values

	noisy_mse = MSE(im,imt,H,W); noisy_psnr = PSNR(noisy_mse); 
	printf ("\nBefore denoising, MSE = %lf, PSNR = %lf",noisy_mse,noisy_psnr); fflush (stdout);

	KLT_V = allocate_2d_double(3,3,'0'); 
	KLT_mu = (double*)calloc(3,sizeof(double));
	// compute the KLT transform
	compute_klt (im,H,W,KLT_V,KLT_mu);
	// convert to KLT space
	rgb2klt (im,H,W,KLT_V,KLT_mu);

	printf ("\nThe KLT transform is:\n");
	printmatrix(KLT_V,3,3);
	fflush (stdout);
	printf ("\nThe mean vector is %lf %lf %lf",KLT_mu[0],KLT_mu[1],KLT_mu[2]);
	fflush (stdout);

	THRESHOLDS[0] = THRESHOLDS[1]  = THRESHOLDS[2] = sqrt(2*log(ps*ps))*SIGMA;

	printf ("\nNoise std. dev. = %lf, patchsize = %d, file %s",SIGMA,ps,argv[2]);

	// set critical global variables. SR is the radius of the search window to look for similar patches
	// K is the number of nearest neighbors.
	SR = 20;
    
	for(i=0;i<3;i++)
	{
		im2[i] = allocate_2d_double(H,W,'0');
		refim[i] = allocate_2d_double(H,W,'0');
		res[i] = allocate_2d_double(H,W,'0');
	}

	if (choiceDCT == 0)
	{
		if (VERBOSE) {printf ("\nDoing the DCT filtering"); fflush (stdout);}
		rgb2ycbcr(im,H,W);
		dctFilter (im,refim,ps,H,W);
		if (VERBOSE) {printf ("\nWriting DCT file output"); fflush (stdout);}
		sprintf (opfname,"DCT_P%dx%d_Sig%.2lf_K%d_iter0_%s",ps,ps,SIGMA,K,argv[2]);
		ycbcr2rgb (refim,H,W);
		writePPM (opfname, refim, H, W);
		mse = MSE(imt,refim,H,W); psnr = PSNR(mse); 
		printf ("\nAfter DCT, MSE = %lf, PSNR = %lf",mse,psnr); fflush (stdout);

		numIters = 0;
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
		//printf ("Time taken = %lf", dtime);
 
		klt2rgb (im2,H,W,KLT_V,KLT_mu);
		klt2rgb (im,H,W,KLT_V,KLT_mu);
		//sprintf (opfname,"smoothed_P%dx%d_Sig%.2lf_K%d_iter%d_%s",ps,ps,SIGMA,K,i,argv[2]);
		sprintf (opfname,"smoothed_P%dx%d_Sig%.2lf_D%.2lf_iter%d_%s",ps,ps,SIGMA,DISTANCE_THRESHOLD,i,argv[2]);
		writePPM (opfname, im2, H, W);

		compute_residuals (res,im,im2,H,W);
		sprintf (opfname,
		"residual_dynamicK_P%dx%d_Sig%.2lf_D%.2lf_iter%d_%s",ps,ps,SIGMA,DISTANCE_THRESHOLD,i,argv[2]);
		writePPM (opfname,res,H,W);

		mse = MSE(im2,imt,H,W); psnr = PSNR(mse); 
		printf ("\nAfter denoising, MSE = %lf, PSNR = %lf",mse,psnr); fflush (stdout);

		copyimage (refim,im2,H,W);
		ycbcr2rgb (refim,H,W);

		fp = fopen ("COLOR_RESULTS_KLT.TXT","a");
	// <FILENAME> <K> <ps> <SR> <sigma> <iteration number> <noisy mse> <noisy psnr> <filtered mse> <filtered psnr>
		fprintf (fp,"\n%s %d %d %d %lf %d %lf %lf %lf %lf",argv[2],K,ps,SR,SIGMA,i,noisy_mse,noisy_psnr,mse,psnr);
		fclose (fp);
	}

}
