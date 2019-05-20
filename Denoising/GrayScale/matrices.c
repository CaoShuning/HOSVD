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

void matrix2vector(double **B, double *vec, int n1, int n2)
{
	int i,j,count=0;

	for(i=0;i<n1;i++)
	{
		for(j=0;j<n2;j++)
		{
			vec[count++] = B[i][j];
		}
	}
}

void vector2matrix (double **B, double *vec, int n1, int n2)
{
	int i,j,count=0;

	for(i=0;i<n1;i++)
	{
		for(j=0;j<n2;j++)
		{
			B[i][j] = vec[count++];
		}
	}
}


void GSL2matrix (gsl_matrix *gF, double **F, int n1, int n2)
{
	int k1,k2;

	for (k1 = 0; k1 < n1; k1++)
	{
	     	for (k2=0;k2 < n2; k2++)
		{		  
			F[k1][k2] = gsl_matrix_get (gF, k1, k2);
		}
	}

}

void matrix2GSL (gsl_matrix *gF, double **F, int n1, int n2)
{
	int k1,k2;

	for (k1 = 0; k1 < n1; k1++)
	 {
       	for (k2=0;k2< n2; k2++)
		{		        
	   		gsl_matrix_set (gF, k1, k2, F[k1][k2]);
		}
	}
}

void setmatrix_zero (double **A, int n1, int n2)
{
	int i,j;

	for(i=0;i<n1;i++)
	{
		for (j=0;j<n2;j++)
		{
			A[i][j] = 0;
		}
	}
}

void kronecker_product (double **B, double **A1, double **A2, int n1, int n2)
{
	int i,j,k,l,ind1,ind2;
	double a;

	for (k=0;k<n1*n2;k=k+n2)
	{
		for(l=0;l<n1*n2;l=l+n2)
		{
			ind1 = k/n2;
			ind2 = l/n2;
			a = A1[ind1][ind2];
			for (i=0;i<n2;i++)
			{
				for(j=0;j<n2;j++)
				{
					B[k+i][l+j] = a*A2[i][j];
				}		
			}
		}
	}
}

void vector_outerproduct (double *a, double *b, double **A, int n)
{
	int i,j;

	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			A[i][j] = a[i]*b[j];
		}
	}
}

void pad2DArray (double **im1, double **im2, int H, int W, int padsize)
{
	int i,j,k;	

	im1 = allocate_2d_double (H+padsize,W+padsize,'0');
	k = padsize/2;

	for(i=k;i<H+k;i++)
	{
		for(j=k;j<W+k;j++)
		{
			im1[i][j] = im2[i-k][j-k];
		}
	}
}

// A has size N1 x N2, src has size N2 x 1, dest has size N2 x 1
void matrix_vector_multiply (double **A, double *src, double *dest, int N1, int N2)
{
	int i,j;

	for(i=0;i<N1;i++)
	{
		dest[i] = 0;
		for(j=0;j<N2;j++)
		{
			dest[i] += src[j]*A[	i][j];
		}
	}
}

void matrix_transpose (double **A, double **B, int r, int c)
{
	int i,j;

	for (i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			B[i][j] = A[j][i];
		}
	}
}

// A -> (r1,c1) B -> (c1,c2)
void matrix_multiply (double **A, double **B, double **C, int r1, int c1, int c2)
{
	int i,j,k;

	for(i=0;i<r1;i++)
	{
		for(j=0;j<c2;j++)
		{
			C[i][j] = 0;
			for (k=0;k<c1;k++)
			{
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

void printmatrix (double **A, int r, int c)
{
	int i,j;

	printf ("\n\n"); fflush(stdout);
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			printf ("%lf ",A[i][j]); fflush(stdout);
		}
		printf ("\n");
		fflush (stdout);
	}
	
}
