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

void getDCTFilter (double **A, int ps)
{
	int i,j;
	double sps = sqrt((double)(ps));
	double SQRT2 = sqrt(2.0);

	for (i=0;i<ps;i++)
	{
		for(j=0;j<ps;j++)
		{
			A[j][i] = cos(M_PI*(i)*(2*(j+1)-1)/(2*ps));
			if (i == 0) 
			{
				A[j][i] = A[j][i]/sps; 
			}
			else 
			{ 
				A[j][i] = A[j][i]*SQRT2/sps;
			}
		}
	}

//	printmatrix (A,ps,ps);

// SEE MATLAB HELP for the exact DCT formula that we implemented.
/*U = zeros(ps,ps);
for i=1:ps
   U(:,i) = cos(pi*(2*[1:ps]-1)*(i-1)/(2*ps));
   if i == 1, U(:,i) = U(:,i)/sqrt(ps); else U(:,i) = U(:,i)*sqrt(2/ps); end
end
*/
}


void dctFilter (double **im[3], double **im2[3], int ps, int H, int W)
{
	int i,j,k1,k2,k;
	double a;
	double **U, **UT,**B, **T1, **S;
	int **numcount[3];

	U = allocate_2d_double(ps,ps,'0');
	UT = allocate_2d_double(ps,ps,'0');
	B = allocate_2d_double(ps,ps,'0');
	T1 = allocate_2d_double(ps,ps,'0');
	S = allocate_2d_double(ps,ps,'0');

	// memory for DCT filtered output is assumed to be allocated.
	for (k=0;k<3;k++)
	{
		numcount[k] = allocate_2d_int(H,W,'0');
	}

	getDCTFilter (U,ps);
	matrix_transpose (U,UT,ps,ps);


	printf ("\n"); fflush(stdout);
	for(i=0;i<H-ps+1;i=i+1)
	{
		if (VERBOSE) {printf ("%d ",i); fflush (stdout);}
		for(j=0;j<W-ps+1;j=j+1)
		{
			for(k=0;k<3;k++)
			{
				getPatch(B,im,i,j,k);
			
				matrix_multiply(UT,B,T1,ps,ps,ps); // T1 = UT*B
				matrix_multiply(T1,U,S,ps,ps,ps); // S = T1*U

				for(k1=0;k1<ps;k1++)
				{
					for(k2=0;k2<ps;k2++)
					{
						if (fabs(S[k1][k2]) < THRESHOLDS[k]) S[k1][k2] = 0;
					}
				}

				matrix_multiply(U,S,T1,ps,ps,ps); // T1 = U*S
				matrix_multiply(T1,UT,B,ps,ps,ps); //B = T1*U'				
						
				addPatch(B,im2,numcount[k],i,j,k);
			}
		}
	}

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			for(k=0;k<3;k++)
			{
				a = im2[k][i][j] /((double)numcount[k][i][j]);

				if (numcount[k][i][j] == 0)
				{
					printf ("\n%lf %d [%d %d]",im2[k][i][j],numcount[k][i][j],i,j); fflush(stdout);
				}
	
				im2[k][i][j] = a;	
			}	
		}		
	}

	for (k=0;k<3;k++) free_2d_int(numcount[k],H);
	free_2d_double(U,ps);
	free_2d_double(UT,ps);
	free_2d_double(B,ps);
	free_2d_double(S,ps);
	free_2d_double(T1,ps);
}
