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
// formulae from book: JPEG2000 standard for image compression: concepts, algorithms and VLSI ...
//  By Tinku Acharya, Ping-Sing Tsai
void rgb2ycbcr (double **im[3], int H, int W)
{
	int i,j,k;
	double y,cb,cr;
	double R,G,B;

	for (i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			R = im[0][i][j];
			G = im[1][i][j];
			B = im[2][i][j];

			y = 0.299*R + 0.587*G + 0.114*B;
			cb = -0.1687*R - 0.3312*G + 0.5*B ;
			cr = 0.5*R - 0.4186*G - 0.0813*B ;

			im[0][i][j] = y;
			im[1][i][j] = cb;
			im[2][i][j] = cr;
		}
	}	
}

void ycbcr2rgb (double **im[3], int H, int W)
{
	int i,j,k;
	double y,cb,cr;
	double R,G,B;

	for (i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			y = im[0][i][j];
			cb = im[1][i][j];
			cr = im[2][i][j];

			R = y +1.402*cr;
			G = y - 0.344*cb - 0.71414*cr ;
			B = y + 1.7718*cb ;

			im[0][i][j] = R;
			im[1][i][j] = G;
			im[2][i][j] = B;
		}
	}	
}


void compute_klt_local (double ***P1, double ***P2, double ***P3, int N, double **V, double *meanval)
{
	int i,j,k,k1,k2;
	double vec[3],vec2[3];
	double **tempC,**C;
	gsl_matrix *gC,*gV;
	gsl_eigen_symmv_workspace *ws;
	gsl_vector *gEv;

	tempC = allocate_2d_double(3,3,'0');
	C = allocate_2d_double(3,3,'0');

	ws = gsl_eigen_symmv_alloc (4*3);
	gC  = gsl_matrix_calloc (3,3);
	gV  = gsl_matrix_calloc (3,3);
	gEv = gsl_vector_calloc (3);

	meanval[0] = meanval[1] = meanval[2] = 0;
	setmatrix_zero(C,3,3);

	for (k=0;k<N;k++)
	{
		for(i=0;i<ps;i++)
		{
			for(j=0;j<ps;j++)
			{
				meanval[0] += P1[k][i][j];
				meanval[1] += P2[k][i][j];
				meanval[2] += P3[k][i][j];
			}
		}
	}
	for (k=0;k<3;k++)
		meanval[k] /= N*ps*ps;

	for(k=0;k<N;k++)
	{
		for(i=0;i<ps;i++)
		{
			for(j=0;j<ps;j++)
			{
				vec[0] = P1[k][i][j];
				vec[1] = P2[k][i][j];
				vec[2] = P3[k][i][j];
	
				vector_outerproduct (vec,vec,tempC,3);
	
				for(k1=0;k1<3;k1++)	
				{
					for(k2=0;k2<3;k2++)
					{
						C[k1][k2] += tempC[k1][k2];
					}
				}
			}
		}
	}

	for(k1=0;k1<3;k1++)	
	{
		for(k2=0;k2<3;k2++)
		{
			C[k1][k2] /= (ps*ps*N);
		}
	}

	// compute the transformation
	matrix2GSL(gC,C,3,3);
	gsl_eigen_symmv (gC,gEv,gV,ws);
	GSL2matrix(gV,V,3,3);

	free_2d_double (tempC,3);
	free_2d_double (C,3);

	gsl_eigen_symmv_free (ws);
	gsl_matrix_free (gC);
	gsl_matrix_free (gV);
	gsl_vector_free (gEv);
}

void compute_klt (double **im[3], int H, int W, double **V, double *meanval)
{
	int i,j,k,k1,k2;
	double vec[3],vec2[3];
	double **tempC,**C;
	gsl_matrix *gC,*gV;
	gsl_eigen_symmv_workspace *ws;
	gsl_vector *gEv;

	tempC = allocate_2d_double(3,3,'0');
	C = allocate_2d_double(3,3,'0');

	ws = gsl_eigen_symmv_alloc (4*3);
	gC  = gsl_matrix_calloc (3,3);
	gV  = gsl_matrix_calloc (3,3);
	gEv = gsl_vector_calloc (3);

	meanval[0] = meanval[1] = meanval[2] = 0;
	setmatrix_zero(C,3,3);
	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			for (k=0;k<3;k++)
				meanval[k] += im[k][i][j];
		}
	}
	for (k=0;k<3;k++)
		meanval[k] /= H*W;

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			vec[0] = im[0][i][j];
			vec[1] = im[1][i][j];
			vec[2] = im[2][i][j];

			vector_outerproduct (vec,vec,tempC,3);

			for(k1=0;k1<3;k1++)	
			{
				for(k2=0;k2<3;k2++)
				{
					C[k1][k2] += tempC[k1][k2];
				}
			}
		}
	}

	for(k1=0;k1<3;k1++)	
	{
		for(k2=0;k2<3;k2++)
		{
			C[k1][k2] /= (H*W);
		}
	}

	// compute the transformation
	matrix2GSL(gC,C,3,3);
	gsl_eigen_symmv (gC,gEv,gV,ws);
	GSL2matrix(gV,V,3,3);

	free_2d_double (tempC,3);
	free_2d_double (C,3);

	gsl_eigen_symmv_free (ws);
	gsl_matrix_free (gC);
	gsl_matrix_free (gV);
	gsl_vector_free (gEv);
}

void rgb2klt (double **im[3], int H, int W, double **V, double *meanval)
{
	int i,j;
	double vec[3],vec2[3];
	double **VT;

	VT = allocate_2d_double(3,3,'0');

	matrix_transpose(V,VT,3,3);
	// apply the transformation 
	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			vec[0] = im[0][i][j]-meanval[0]; 
			vec[1] = im[1][i][j]-meanval[1]; 
			vec[2] = im[2][i][j]-meanval[2];

			matrix_vector_multiply (VT,vec,vec2,3,3);

			im[0][i][j] = vec2[0]+meanval[0];
			im[1][i][j] = vec2[1]+meanval[1];
			im[2][i][j] = vec2[2]+meanval[2];		
		}
	}

	free_2d_double (VT,3);
}


void klt2rgb (double **im[3], int H, int W, double **V, double *meanval)
{
	int i,j;
	double vec[3],vec2[3];	

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			vec[0] = im[0][i][j]-meanval[0]; 
			vec[1] = im[1][i][j]-meanval[1]; 
			vec[2] = im[2][i][j]-meanval[2];

			matrix_vector_multiply (V,vec,vec2,3,3);

			im[0][i][j] = vec2[0]+meanval[0];
			im[1][i][j] = vec2[1]+meanval[1];
			im[2][i][j] = vec2[2]+meanval[2];					
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void rgb2klt_local (double ***P1, double ***P2, double ***P3, int N, double **V, double *meanval)
{
	int i,j,k;
	double vec[3],vec2[3];
	double **VT;

	VT = allocate_2d_double(3,3,'0');

	matrix_transpose(V,VT,3,3);
	// apply the transformation 
	for(k=0;k<N;k++)
	{
		for(i=0;i<ps;i++)
		{
			for(j=0;j<ps;j++)
			{
				vec[0] = P1[k][i][j]-meanval[0]; 
				vec[1] = P2[k][i][j]-meanval[1]; 
				vec[2] = P3[k][i][j]-meanval[2];
	
				matrix_vector_multiply (VT,vec,vec2,3,3);
	
				P1[k][i][j] = vec2[0]+meanval[0];
				P2[k][i][j] = vec2[1]+meanval[1];
				P3[k][i][j] = vec2[2]+meanval[2];		
			}
		}
	}

	free_2d_double (VT,3);
}


void klt2rgb_local (double ***P1, double ***P2, double ***P3, int N, double **V, double *meanval)
{
	int i,j,k;
	double vec[3],vec2[3];	

	for(k=0;k<N;k++)
	{
		for(i=0;i<ps;i++)
		{
			for(j=0;j<ps;j++)
			{
				vec[0] = P1[k][i][j]-meanval[0]; 
				vec[1] = P2[k][i][j]-meanval[1]; 
				vec[2] = P3[k][i][j]-meanval[2];

				matrix_vector_multiply (V,vec,vec2,3,3);
	
				P1[k][i][j] = vec2[0]+meanval[0];
				P2[k][i][j] = vec2[1]+meanval[1];
				P3[k][i][j] = vec2[2]+meanval[2];					
			}
		}
	}
}




