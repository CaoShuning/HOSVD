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


// unfold A -> (m1,m2,m3) to B(m1,m2*m3) [m2 copies of m1 x m3 matrices]
void getFolding1 (double ***A, double **B, int m1, int m2, int m3)
{
	int i1,i2,i3;

/*
	for(i=0;i<m2;i++)
	{
		for(j=0;j<m1;j++)
		{
			for(k=0;k<m3;k++)
			{
				B[j][k+i*m3] = A[j][i][k];
			}
		}
	}
*/


	for(i1=0;i1<m1;i1++)
	{
		for(i2=0;i2<m2;i2++)
		{
			for(i3=0;i3<m3;i3++)
			{
				B[i1][i2*m3+i3] = A[i1][i2][i3];
			}
		}
	} 
}


void getReverseFolding1 (double ***A, double **B, int m1, int m2, int m3)
{
	int i1,i2,i3;

	for(i1=0;i1<m1;i1++)
	{
		for(i2=0;i2<m2;i2++)
		{
			for(i3=0;i3<m3;i3++)
			{
				A[i1][i2][i3] = B[i1][i2*m3+i3];
			}
		}
	} 
}

// unfold A -> (m1,m2,m3) to B(m2,m1*m3) [m3 copies of m2 x m1 matrices]
void getFolding2 (double ***A, double **B, int m1, int m2, int m3)
{
	int i1,i2,i3;

	/*
	for(i=0;i<m3;i++)
	{
		for(j=0;j<m2;j++)
		{
			for(k=0;k<m1;k++)
			{
				B[j][k+i*m1] = A[j][k][i];
			}
		}
	}
	*/

	for(i1=0;i1<m1;i1++)
	{
		for(i2=0;i2<m2;i2++)
		{
			for(i3=0;i3<m3;i3++)
			{
				B[i2][i3*m1+i1] = A[i1][i2][i3];
			}
		}
	} 

}


// unfold A -> (m1,m2,m3) to B(m3,m2*m1) [m1 copies of m3 x m2 matrices]
void getFolding3 (double ***A, double **B, int m1, int m2, int m3)
{
	int i1,i2,i3;

/*
	for(i=0;i<m1;i++)
	{
		for(j=0;j<m3;j++)
		{
			for(k=0;k<m2;k++)
			{
				B[j][k+i*m2] = A[i][j][k];
			}
		}
	}
*/
	for(i1=0;i1<m1;i1++)
	{
		for(i2=0;i2<m2;i2++)
		{
			for(i3=0;i3<m3;i3++)
			{
				B[i3][i1*m2+i2] = A[i1][i2][i3];
			}
		}
	} 
}


