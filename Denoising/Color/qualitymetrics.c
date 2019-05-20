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

double MSE (double **im1[3], double **im2[3], int H, int W)
{
	int i,j;
	double mse = 0,a,b,c;

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			a = im1[0][i][j] - im2[0][i][j];
			b = im1[1][i][j] - im2[1][i][j];
			c = im1[2][i][j] - im2[2][i][j];
			mse += a*a + b*b + c*c;
		}
	}

	mse /= (H*W*3);

	return mse;
}

double PSNR (double mse)
{
	return 10*log10(255*255/mse);
}

void compute_residuals (double **res[3], double **imnoisy[3], double **imdenoised[3], int H, int W)
{
	int i,j,k;
	double minval,maxval;

	for(k=0;k<3;k++)
	{
		maxval = -MY_INFINITY;
		minval = MY_INFINITY;
		for(i=0;i<H;i++)
		{
			for(j=0;j<W;j++)
			{
				res[k][i][j] = imnoisy[k][i][j]-imdenoised[k][i][j];
				if (res[k][i][j] < minval) minval = res[k][i][j];
				if (res[k][i][j] > maxval) maxval = res[k][i][j];
			}
		}
	
		if (minval == maxval) return;

		for(i=0;i<H;i++)
		{
			for(j=0;j<W;j++)
			{
				res[k][i][j] = 255*(res[k][i][j]-minval)/(maxval-minval);
			}
		}	
	}
}



