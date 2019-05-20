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

double MSE (double **im1, double **im2, int r, int c)
{
	int i,j;
	double mse = 0,a;

	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			a = im1[i][j] - im2[i][j];
			mse += a*a;
		}
	}

	mse /= r*c;

	return mse;
}

double PSNR (double mse)
{
	return 10*log10(255*255/mse);
}


double MSSIM (double **im1, double **im2, int r, int c)
{
	int i,j,k1,k2,count = 0;
	double MSSIM = 0,s1,s2,s12,mu1,mu2;
	double c1,c2;
	double tm,ts;

	c1 = 0.01*255; c1 = c1*c1;
	c2 = 0.03*255; c2 = c2*c2;
	for(i=0;i<=r-ps;i++)
	{
		for(j=0;j<=c-ps;j++)
		{
			count++;
			mu1 = mu2 = 0;
			s12 = s1 = s2 = 0;

			for(k1=i;k1<i+ps;k1++)
			{
				for(k2=j;k2<j+ps;k2++)
				{
					mu1 += im1[k1][k2]; mu2 += im2[k1][k2];
				}
			}
			mu1 /= (ps*ps); mu2 /= (ps*ps);

			for(k1=i;k1<i+ps;k1++)
			{
				for(k2=j;k2<j+ps;k2++)
				{
					s1 += (im1[k1][k2] - mu1)* (im1[k1][k2] - mu1);
					s2 += (im2[k1][k2] - mu2)* (im2[k1][k2] - mu2);
					s12 += (im1[k1][k2] - mu1)* (im2[k1][k2] - mu2);
				}
			}
			s1 /= (ps*ps); s2 /= (ps*ps); s12 /= (ps*ps);
		
			tm = (2*mu1*mu2 + c1)/(mu1*mu1+mu2*mu2 + c1);
			ts = (2*s12 + c2)/(s1+s2 + c2);
			MSSIM += tm*ts;
		}
	}

	MSSIM /= count;
	return MSSIM;
}


void compute_residuals (double **res, double **imnoisy, double **imdenoised, int H, int W)
{
	int i,j;
	double maxval = -MY_INFINITY, minval = MY_INFINITY;

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			res[i][j] = imnoisy[i][j]-imdenoised[i][j];
			if (res[i][j] < minval) minval = res[i][j];
			if (res[i][j] > maxval) maxval = res[i][j];
		}
	}

	if (minval == maxval) return;

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			res[i][j] = 255*(res[i][j]-minval)/(maxval-minval);
		}
	}	
}


void compute_residuals_nonnormalized (double **res, double **imnoisy, double **imdenoised, int H, int W)
{
	int i,j;

	for(i=0;i<H;i++)
	{
		for(j=0;j<W;j++)
		{
			res[i][j] = imnoisy[i][j]-imdenoised[i][j];
		}
	}
}

