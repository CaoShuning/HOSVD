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
#define EPS1 0.001
#define EPS2 1.0e-8

int mycompare (const void *a, const void *b)
{
	double da,db;
	da = *((double * )a);
	db = *((double * )b);
	if (da > db) return 1;
	if (da < db) return -1;
	return 0;
}

double probks(double alam)
{
	int j;
	double a2,fac=2.0,sum=0.0,term,termbf=0.0;

	a2 = -2*alam*alam;
	for(j=1;j<=100;j++)
	{
		term = fac*exp(a2*j*j);
		sum += term;
		if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
		fac = -fac;
		termbf = fabs(term);
	}
	return 1.0;
}


void ksLaplacian (double *data, int n, double *d, double *prob, double invlambda)
{
	int j;
	double dt,en,ff,fn,fo=0.0,f1,f2,s;

	qsort(data,n,sizeof(double),mycompare);
	en = n;
	*d = 0.0;

	for(j=0;j<n;j++)
	{
		fn = (j+1)/en;
		if (data[j]>0) s = 1.0; else if (data[j] == 0) s = 0.0; else s = -1.0;
		ff = 0.5*(1+s*(1-exp(-fabs(data[j])*invlambda)));
		f1 = fabs(fo-ff); f2 = fabs(fn-ff); 
		if (f1 > f2) dt = f1; else dt = f2;
		if (dt > *d) *d = dt;
		fo = fn;
	}
	en = sqrt(en);
	*prob = probks((en+0.12+0.11/en)*(*d));
}



void ksGaussian (double *data, int n, double *d, double *prob, double sigma)
{
	int j;
	double dt,en,ff,fn,fo=0.0,f1,f2;

	qsort(data,n,sizeof(double),mycompare);
	en = n;
	*d = 0.0;

	for(j=0;j<n;j++)
	{
		fn = (j+1)/en;
		ff = 0.5*(1+erf(data[j]/(sqrt(2.0)*sigma)));
		f1 = fabs(fo-ff); f2 = fabs(fn-ff); 
		if (f1 > f2) dt = f1; else dt = f2;
		if (dt > *d) *d = dt;
		fo = fn;
	}
	en = sqrt(en);
	*prob = probks((en+0.12+0.11/en)*(*d));
}

// code from Numerical Recipes in C, page 624-626, section 14.3
void kstwo (double *data1, int n1, double *data2, int n2, double *d, double *prob)
{
	int j1=0,j2=0;
	double d1,d2,dt,en1,en2,en,fn1=0.0,fn2=0.0;

	qsort(data1,n1,sizeof(double),mycompare);
	qsort(data2,n2,sizeof(double),mycompare);
	en1 = n1; en2 = n2;
	*d = 0.0;

	while (j1 < n1 && j2 < n2)
	{
		d1 = data1[j1]; d2 = data2[j2];
		if (d1 <= d2) fn1 = j1++/en1;
		if (d2 <= d1) fn2 = j2++/en2;
		dt = fabs(fn2-fn1);
		if (dt > *d) *d = dt;
	}

	en = sqrt(en1*en2/(en1+en2));
	*prob = probks((en+0.12+0.11/en)*(*d));
}

// code from Numerical Recipes in C, page 624-626, section 14.3
void kstwo_nosort (double *data1, int n1, double *data2, int n2, double *d, double *prob)
{
	int j1=0,j2=0;
	double d1,d2,dt,en1,en2,en,fn1=0.0,fn2=0.0;

	en1 = n1; en2 = n2;
	*d = 0.0;

	while (j1 < n1 && j2 < n2)
	{
		d1 = data1[j1]; d2 = data2[j2];
		if (d1 <= d2) fn1 = j1++/en1;
		if (d2 <= d1) fn2 = j2++/en2;
		dt = fabs(fn2-fn1);
		if (dt > *d) *d = dt;
	}

	en = sqrt(en1*en2/(en1+en2));
	*prob = probks((en+0.12+0.11/en)*(*d));
}

void KSTest_Noiseness (double **res, int H, int W, int ps, double *Kval, double *logPval)
{
	int i,j,k,k1,k2,numpatches,count;
	double dval,pval,**patches;	

	numpatches = floor(H*W/(ps*ps));
	patches = allocate_2d_double(numpatches,ps*ps,'0');	

	k = 0;
	for(i=0;i<=H-ps;i=i+ps)
	{
		for(j=0;j<=W-ps;j=j+ps)
		{
			count = 0;
			for(k1=i;k1<i+ps;k1++)
			{
				for(k2=j;k2<j+ps;k2++)
				{
					patches[k][count++] = res[k1][k2];
				}
			}	
			qsort(patches[k],count,sizeof(double),mycompare);		
			k++;
		}
	}

	*Kval = 0;
	*logPval = 0;
	count = 0;
	numpatches = k;
	for(i=0;i<numpatches-1;i++)
	{
		for(j=i+1;j<numpatches;j++)
		{
			kstwo_nosort (patches[i],ps*ps,patches[j],ps*ps,&dval,&pval);
			*Kval = *Kval + dval;
			*logPval = *logPval + log(pval);
			count++;
		}
	}

	*Kval /= count;
	*logPval = -*logPval/count;

	//printf ("\nScale = %d, K = %lf, Pval = %lf",ps,*Kval,*logPval); fflush(stdout);
	free_2d_double(patches,numpatches);
}

void KSTest_Noiseness_multiscale (double **res, int H, int W, int lowest_scale, int biggest_scale, double *Kval, double *logPval)
{
	int s;
	double a,b;

	*Kval = *logPval = 0;
	for(s=lowest_scale;s<=biggest_scale;s++)
	{
		KSTest_Noiseness (res,H,W,s,&a,&b);
		*Kval = *Kval + a;
		*logPval = *logPval + b;
	}

	*Kval /= (biggest_scale-lowest_scale+1);
	*logPval /= (biggest_scale-lowest_scale+1);
}

double corrcoeff (double *data1, double *data2, int n)
{
	int i,j;
	double mu1=0,mu2=0,c=0,s1=0,s2=0;

	for(i=0;i<n;i++)
	{
		mu1 += data1[i];
		mu2 += data2[i];
	}
	mu1 /= n;
	mu2 /= n;

	for(i=0;i<n;i++)
	{
		s1 += (mu1-data1[i])*(mu1-data1[i]);
		s2 += (mu2-data2[i])*(mu2-data2[i]);
	}
	s1 /= n;
	s2 /= n;
	s1 = sqrt(s1);
	s2 = sqrt(s2);


	for(i=0;i<n;i++)
	{
		c += (data1[i]-mu1)*(data2[i]-mu2);
	}
	c /= (s1*s2*n);
	
	return c;
}

double corrcoeff_Noiseness (double **res, int H, int W, int ps)
{
	int i,j,k,k1,k2,numpatches,count;
	double c=0.0,**patches;	

	numpatches = floor(H*W/(ps*ps));
	patches = allocate_2d_double(numpatches,ps*ps,'0');	

	k = 0;
	for(i=0;i<=H-ps;i=i+ps)
	{
		for(j=0;j<=W-ps;j=j+ps)
		{
			count = 0;
			for(k1=i;k1<i+ps;k1++)
			{
				for(k2=j;k2<j+ps;k2++)
				{
					patches[k][count++] = res[k1][k2];
				}
			}	
			k++;
		}
	}

	count = 0;
	numpatches = k;
	for(i=0;i<numpatches-1;i++)
	{
		for(j=i+1;j<numpatches;j++)
		{
			c += fabs(corrcoeff(patches[i],patches[j],ps*ps));
			count++;
		}
	}
	free_2d_double(patches,k);

	c /= count;
	return c;
}

double corrcoeff_Noiseness_multiscale (double **res, int H, int W, int lowest_scale, int biggest_scale)
{
	int s;
	double c = 0.0;

	for(s=lowest_scale;s<=biggest_scale;s++)
	{
		c += corrcoeff_Noiseness (res,H,W,s);
	}
	c /= (biggest_scale-lowest_scale+1);

	return c;
}

