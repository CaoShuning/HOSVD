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

float **allocate_2d_float(int N,int M,char zero)
{
	int i;
	float **mymat;
	
	mymat=(float **)malloc(N*sizeof(float *));
	check_ptr(mymat,"allocate_2d_float");
	if(!zero)
		for(i=0;i<N;i++) {
			mymat[i]=(float *)malloc(M*sizeof(float));
			check_ptr(mymat[i],"allocate_2d_float");
	}
	else
		for(i=0;i<N;i++) {
			mymat[i]=(float *)calloc(M,sizeof(float));
			check_ptr(mymat[i],"allocate_2d_float");
	}
	return(mymat);
}

void free_2d_float(float **a,int N)
{
	int i;
	
	for(i=0;i<N;i++)
		free((void *)a[i]); 
	free((void *)a);
}


char **allocate_2d_char(int N,int M,char zero)
{
	int i;
	char **mymat;
	
	mymat=(char **)malloc(N*sizeof(char *));
	check_ptr(mymat,"allocate_2d_float");
	if(!zero)
		for(i=0;i<N;i++) {
			mymat[i]=(char *)malloc(M*sizeof(char));
			check_ptr(mymat[i],"allocate_2d_float");
	}
	else
		for(i=0;i<N;i++) {
			mymat[i]=(char *)calloc(M,sizeof(char));
			check_ptr(mymat[i],"allocate_2d_float");
	}
	return(mymat);
}




double **allocate_2d_double(int N,int M,char zero)
{
	int i;
	double **mymat;
	
	mymat=(double **)malloc(N*sizeof(double *));
	check_ptr(mymat,"allocate_2d_double");
	if(!zero)
		for(i=0;i<N;i++) {
			mymat[i]=(double *)malloc(M*sizeof(double));
			check_ptr(mymat[i],"allocate_2d_double");
		}
	else
		for(i=0;i<N;i++) {
			mymat[i]=(double *)calloc(M,sizeof(double));
			check_ptr(mymat[i],"allocate_2d_double");
		}
	return(mymat);
}


double ***allocate_3d_double(int N1, int N2, int N3, char zero)
{
	int i;
	double ***mymat;

	mymat=(double ***)malloc(N1*sizeof(double **));
	check_ptr(mymat,"allocate_3d_double");	

	for(i=0;i<N1;i++)
	{
		mymat[i] = allocate_2d_double(N2,N3,zero);
		check_ptr(mymat[i],"allocate_3d_double");
	}

	return mymat;
}

int ***allocate_3d_int(int N1, int N2, int N3, char zero)
{
	int i;
	int ***mymat;

	mymat= (int ***)malloc(N1*sizeof(int **));
	check_ptr(mymat,"allocate_3d_int");	

	for(i=0;i<N1;i++)
	{
		mymat[i] = allocate_2d_int(N2,N3,zero);
		check_ptr(mymat[i],"allocate_3d_int");
	}

	return mymat;
}


void free_2d_double(double **a,int N)

{
	int i;
	
	for(i=0;i<N;i++)
		free((void *)a[i]); 
	free((void *)a);
}


void free_3d_double (double ***a, int N1, int N2)
{
	int i;
	
	for(i=0;i<N1;i++)
		free_2d_double(a[i],N2);
	free((void *)a);	
}

void free_3d_int (int ***a, int N1, int N2)
{
	int i;
	
	for(i=0;i<N1;i++)
		free_2d_int(a[i],N2);
	free((void *)a);	
}


int **allocate_2d_int(int N,int M,char zero)

{
	int i;
	int **mymat;
	
	mymat=(int **)malloc(N*sizeof(int *));
	check_ptr(mymat,"allocate_2d_int");
	if(!zero)
		for(i=0;i<N;i++) {
			mymat[i]=(int *)malloc(M*sizeof(int));
			check_ptr(mymat[i],"allocate_2d_int");
	}
	else
		for(i=0;i<N;i++) {
			mymat[i]=(int *)calloc(M,sizeof(int));
			check_ptr(mymat[i],"allocate_2d_int");
	}
	return(mymat);
}

void free_2d_int(int **a,int N)

{
	int i;
	
	for(i=0;i<N;i++)
		free((void *)a[i]); 
	free((void *)a);
}

float *allocate_1d_float(int N,char zero)
{
	float *arr;
	
	if(!zero)
		arr=(float *)malloc(N*sizeof(float));
	else
		arr=(float *)calloc(N,sizeof(float));
	check_ptr(arr,"allocate_1d_float");
	return(arr);
}

double *allocate_1d_double(int N,char zero)
{
	double *arr;
	
	if(!zero)
		arr=(double *)malloc(N*sizeof(double));
	else
		arr=(double *)calloc(N,sizeof(double));
	check_ptr(arr,"allocate_1d_double");
	return(arr);
}

int *allocate_1d_int(int N,char zero)
{
	int *arr;
	
	if(!zero)
		arr=(int *)malloc(N*sizeof(int));
	else
		arr=(int *)calloc(N,sizeof(int));
	check_ptr(arr,"allocate_1d_int");
	return(arr);
}


unsigned char *allocate_1d_uchar(int N,char zero)
{
	unsigned char *arr;
	
	if(!zero)
		arr=(unsigned char *)malloc(N*sizeof(unsigned char));
	else
		arr=(unsigned char *)calloc(N,sizeof(unsigned char));
	check_ptr(arr,"allocate_1d_uchar");
	return(arr);
}

