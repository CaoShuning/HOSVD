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

void readPGM_Dump (double **im, char *fname, int height, int width)
{
	FILE *fp;
	int i,j,k;
	double c;

	fp = fopen (fname,"rb");
	for (i=0;i<height;i++)
	{
		for(j=0;j<width;j++)
		{
			k = fscanf (fp,"%lf",&c);
			im[i][j] = (double)c;
		}
	}
	fclose (fp);
	printf ("\nFinished reading %s (dump)",fname); fflush (stdout);
}

double ** readPGM (char *fname, int *height, int *width)
{
	FILE *fp;
	char currstring[500],*p;
	int i,j,k,maxval;
	unsigned int c;
	double **im;
	char *cc;

	fp = fopen (fname,"rb");
	cc = fgets(currstring,100,fp); /* magic number P2 */
	do
	{
		cc = fgets(currstring,100,fp);
	} while (currstring[0] == '#');

	p = strtok(currstring," \t\n");
	*width = atoi(p); /* height and width */
	p = strtok(NULL," \t\n");
	*height = atoi(p);
	im = allocate_2d_double(*height,*width,'0');

	cc = fgets(currstring,100,fp);
	p = strtok(currstring," \t\n");
	maxval = atoi(p); /* height and width */

	if (VERBOSE) {printf ("\nHeight = %d, width = %d, max value = %d",*height,*width,maxval); fflush(stdout);}
	for (i=0;i<*height;i++)
	{
		for(j=0;j<*width;j++)
		{
			k = fscanf (fp,"%d",&c);
			im[i][j] = (double)c;
		}
	}

	fclose (fp);
	printf ("\nFinished reading %s",fname); fflush (stdout);

	return im;
}

void writePGM (char *fname, double **im, int height, int width)
{
	FILE *fp;
	int i,j,count=0;
	double a;

	fp = fopen (fname,"wb");

	fprintf (fp,"P2\n");
	fprintf (fp,"%d %d\n255\n",width,height);

	if (VERBOSE) {printf  ("\nWriting to %s of height %d and width %d",fname,height,width);}
	fflush (stdout);

	for (i=0;i<height;i++)
	{
		for(j=0;j<width;j++)
		{
			a = im[i][j]; if (a < 0) a = 0; if (a > 255.0) a = 255.0;
			fprintf (fp,"%d ",(unsigned int)a);
			//printf ("%.2lf ",im[i][j]);
			count++;
			if (count % 20 == 0) fprintf (fp,"\n");
		}
	}

	fclose (fp);
}

