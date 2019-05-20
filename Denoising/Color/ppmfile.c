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

void readPPMBinary (char *fname, int *height, int *width, double **im[3])
{
	FILE *fp;
	char currstring[500],*p;
	int i,j,maxval;
	unsigned char c1,c2,c3;

	fp = fopen (fname,"rb");
	fgets(currstring,100,fp); /* magic number P3 */
	do
	{
		fgets(currstring,100,fp);
	} while (currstring[0] == '#');

	p = strtok(currstring," \t\n");
	*width = atoi(p); /* height and width */
	p = strtok(NULL," \t\n");
	*height = atoi(p);

	fgets(currstring,100,fp);
	p = strtok(currstring," \t\n");
	maxval = atoi(p); 

	if (VERBOSE) {printf ("\nHeight = %d, width = %d, max value = %d",*height,*width,maxval); fflush(stdout);}

	im[0] = allocate_2d_double(*height,*width,'0');
	im[1] = allocate_2d_double(*height,*width,'0');	
	im[2] = allocate_2d_double(*height,*width,'0');

	for (i=0;i<*height;i++)
	{
		for(j=0;j<*width;j++)
		{
			//fscanf (fp,"%d %d %d",&c1,&c2,&c3);
			fread(&c1,sizeof(unsigned char),1,fp);
			fread(&c2,sizeof(unsigned char),1,fp);
			fread(&c3,sizeof(unsigned char),1,fp);

			im[0][i][j] = (double)c1;
			im[1][i][j] = (double)c2;
			im[2][i][j] = (double)c3;
		}
	}

	fclose (fp);
	printf ("\nFinished reading %s",fname); fflush (stdout);
}



void readPPM (char *fname, int *height, int *width, double **im[3])
{
	FILE *fp;
	char currstring[500],*p;
	int i,j,maxval;
	unsigned int c1,c2,c3;

	fp = fopen (fname,"rb");
	fgets(currstring,100,fp); /* magic number P3 */
	do
	{
		fgets(currstring,100,fp);
	} while (currstring[0] == '#');

	p = strtok(currstring," \t\n");
	*width = atoi(p); /* height and width */
	p = strtok(NULL," \t\n");
	*height = atoi(p);

	fgets(currstring,100,fp);
	p = strtok(currstring," \t\n");
	maxval = atoi(p); /* height and width */

	if (VERBOSE) {printf ("\nHeight = %d, width = %d, max value = %d",*height,*width,maxval); fflush(stdout);}

	im[0] = allocate_2d_double(*height,*width,'0');
	im[1] = allocate_2d_double(*height,*width,'0');	
	im[2] = allocate_2d_double(*height,*width,'0');

	for (i=0;i<*height;i++)
	{
		for(j=0;j<*width;j++)
		{
			fscanf (fp,"%d %d %d",&c1,&c2,&c3);
			im[0][i][j] = (double)c1;
			im[1][i][j] = (double)c2;
			im[2][i][j] = (double)c3;
		}
	}

	fclose (fp);
	printf ("\nFinished reading %s",fname); fflush (stdout);
}

void writePPM (char *fname, double **im[3], int height, int width)
{
	FILE *fp;
	int i,j,k,count=0;
	double a;

	fp = fopen (fname,"wb");

	fprintf (fp,"P3\n");
	fprintf (fp,"%d %d\n255\n",width,height);

	if (VERBOSE) {printf  ("\nWriting to %s of height %d and width %d",fname,height,width);}
	fflush (stdout);

	for (i=0;i<height;i++)
	{
		for(j=0;j<width;j++)
		{
			for(k=0;k<3;k++)
			{
				a = im[k][i][j]; if (a < 0) a = 0; if (a > 255.0) a = 255.0;
				fprintf (fp,"%d ",(unsigned int)a);
				//printf ("%.2lf ",im[i][j]);
				count++;
			}
		}
		fprintf (fp,"\n");
	}

	fclose (fp);
}

void writePPMBinary (char *fname, double **im[3], int height, int width)
{
	FILE *fp;
	int i,j,k,count=0;
	double a;
	unsigned char c;

	fp = fopen (fname,"wb");

	fprintf (fp,"P6\n");
	fprintf (fp,"%d %d\n255\n",width,height);

	if (VERBOSE) {printf  ("\nWriting to %s of height %d and width %d",fname,height,width);}
	fflush (stdout);

	for (i=0;i<height;i++)
	{
		for(j=0;j<width;j++)
		{
			for(k=0;k<3;k++)
			{
				a = im[k][i][j]; if (a < 0) a = 0; if (a > 255.0) a = 255.0;
				c = (unsigned char)a;
				fwrite(&c,sizeof(unsigned char),1,fp);
			}
		}
	}

	fclose (fp);
}


