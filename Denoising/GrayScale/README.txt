What does each file do?

1) denoise_dynamicK.c : NL-SVD algorithm
2) denoise_dynamicK_nonrobust.c : NL-SVD algorithm without the robust covariance matrix using the hypothesis (KS) test
3) denoise_hosvd_dynamicK.c : HOSVD algorithm for denoising
4) denoise_DCT_dynamicK.c : our implementation of the BM3D algorithm using 3D-DCT (compare with the original BM3D paper)
5) generateNoise.c : a piece of code to generate noise from different noise models (using functions from GSL)
6) bestpatchsize_denoise_dynamicK.c : an NL-SVD implementation in which the best patch size is chosen using a criterion based on the residual image (correlations between patches in the residual)
7) bestsigma_denoise_dynamicK.c : an NL-SVD implementation in which the best sigma is chosen using the same criterion
8) denoise_hosvd_dynamicK_wiener.c : HOSVD and HOSVD2 algorithms for denoising
(HOSVD2 involves the subsequent Wiener filtering step)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
COMPILING THE CODE:

To compile these files, execute ./makefile at the command prompt. You must run this makefile on a machine which has the GSL (GNU SCIENTIFIC LIBRARY) installed. 

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
GENERATING NOISE INSTANCES:

./generateNoise.o <clean file name> <sigma> <noise model>
For noisemodel: (1) Gaussian, (2) Poisson, (3) Negative exponential. SIGMA is not relevant for Poisson model

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
EXECUTING THE FILES:

./<denoiser>.o <noisy file> <clean file> <SIGMA> <patchsize> <choiceDCT (1)/(0)
where '<denoiser>' is the name of the program.

where 
-<noisyfile> = name of the noisy file (in ASCII pgm format)
-<clean file> = name of the clean file (in ASCII pgm format)
-<SIGMA> = standard deviation of the noise (Gaussian noise)
-<patchsize> = width/height of the square patch
-<choiceDCT (1/0)> = use 0 if you want to initialize the procedure with a DCT-filtered output, and 1 othewise (in all results, we did not initialize the procedure with a DCT, i.e. we always chose 1 (or anything nonzero). Also, the procedure was always run for only one iteration).
-<MAXK - usually set to 30> - the maximum number of "similar patches" per reference patch
(usually set to 30) - NOTE: THIS PARAMETER IS FOR denoise_hosvd_dynamicK.o
AND denoise_hosvd_dynamicK_wiener.o ONLY.

Here is a sample run:
./denoise_hosvd_dynamicK.o <noisy file> <clean file> <SIGMA> <patchsize> <choiceDCT (1)/(0) <MAXK>
./denoise_hosvd_dynamicK_wiener.o <noisy file> <clean file> <SIGMA> <patchsize> <choiceDCT (1)/(0) <MAXK>

A typical sample run for HOSVD is as follows:
./denoise_hosvd_dynamicK.o  noisy_Sig20.00_barbara.pgm barbara.pgm 20 8 1 30
./denoise_hosvd_dynamicK_wiener.o  noisy_Sig20.00_barbara.pgm barbara.pgm 20 8 1 30

A typical sample run for NL-SVD is as follows:
./denoise_dynamicK.o  noisy_Sig20.00_barbara.pgm barbara.pgm 20 8 1

The output image will be stored in a file whose name will be printed on the screen. The format for the filename is as follows for
HOSVD/HOSVD2:
hosvd_P8x8_Sig20.00_D76800.00_K30_iter0_barbara.pgm and
hosvd2_P8x8_Sig20.00_D76800.00_K30_iter0_barbara.pgm

for patchsize 8 x 8, noise level 20, D = 3 X 8 X 8 X 20 X 20 = 76800, K = 30 and original file barbara.pgm.

Note that we supply the name of the clean file solely for the purpose of calculation of the PSNR or MSE (the true file is in no way used in the code :-)). This can be easily verified by purposefully supplying an incorrect file name. The PSNR values will no more make sense, but the output file will be the same).

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Note that though the program spits out SSIM values (in addition to PSNR), those figures are not valid as 
the SSIM implementation does not conform to the standards proposed by the authors of that quality metric. 
For comparison, execute the MATLAB script results_psnr_ssim_script.m (located in the parent directory) which in turn uses the SSIM implementation 
ssim_index.m (which is the original implementation by simoncelli and zhou wang). 
The comparison is between the output image files (stored on disk) in the pgm format. 



DATABASE:
For testing our algorithm, download images from the Lansel benchmark. Further details about this are included in the directory
./OrigImages/source_of_images.txt
