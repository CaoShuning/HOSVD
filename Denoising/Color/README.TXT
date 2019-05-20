WHAT DO THESE FILES DO?

generateGaussianNoise.c : Generate a corrupted version of a clean underlying image
denoise_color_dynamicK.c : NL-SVD algorithm on color images in the Y-cb-cr representation
denoise_color_dynamicK_KLT.c : NL-SVD algorithm on color images in the KLT representation (learned from the noisy pixels)
denoise_color_hosvd_dynamicK_KLT.c: 3D-IHOSVD algorithm on color images with the KLT representation (learned from the noisy pixels)
color_hosvd_denoise.m: This code implements the 4D-HOSVD algorithm. 

TO COMPILE C CODE, GSL (i.e. the GNU scientific library) must be installed on the machine:
chmod +x makefile
./makefile

TO RUN THE MATLAB CODE:
To run the files color_hosvd_denoise.m (which executes the 4D-HOSVD algorithm) and color_hosvd_denoise_wiener.m (which executes 4D-HOSVD with a Wiener filter), you need to download and install the Tensor Algebra Library
from here below:
http://tptool.sztaki.hu/uploads/Main/tptool_090831.zip
Specifically, we need two routines 'hosvd.m' and 'tprod.m' from the subdirectory /tptools/array which is created when you unzip the downloaded zip file.
Make sure this subdirectory is included in your MATLAB path before you run the files color_hosvd_denoise.m or
color_hosvd_denoise_wiener.m .

The routine is to be run as follows:
[im2,psnr] = color_hosvd_denoise (im1,im,sigma,ps)
where im2 = output image, psnr = output PSNR, im1 = noisy image, im = clean image, sigma = std. dev. of the zero-mean Gaussian noise, ps = patchsize.

AND

[im2,im3,psnr,psnr2] = color_hosvd_denoise_wiener (im1,im,sigma,ps)
where im2 and im3 are the two output images (the rest of the output variables and parameters are the same)

INPUT FORMAT:
Read the README.txt file from the GrayScale directory. The format of execution is very similar. 
The one major difference is that all inputs to the color image denoising programs in this folder are in ASCII .ppm (as opposed to .pgm) format. 
The outputs are also produced in the .ppm format.


DATABASE FOR TESTING:
Download data from the Kodak gallery online. For further instructions refer to the file
./OrigColorImages/source_of_images.txt
