%
%Code for submission to TPAMI: Image Denoising using the Higher Order Singular Value Decomposition, (version %0.0.1):
%---------------------------------------------------
%Copyright (C) 2011 Ajit Rajwade, Anand Rangarajan and Arunava Banerjee
% Authors: Ajit Rajwade, Anand Rangarajan and Arunava Banerjee
% Date:    Jan 21st 2011
%
% Contact Information:
%
%Ajit Rajwade:	avr@cise.ufl.edu
%Anand Rangarajan: anand@cise.ufl.edu
%Arunava Banerjee: arunava@cise.ufl.edu
% Terms:	  
%
%The source code is provided under the
%terms of the GNU General Public License (version 3).
%

%
%   This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>. */
%

clear;

% NOTE: PRIOR TO running this file, download ssim_index.m from the following link:
% https://ece.uwaterloo.ca/~z70wang/research/ssim/ssim_index.m



% ORACLE SVD SCRIPT
ippath = './GrayScale/OrigImages';

ps = 8;
choice = 1; % choice = 1 (PSNR), 2 (SSIM)

for sigma = 10:5:35
    D = 3*ps*ps*sigma*sigma;    
    files = dir([ippath '/*.png']);
    n = length(files);
    fprintf ('\n\nSIGMA = %d',sigma);
    
    for i=1:n
        filename = sprintf ('%s/%s',ippath,files(i).name);
        im = double((imread(filename))); 
        im = im(:,:,1);
        [H,W] = size(im);
        
        fname = files(i).name; fname(length(fname)-3:length(fname)) = []; 
        fname2 = [fname '.pgm'];
        
        % our method
        filename = sprintf ('./CCode/smoothed_P%dx%d_Sig%.2f_D%.2f_iter0_%s.pgm',ps,ps,sigma,D,fname);
        im2 = double(imread(filename));
        
        MSE2 = sum((im(:)-im2(:)).^2)/(H*W);
        PSNR_our = 10*log10(255*255/MSE2);        
        SSIM_our = ssim_index(im,im2);
        MSSIM_our = ssim(im,im2);
        if choice == 1, fprintf ('\n%s & %.3f',fname,PSNR_our);  else fprintf ('\n%s & %.3f',fname,SSIM_our); end
             
        % HOSVD method
        %hosvd_P8x8_Sig25.00_D120000.00_iter0_airplane.pgm
        filename = sprintf ('./CCode/hosvd_P%dx%d_Sig%.2f_D%.2f_iter0_%s.pgm',ps,ps,sigma,D,fname);
        im2 = double(imread(filename));        
        MSE2 = sum((im(:)-im2(:)).^2)/(H*W);
        PSNR_HOSVD = 10*log10(255*255/MSE2);        
        SSIM_HOSVD = ssim_index(im,im2);
        MSSIM_HOSVD = ssim(im,im2);
        if choice == 1, fprintf (' & %.3f',PSNR_HOSVD);  end;
        if choice == 2, fprintf (' & %.3f',SSIM_HOSVD);  end
        if choice == 3, fprintf (' & %.3f',MSSIM_HOSVD);  end
        
         
         % 3D-DCT (our implementation)
         % DCTdynamicK_P8x8_Sig25.00_D120000.00_iter0_airplane.pgm
         filename = sprintf ('./CCode/DCTdynamicK_P%dx%d_Sig%.2f_D%.2f_iter0_%s.pgm',ps,ps,sigma,D,fname);
         im2 = double(imread(filename));        
         MSE2 = sum((im(:)-im2(:)).^2)/(H*W);
         PSNR_3DDCT = 10*log10(255*255/MSE2);        
         SSIM_3DDCT = ssim_index(im,im2);
         MSSIM_3DDCT = ssim(im,im2);
         if choice == 1, fprintf (' & %.3f',PSNR_3DDCT); end
         if choice == 2, fprintf (' & %.3f',SSIM_3DDCT); end 
         if choice == 3, fprintf (' & %.3f',MSSIM_3DDCT); end 
        
        fprintf ('\\\\');
    end
end




