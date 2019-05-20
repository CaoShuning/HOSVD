%
%Code for submission to TPAMI: Image Denoising using the Higher Order Singular Value Decomposition, (version %0.0.1):
%---------------------------------------------------
%Copyright (C) 2012 Ajit Rajwade, Anand Rangarajan and Arunava %Banerjee
% Authors: Ajit Rajwade, Anand Rangarajan and Arunava Banerjee
% Date:    June 4th 2012
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


for i=1:24
	if i < 10, filename = sprintf ('kodim0%d.ppm',i);
	else filename = sprintf ('kodim%d.ppm',i);
	end;

	im = double(imread(filename));
	filename2 = sprintf('colorhosvd_P8x8_Sig30.00_%s',filename);
	im2 = double(imread(filename2));
 
	[H,W] = size(im);
	MSE = sum((im(:)-im2(:)).^2)/(H*W);
	PSNR = 10*log10(255*255/MSE);
	fprintf ('\n%.3f',PSNR);

	im = double(imread(filename));
	filename2 = sprintf('colorhosvd2_P8x8_Sig30.00_%s',filename);
	im2 = double(imread(filename2));
 
	[H,W] = size(im);
	MSE = sum((im(:)-im2(:)).^2)/(H*W);
	PSNR = 10*log10(255*255/MSE);
	fprintf ('\t%.3f',PSNR);
end
