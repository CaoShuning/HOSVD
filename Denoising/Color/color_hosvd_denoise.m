%
%Code for submission to TPAMI: Image Denoising using the Higher Order Singular Value Decomposition, (version %0.0.1):
%---------------------------------------------------
%Copyright (C) 2012 Ajit Rajwade, Anand Rangarajan and Arunava Banerjee
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


function [im2,psnr] = color_hosvd_denoise (im1,im,sigma,ps)

[H,W,D] = size(im);
im2 = zeros(H,W,D);
numcount = zeros(H,W);

SR = 20;
dist_threshold = 3*3*sigma*sigma*ps*ps;
maxK = 30;

for i=104:H-ps+1
    fprintf (' %d',i);
    for j=49:W-ps+1
        refpatch = im1(i:i+ps-1,j:j+ps-1,:);
       
        sr_top = max([i-SR 1]);
        sr_left = max([j-SR 1]);
        sr_right = min([j+SR W-ps+1]);
        sr_bottom = min([i+SR H-ps+1]);
        
        count = 0;
        similarity_indices = zeros((2*SR+1)^2,2); distvals = similarity_indices(:,1);
        for i1=sr_top:sr_bottom
            for j1=sr_left:sr_right
                currpatch = im1(i1:i1+ps-1,j1:j1+ps-1,:);
                dist = sum((refpatch(:)-currpatch(:)).^2);
                if dist < dist_threshold, count = count+1; distvals(count) = dist; similarity_indices(count,:) = [i1 j1]; end;
            end
        end
        
        similarity_indices = similarity_indices(1:count,:);
        distvals = distvals(1:count);
        
        if count > maxK
            [~,sortedindices] = sort(distvals,'ascend');
            similarity_indices = similarity_indices(sortedindices(1:maxK),:);
            count = maxK;
        end
        
        A = zeros(ps,ps,D,count);
        for k=1:count
            yindex = similarity_indices(k,1);
            xindex = similarity_indices(k,2);
            A(:,:,:,k) = im1(yindex:yindex+ps-1,xindex:xindex+ps-1,:);
        end
        
        coeff_threshold = sigma*sqrt(2*log(ps*ps*count*3));
        [S U] = hosvd(A);
        S(abs(S(:)) < coeff_threshold) = 0;
        A = tprod(S,U);
        
        for k=1:count
            yindex = similarity_indices(k,1);
            xindex = similarity_indices(k,2);
            im2(yindex:yindex+ps-1,xindex:xindex+ps-1,:) = im2(yindex:yindex+ps-1,xindex:xindex+ps-1,:)+A(:,:,:,k);
            numcount(yindex:yindex+ps-1,xindex:xindex+ps-1) = numcount(yindex:yindex+ps-1,xindex:xindex+ps-1)+1;
        end        
    end
end

im2(:,:,1) = im2(:,:,1)./numcount(:,:);
im2(:,:,2) = im2(:,:,2)./numcount(:,:);
im2(:,:,3) = im2(:,:,3)./numcount(:,:);

mse = sum((im2(:)-im(:)).^2)/(H*W*D);
psnr = 10*log10(255*255/mse);

    
