% ========================================================================
% Learning a Non-Locally Regularized Convolutional Sparse Representation
% for Joint Chromatic and Polarimetric Demosaicking, Version 1.0
% Copyright(c) 2024 Yidong Luo, Junchao Zhang, Jianbo Shao, Jiandong Tian,
% Jiayi Ma, All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%
%----------------------------------------------------------------------
% This is an demo of "Learning a Non-Locally Regularized Convolutional 
% Sparse Representation for Joint Chromatic and Polarimetric Demosaicking."
%
% Please cite the following paper when you use it:
%
% Yidong Luo, Junchao Zhang, Jianbo Shao, Jiandong Tian, Jiayi Ma, "Learn- 
% -ing a Non-Locally Regularized Convolutional Sparse Representation for 
% Joint Chromatic and Polarimetric Demosaicking." IEEE Trans. on Image 
% Processing, 33, 5029-5044 (2024).
%----------------------------------------------------------------------

function S_rec = patch2image(S_recblk,para)
f = para.f;
ch = para.ch;
h = para.h;
w = para.w;

blk = zeros(h,w,ch);
cnt = zeros(h,w,ch);
nn = 0;
for j = para.c
    for i = para.r
        nn = nn + 1;
        temp = S_recblk(:,:,:,nn);
        blk(i:i+f-1,j:j+f-1,:) = blk(i:i+f-1,j:j+f-1,:) + temp;
        cnt(i:i+f-1,j:j+f-1,:) = cnt(i:i+f-1,j:j+f-1,:) + ones(f,f,ch);
    end
end

S_rec = double(blk./cnt);
            