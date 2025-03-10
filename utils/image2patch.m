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

function [S_blk,para] = image2patch(S,para)
[h,w,ch] = size(S);
f = para.f;
s = para.step;
H = h-f+1;
W = w-f+1;
r = 1:s:H; r = [r r(end)+1:H]; para.r = r;
c = 1:s:W; c = [c c(end)+1:W]; para.c = c;
H1 = length(r);
W1 = length(c);
L1 = H1*W1;

S_blk = zeros(f,f,ch,L1);
a = 0;
for j = c
    for i = r
        a = a + 1;
        S_blk(:,:,:,a) = S(i:i+f-1,j:j+f-1,:);
    end
end