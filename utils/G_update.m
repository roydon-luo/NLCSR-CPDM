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

function G = G_update(X,W,S,Sim,sig,eta,para)
C = conj(fft2(X))./((1+eta).*sum(abs(fft2(X)).^2,4)+sig);
Rf = fft2(S) + eta.*fft2(Sim) - sum(fft2(X).*fft2(fft2(W)),4);
Gf = fft2(repmat(W,[1,1,para.ch,1,1])) + C.*Rf;
clear X W S Sim eta C Rf
if para.gpu_num > 0
    Gf = gather(Gf);
    G = ifft2(Gf,'symmetric');
    G = sum(G,3)/para.ch;
    G = gpuArray(G);
else
    G = ifft2(Gf,'symmetric');
    G = sum(G,3)/para.ch;
end
end