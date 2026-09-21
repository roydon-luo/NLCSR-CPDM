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

function B = get_simpara(S_pat, Sim_pat, Pat_arr, para)
% L2 parameters for NLSS constraint
f    = para.f;
ch   = para.ch;
L    = size(S_pat,4);
ns   = para.ns;
nSig = para.nSig;

if para.gpu_num > 0
    Cu0  = zeros(f,f,ch,L,'gpuArray');
    b0   = zeros(f,f,ch,L,'gpuArray');
else
    Cu0  = zeros(f,f,ch,L);
    b0   = zeros(f,f,ch,L);
end

for i = 1:L
    coe = S_pat(:,:,:,Pat_arr(:,i)) - repmat(Sim_pat(:,:,:,i),[1,1,1,ns]);
    Cu0(:,:,:,i) = mean(coe.^2,4);
    b0(:,:,:,i)  = S_pat(:,:,:,i) - Sim_pat(:,:,:,i);
end
Cu0 = max(0, Cu0-nSig^2);
b0  = (Cu0.*(b0.^2)).^(1/4);
B   = b0.^2 + para.Eps;
end
