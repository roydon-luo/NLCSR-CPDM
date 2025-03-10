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

function [para] = para_set(S)
%% image
[h,w,ch] = size(S);
para.h   = h;       % height
para.w   = w;       % width
para.ch  = ch;      % channels
%% dictionary
para.m   = 5;       % dictionary size
para.K   = 16;      % number of filters     
%% patches and similar patches
para.f                = 6;     % patch size
para.ns               = 10;     % number of similar patches for each patch
para.step             = 3;      % step
para.nSig             = 0;
para.C                = 0.28;    % NLSS regularization
para.C_min            = 1e-4;
para.min_dis          = 0.1;
para.Eps              = 1e-6;
%% cdl (d update)
para.relaxParam_d     = 1.8;
para.sigma            = 1;
para.AutoSigma        = 1;
para.SigmaUpdateCycle = 1;
para.cdl_iters        = 1;
%% csc (x update)
para.relaxParam_x     = 1.8;
para.rho              = 1;
para.AutoRho          = 1;
para.RhoUpdateCycle   = 1;
para.csc_iters        = 1;
para.lamb             = 1e-4;  % L1 regularization
%% s update
para.relaxParam_s     = 1.8;
para.gamma            = 0.1;   % down-sampled regularization
para.tol              = 1e-5;
para.rho_hub          = repmat(0.05,[12,1]);
para.delta_hub        = repmat(0.001,[12,1]);
para.beta_hub         = repmat(0.95,[12,1]);
%% the general
para.main_iternum     = 1;
para.mu               = 5.0;
para.itar             = 1.2;
para.maxiter          = 20;
para.eAbs             = 1e-4;
para.eRel             = 1e-4;
para.patch_num        = 1000;  % the number of image patches processed at once
para.gpu_num          = gpuDeviceCount;
end