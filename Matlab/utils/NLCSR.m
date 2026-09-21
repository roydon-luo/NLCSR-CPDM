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

function [DX, D] = NLCSR(D0, S_pat, para)
%% load parameters
maxiter = para.maxiter;
eAbs = para.eAbs;
eRel = para.eRel;
ch = para.ch;
m = para.m; K = para.K;
mu = para.mu; itar = para.itar; % varying rho/sigma/omega
vec = @(x) x(:);
% cdl (d update)
alphad = para.relaxParam_d;
SigUpdateCycle = para.SigmaUpdateCycle;
cdl_iters = para.cdl_iters;
rho = para.rho;
% csc (x update)
f = para.f;
C = para.C;
alphax = para.relaxParam_x;
lamb = para.lamb;
RhoUpdateCycle = para.RhoUpdateCycle;
csc_iters = para.csc_iters;
sig = para.sigma;

%% GPU
gpu_num = para.gpu_num;
if gpu_num > 0
    D0 = gpuArray(D0);
    S_pat = gpuArray(S_pat);
end

%% find similar patches
L = size(S_pat,4);
[Sim_pat, Pat_arr] = get_simblk(S_pat, para);
B = reshape(get_simpara(S_pat, Sim_pat, Pat_arr, para),f,f,ch,1,L);

%% optimal iteration
D = padarray(reshape(D0,m,m,1,K),[f-m,f-m],'post');
S_pat = reshape(S_pat,f,f,ch,1,L);
Sim_pat = reshape(Sim_pat,f,f,ch,1,L);
X = zeros(f,f,ch,K,L);
U = zeros(f,f,ch,K,L); % residual for X
V = zeros(f,f,1,K,L);  % residual for D
Nx = numel(X);
Nd = numel(D);

% D and X updating
itr_dx = 0;
eprix = 0; eduax = 0; rx = inf; sx = inf;
eprid = 0; eduad = 0; rd = inf; sd = inf;
while itr_dx <= maxiter && (rx > eprix || sx > eduax || rd > eprid || sd > eduad)
    itr_dx = itr_dx + 1;
    eta = C./B;
    if rem(itr_dx,1) == 0 && C >= para.C_min
        C = C*0.92;
    end
    % X update
    for tx = 1:csc_iters
        Xprv = X;
        Z = Z_update(X-U,D,S_pat,Sim_pat,rho,eta,para);
        Zr = alphax*Z + (1-alphax)*X;
        X = sfthrsh(Zr+U, lamb/rho);
        U = Zr - X + U;
    end
    % residuals X
    nX = norm(X(:)); nZ = norm(Z(:)); nU = norm(U(:));
    rx = norm(vec(X-Z));    % primal residual
    sx = rho*norm(vec(Xprv-X)); % dual residual
    eprix = sqrt(Nx)*eAbs+max(nX,nZ)*eRel;
    eduax = sqrt(Nx)*eAbs+rho*nU*eRel;
    if para.AutoRho && rem(itr_dx,RhoUpdateCycle) == 0
        [rho, U] = par_update(rho,rx,sx,mu,itar,U);
    end
    % D update
    for td = 1:cdl_iters
        Dprv = D;
        G = G_update(X,D-V,S_pat,Sim_pat,sig,eta,para);
        Gr = alphad*G + (1-alphad)*D;
        D = D_proj(sum(Gr+V,5)/L);
        V = Gr - D + V;
    end
    % residual D
    nG = norm(G(:)); nD = norm(D(:))*sqrt(L); nV = norm(V(:));
    rd = norm(vec(G-D));
    sd = sig*norm(vec(Dprv-D));
    eprid = sqrt(Nd)*eAbs+max(nD,nG)*eRel;
    eduad = sqrt(Nd)*eAbs+sig*nV*eRel;
    if para.AutoSigma && rem(itr_dx,SigUpdateCycle) == 0
        [sig,V] = par_update(sig,rd,sd,mu,itar,V);
    end
end

if gpu_num > 0
    D = gather(D);
    X = gather(X);
    DX = ifft2(fft2(D,f,f).*fft2(X),'symmetric');
    DX = gpuArray(DX);
    DX = gather(reshape(sum(DX,4),f,f,ch,size(DX,5)));
    D = gather(reshape(D(1:m,1:m,:),[m,m,K]));
else
    DX = ifft2(fft2(D,f,f).*fft2(X),'symmetric');
    DX = reshape(sum(DX,4),f,f,ch,size(DX,5));
    D = reshape(D(1:m,1:m,:),[m,m,K]);
end
end

function y = sfthrsh(x, kappa)
y = sign(x).*max(0,abs(x)-kappa);
end

function D = D_proj(D)
D = D./max(sqrt(sum(D.^2,1:2)),0);
end

function [par,Y] = par_update(par,r,s,mu,itar,Y)
a = 1;
if r > mu*s
    a = itar;
end
if s > mu*r
    a = 1/itar;
end
par_ = a*par;
if par_>1e-4
    par = par_;
    Y = Y/a;
end
end