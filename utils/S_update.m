% This code is the customized version from "Linear polarization demosaicking
% for monochrome and colour polarization focal plane arrays, in: Computer
% Graphics Forum, Vol. 40, Wiley Online Library, 2021, pp. 77–89"

function S = S_update(DX,Mask,para)
% parameters
h         = para.h;
w         = para.w;
ch        = para.ch;
beta_hub  = para.beta_hub;
rho_hub   = para.rho_hub;
delta_hub = para.delta_hub;
maxiter   = para.maxiter;
tol       = para.tol;
S_mask    = para.Smask;
gamma     = para.gamma;

% functions defining
norm2  = @(x) norm(x(:));
rho_hub_update = @(x) bsxfun(@times,x,permute(rho_hub,[3 2 1]));
huber  = @(x,delta,rho,lamb) ...
         (abs(x) <= (lamb+rho)/(rho+eps)*delta) .* (rho/(lamb+rho+eps)*x) + ...
         (x > (lamb+rho)/(rho+eps)*delta) .* (x-delta/(rho+eps)*lamb) + ...
         (x < -(lamb+rho)/(rho+eps)*delta) .* (x+delta/(rho+eps)*lamb) ;

% 2D convolutional operators defining
bond_condition = 'circular';
dx = [0 -1 1];
dy = [0 -1 1]';
dxx = conv2(dx,dx);
dyy = conv2(dy,dy);
dxy = conv2(dx,dy);

Dx    = @(x) imfilter(x,dx, bond_condition);
Dy    = @(x) imfilter(x,dy, bond_condition);
Dxx   = @(x) imfilter(x,dxx,bond_condition);
Dyy   = @(x) imfilter(x,dyy,bond_condition);
Dxy   = @(x) imfilter(x,dxy,bond_condition);
Diff  = @(v) cat(4, Dx(v), Dy(v), Dxx(v), Dxy(v), Dyy(v));

DxT   = @(x) imfilter(x, rot90(dx,2),  bond_condition);
DyT   = @(x) imfilter(x, rot90(dy,2),  bond_condition);
DxxT  = @(x) imfilter(x, rot90(dxx,2), bond_condition);
DyyT  = @(x) imfilter(x, rot90(dyy,2), bond_condition);
DxyT  = @(x) imfilter(x, rot90(dxy,2), bond_condition);
DiffT = @(w) DxT(w(:,:,:,1)) + DyT(w(:,:,:,2)) + ...
             DxxT(w(:,:,:,3)) + DxyT(w(:,:,:,4)) + DyyT(w(:,:,:,5));

s     = [0 0 1 0 0; 0 1 -7 1 0; 1 -7 20 -7 1; 0 1 -7 1 0; 0 0 1 0 0];
DTD   = @(v) imfilter(v,s,bond_condition);
b1    = DX + gamma*S_mask;
Fi    = @(x) (gamma*Mask+1).*x;
A     = @(v) Fi(v) + rho_hub_update(DTD(v));

% initialization
Z_hub = zeros(h,w,ch,5);
U_hub = zeros(h,w,ch,5);
obj_total = zeros(maxiter+1, 3);
obj_total(1,:) = Inf;

% S update
itr = 0;
while itr <= maxiter
    %
    itr = itr + 1;
    b2  = rho_hub_update(DiffT(Z_hub-U_hub));
    b   = b1 + b2;
    x   = pcg_ND(A,b,tol);
    d   = Diff(x);
    dU_hub = d + U_hub;
    for ch_hub = 1:ch
        Z_hub(:,:,ch_hub,:) = huber(dU_hub(:,:,ch_hub,:), delta_hub(ch_hub), rho_hub(ch_hub), beta_hub(ch_hub));
    end
    U_hub  = dU_hub - Z_hub;
    % objective function
    obj_data = 0.5*(norm2(x - b1)^2 + gamma*norm2(S_mask-Mask.*x)^2);
    obj_Hub  = 0;
    for ch_hub = 1:ch
        obj_Hub = obj_Hub + beta_hub(ch_hub) * objective_hub(d(:,:,ch_hub,:), delta_hub(ch_hub));
    end
    obj_total(itr+1, 1) = obj_data + obj_Hub;
    obj_total(itr+1, 2) = obj_data;
    obj_total(itr+1, 3) = obj_Hub;  
    rel_err = abs(obj_total(itr+1) - obj_total(itr));
    if rel_err < tol
        break
    else
        gamma = gamma*0.9;
    end
end
S = x;
end

function p = objective_hub(z, delta)
    p = sum(HuberLoss(z(:), delta));
end

function g_hub = HuberLoss(z,k)
g_hub = zeros(size(z));
if k ~= 0
    g_hub(abs(z)<=k) = (1/2)*(z(abs(z)<=k).^2);
    g_hub(abs(z)>k)  = k*((abs(z(abs(z) > k))) - (1/2)*k);
else
    g_hub  = abs(z);    
end
end