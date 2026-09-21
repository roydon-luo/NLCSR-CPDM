function [S_ini,mask] = init_interp(CPFA)
h = size(CPFA,1);
w = size(CPFA,2);
pattern    = 'rggb';
col_method = "RI";
pol_method = "PCDP";
mask = get_cpdmask(h,w,pattern);

% color demosaicking
Bayer_90  = CPFA(1:2:end,1:2:end);
Bayer_45  = CPFA(1:2:end,2:2:end);
Bayer_135 = CPFA(2:2:end,1:2:end);
Bayer_0   = CPFA(2:2:end,2:2:end);

if col_method == "RI"
    sigma = 1; 
    eps = 1e-32;
    BayerDem_90  = bayer_residual(repmat(Bayer_90,[1,1,3]),pattern,sigma,eps);
    BayerDem_45  = bayer_residual(repmat(Bayer_45,[1,1,3]),pattern,sigma,eps);
    BayerDem_135 = bayer_residual(repmat(Bayer_135,[1,1,3]),pattern,sigma,eps);
    BayerDem_0   = bayer_residual(repmat(Bayer_0,[1,1,3]),pattern,sigma,eps);  
elseif col_method == "BI"
    eps = 1e-32;
    mask_bayer(1:2:h/2, 1:2:w/2, 1) = 1;
    mask_bayer(1:2:h/2, 2:2:w/2, 2) = 1;
    mask_bayer(2:2:h/2, 1:2:w/2, 2) = 1;
    mask_bayer(2:2:h/2, 2:2:w/2, 3) = 1;
    
    mosaic_bayer90 = Bayer_90.*mask_bayer;
    mosaic_bayer45 = Bayer_45.*mask_bayer;
    mosaic_bayer135 = Bayer_135.*mask_bayer;
    mosaic_bayer0 = Bayer_0.*mask_bayer;
    
    BayerDem_90  = bayer_bilinear(mosaic_bayer90);
    BayerDem_45  = bayer_bilinear(mosaic_bayer45);
    BayerDem_135 = bayer_bilinear(mosaic_bayer135);
    BayerDem_0   = bayer_bilinear(mosaic_bayer0);   
end

% polarization demosaicking
BayerDem_RGB = zeros(h,w,3);
BayerDem_RGB(1:2:end,1:2:end,:) = BayerDem_90;
BayerDem_RGB(1:2:end,2:2:end,:) = BayerDem_45;
BayerDem_RGB(2:2:end,1:2:end,:) = BayerDem_135;
BayerDem_RGB(2:2:end,2:2:end,:) = BayerDem_0;  

if pol_method == "PCDP"
    mask_dofp = zeros(h,w,4);
    mask_dofp(1:2:h, 1:2:w, 1) = 1;
    mask_dofp(1:2:h, 2:2:w, 2) = 1;
    mask_dofp(2:2:h, 1:2:w, 4) = 1;
    mask_dofp(2:2:h, 2:2:w, 3) = 1;
    for k = 1:3
        mosaic = padarray(BayerDem_RGB(:,:,k).*mask_dofp,[8,8,0],'symmetric');
        mask_dofp_in = padarray(mask_dofp,[8,8,0],'symmetric');
        [i90, i45, i0, i135] = PCDP(mosaic,mask_dofp_in);
        Dem_0(:,:,k)   = i0; 
        Dem_45(:,:,k)  = i45; 
        Dem_90(:,:,k)  = i90; 
        Dem_135(:,:,k) = i135;
    end
elseif pol_method == "EARI"
    mask_p90  = mask(:,:,7) + mask(:,:,8) + mask(:,:,9);
    mask_p45  = mask(:,:,4) + mask(:,:,5) + mask(:,:,6);
    mask_p135 = mask(:,:,10) + mask(:,:,11) + mask(:,:,12);
    mask_p0   = mask(:,:,1) + mask(:,:,2) + mask(:,:,3);
    [Dem_0, Dem_45, Dem_90, Dem_135] = EARI(BayerDem_RGB,eps,mask_p0,mask_p45,mask_p90,mask_p135);
end

S_ini = cat(3,Dem_0,Dem_45,Dem_90,Dem_135);
end
