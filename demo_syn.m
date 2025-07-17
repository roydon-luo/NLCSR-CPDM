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

clear
clc
close all
addpath utils
addpath utils\functions
dbstop if error

saveornot = 0;
dataset_name = 'Monno'; % others: Qiu, Wen

%% load images
imagedir = dir(fullfile('data', 'synthetic', dataset_name));
imagedir = imagedir(3:end,:);

for nn = 1
    image_name = imagedir(nn).name;
    im0   = im2double(imread(fullfile('data', 'synthetic', dataset_name, image_name ,'0.png')));
    im45  = im2double(imread(fullfile('data', 'synthetic', dataset_name, image_name ,'45.png')));
    im90  = im2double(imread(fullfile('data', 'synthetic', dataset_name, image_name ,'90.png')));
    im135 = im2double(imread(fullfile('data', 'synthetic', dataset_name, image_name ,'135.png')));
    S_ori = cat(3,im0,im45,im90,im135);
    S_ori = S_ori(301:516,401:616,:); % for quicker test
    Mask = get_cpdmask(size(S_ori,1),size(S_ori,2),'rggb');
    CPFA = sum(S_ori.*Mask,3);
    [S_ini,Mask] = init_interp(CPFA);
    para = para_set(S_ini);
    para.Smask = CPFA.*Mask;
    para.C = 0.01;
    para.C_min = 1e-5; % lower than in processing real-world image due to it has fewer noises

    %% main iteration
    S = S_ini;
    D0 = init_dict(para);
    D = D0;
    
    for itr = 1:para.main_iternum
        [S_pat,para] = image2patch(S,para);
        L = size(S_pat,4);
        if L >= 1000
            patch_num = para.patch_num;
        else
            patch_num = L;
        end
        n_patches = ceil(L/patch_num);
        D = repmat(D,[1,1,1,n_patches]);
        DX = zeros(size(S_pat));
        for n = 1:n_patches
            idx_start = (n-1)*patch_num + 1;
            idx_end   = min(n*patch_num, L);
            [DX_cur, D_cur] = NLCSR(D(:,:,:,n),S_pat(:,:,:,idx_start:idx_end),para);
            DX(:,:,:,idx_start:idx_end) = DX_cur;
            D(:,:,:,n) = D_cur;
        end
        S_dx  = patch2image(DX,para);
        S_rec = S_update(S_dx,Mask,para);
        delta_S = init_interp(CPFA-sum(S_rec.*Mask,3));
        delta_S(delta_S >= 0.2) = 0;
        S_rec = S_rec + delta_S;
    end
 
    %% visualization
    [imS0,imDoLP,imAoLP] = cal_Stokes(S_ori);
    [iniS0,iniDoLP,iniAoLP] = cal_Stokes(S_ini);
    [imOutS0,imOutDoLP,imOutAoLP] = cal_Stokes(S_rec);
    imDoLP_str = contrast_str(imDoLP,0.3);
    iniDoLP_str = contrast_str(iniDoLP,0.3);
    imOutDoLP_str = contrast_str(imOutDoLP,0.3);
    figure,
    subplot(1,3,1);imshow(imS0,[]);title('GT');
    subplot(1,3,2);imshow(iniS0,[]);title('initial S0');
    subplot(1,3,3);imshow(imOutS0,[]);title('output S0');
    figure,
    subplot(1,3,1);imshow(colorjetmap(imDoLP_str),[]);title('GT');
    subplot(1,3,2);imshow(colorjetmap(iniDoLP_str),[]);title('initial DoLP');
    subplot(1,3,3);imshow(colorjetmap(imOutDoLP_str),[]);title('output DoLP');
    figure,
    subplot(1,3,1);imshow(colorjetmap((imAoLP+pi/2)/pi),[]);title('GT');
    subplot(1,3,2);imshow(colorjetmap((iniAoLP+pi/2)/pi),[]);title('initial AoLP');
    subplot(1,3,3);imshow(colorjetmap((imOutAoLP+pi/2)/pi),[]);title('output AoLP');

    %% measurement
    res_ini_S = cal_RMSE_PSNR_SSIM(S_ini, S_ori, 5, 5, 1);
    res_ini_S0 = cal_RMSE_PSNR_SSIM(iniS0, imS0, 5, 5, 1);
    res_ini_DoLP = cal_RMSE_PSNR_SSIM(iniDoLP, imDoLP, 5, 5, 1);
    res_ini_AoLP = cal_mae(iniAoLP, imAoLP);
    res_S = cal_RMSE_PSNR_SSIM(S_rec, S_ori, 5, 5, 1);
    res_S0 = cal_RMSE_PSNR_SSIM(imOutS0, imS0, 5, 5, 1);
    res_DoLP = cal_RMSE_PSNR_SSIM(imOutDoLP, imDoLP, 5, 5, 1);
    res_AoLP = cal_mae(imOutAoLP, imAoLP);
    
    %% saving
    if saveornot == 1
        file_name = 'synthetic_results';
        [~,save_name,~] = fileparts(imagedir(nn).name);
        if ~exist(file_name,'dir')
            mkdir(file_name)
        end
        imwrite(imOut90,sprintf([file_name '/' image_name '_90.png']));
        imwrite(imOut45,sprintf([file_name '/' image_name '_45.png']));
        imwrite(imOut135,sprintf([file_name '/' image_name '_135.png']));
        imwrite(imOut0,sprintf([file_name '/' image_name '_0.png']));
        imwrite(imOutS0,sprintf([file_name '/' image_name '_S0.png']));
        imwrite(imOutDoLP_str,sprintf([file_name '/' image_name '_DoLP.png']));
        imwrite(imOutAoLP,sprintf([file_name '/' image_name '_AoLP.png']));
    end
end