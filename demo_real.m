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
dataset_name = 'real';

%% load images
imagedir = dir(fullfile('data', dataset_name, '*.png'));

for nn = 1
    image_name = imagedir(nn).name;
    CPFA = load_cpfa(fullfile('data', dataset_name, image_name));
    CPFA = CPFA(301:516,301:516); % for quicker test
    [S_ini,Mask] = init_interp(CPFA);
    para = para_set(S_ini);
    para.Smask = CPFA.*Mask;

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
    end
 
    %% visualization
    [iniS0,iniDoLP,iniAoLP] = cal_Stokes(S_ini);
    [imOutS0,imOutDoLP,imOutAoLP] = cal_Stokes(S_rec);
    iniDoLP_str = contrast_str(iniDoLP,0.3);
    imOutDoLP_str = contrast_str(imOutDoLP,0.3);
    figure,
    subplot(1,2,1);imshow(iniS0,[]);title('initial S0');
    subplot(1,2,2);imshow(imOutS0,[]);title('output S0');
    figure,
    subplot(1,2,1);imshow(colorjetmap(iniDoLP_str),[]);title('initial DoLP');
    subplot(1,2,2);imshow(colorjetmap(imOutDoLP_str),[]);title('output DoLP');
    figure,
    subplot(1,2,1);imshow(colorjetmap((iniAoLP+pi/2)/pi),[]);title('initial AoLP');
    subplot(1,2,2);imshow(colorjetmap((imOutAoLP+pi/2)/pi),[]);title('output AoLP');
    
    %% saving
    if saveornot == 1
        file_name = 'real_results';
        [~,save_name,~] = fileparts(imagedir(nn).name);
        if ~exist(file_name,'dir')
            mkdir(file_name)
        end
        imwrite(imOut90,sprintf([file_name '/' save_name '_90.png']));
        imwrite(imOut45,sprintf([file_name '/' save_name '_45.png']));
        imwrite(imOut135,sprintf([file_name '/' save_name '_135.png']));
        imwrite(imOut0,sprintf([file_name '/' save_name '_0.png']));
        imwrite(imOutS0,sprintf([file_name '/' save_name '_S0.png']));
        imwrite(imOutDoLP_str,sprintf([file_name '/' save_name '_DoLP.png']));
        imwrite(imOutAoLP,sprintf([file_name '/' save_name '_AoLP.png']));
    end
end