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

function [Sim_pat, Pat_arr] = get_simblk(S_pat, para)
f = para.f;
ch = para.ch;
L = size(S_pat,4);
ns = para.ns;
S_vec = reshape(S_pat,f*f*ch,L);
distances = pdist2(S_vec',S_vec','euclidean');

distances(logical(eye(size(distances)))) = inf;
[sorted_dis, sorted_ind] = sort(distances, 2);
min_dis = sorted_dis(:, 1:ns);
min_ind = sorted_ind(:, 1:ns);

Pat_arr = min_ind';
valid_min_dis = min_dis < para.min_dis;
valid_dis = valid_min_dis.*min_dis;
wei = exp(-valid_dis/10);
wei_arr = bsxfun(@rdivide, wei, sum(wei, 2) + eps);
wei_arr = wei_arr';

Sim_pat = zeros(f,f,ch,L);

for i = 1:ns
    wei = wei_arr(i,:);
    wei = reshape(wei,[1,1,1,L]);
    wei = repmat(wei,[f,f,ch,1]);
    sim_cur = S_pat(:,:,:,Pat_arr(i,:)).*wei;
    Sim_pat = Sim_pat + sim_cur;
end

