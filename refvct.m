% function [Vi, msi] = refvct(V1, K1, K2, K3, K4, Vi, msi, i)
% REFVCT to evaluate the effects of AR, rAR, REST, rREST with volume conduction tested
% Input
%      V1  ---- EEG potentials with noise
%      K1  ---- Individual lead field
%      K2  ---- Sparse Individual lead field
%      K3  ---- Avgerged Individual lead field by trace
%      K4  ---- Sphercial lead field based on yao 3-layers concentric spherical head model
%       i  ---- the index for lead field in 89 subjects database
%       j  ---- the index for sensor noise test (see mksnr.m)
% Output
%      Vi  ---- reconstructed EEG potentials at infinity defined in regrefsim.m
%      msi ---- model selection criteria infinity defined in regrefsim.m
% See also regrefsim.m refbys.m

% Andy Hu, Nov, 9, 2017

Nc = size(V1,1); Vr = Hsc(Nc)*V1;

[Vi(:,:,:,1),msi(:,:,1)] = refbys(Vr,[],'ar');
[Vi(:,:,:,2),msi(:,:,2)] = refbys(Vr,K1,'rt');
[Vi(:,:,:,3),msi(:,:,3)] = refbys(Vr,K2,'rt');
[Vi(:,:,:,4),msi(:,:,4)] = refbys(Vr,K3,'rt');
[Vi(:,:,:,5),msi(:,:,5)] = refbys(Vr,K4,'rt');

% end