function [V0,V1,K1, K2, K3, K4]=rseegsim(sjidx,pt,snr)
% RSEEGSIM is using multivariate autoregression model for EEG simulation
% Input: 
%           sjidx ---- subject index for the lead field
%           pt ---- different patch type to generate potentials
%           snr ---- signal to nosie variance ratio
%           mode ---- output different K for REST (see trsK)
% output:
%           V0 ---- pure signal
%           V1 ---- pure signal plus noise
%           lf   ---- lead fields
%           K1 ---- individual lead field
%           K2 ---- sparse individual  lead field
%           K3 ---- averaged lead field
%           K4 ---- Yao three layer spherical lead field
% See also eegarsim

% Andy Hu, July, 2017

kd = dir('.\subject\LF'); kd(1:2)=[];
iK = importdata(kd(sjidx).name); %  lead field 80k
indices = importdata('indices-red-6000.mat'); % 80k to 6k indices
ctx = importdata('MC0000005-cortex.mat'); vertices = ctx.vertices; % 80K vertices
K = iK(:,indices); vertices = vertices(indices,:); % downsampling
faces = importdata('faces-red-6000.mat');  % 6k faces

% Defining cross-spectrum parameter
xp.nfft = 512;         % window length
xp.nv = size(K,2);     % Number of vertices
xp.nseg=12;   % changed from 102
xp.nt = xp.nfft*xp.nseg;     % nfft*X time samples were generated

[Source,A] = genarsos(pt,vertices,faces);
[V0,V1, J] = eegarsim(K,Source,xp,A,snr);
% viscs(J,vertices,faces);% visualing cortical surface

% volume conduction model matching
[K1, K2] = trsK(kd, [], sjidx, J);
[~, ~, K3, K4] = trsK( kd, []);

V0=V0./norm(V0,'fro');
V1=V1./norm(V1,'fro');
end