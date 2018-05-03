function [V0, K, sig1, sig2] = eegsim(nt , St, snr)
% EEGSIM is to generate EEG potentials based on the SNR uniform distribution model
% sig2s and sig1s are the square of sig2 and sig1
% nt: number of time points
% St: Samle total variance across all the electrodes St=ne*sig1s+sig2s*trace(K*K')
% snr: recording pure signal to recording noise
% see also RSEEGSIM

% Andy Hu, 31/08/2017

if nargin==0
    nt=1000;
    St=2000;
    snr=5;
end;

% lead field
iK = importdata('int-MC0000005-CHBM-EEG.mat');  % loading Leadfield...
indices = importdata('indices-red-6000.mat');
K = iK(1:19,indices); % loading indices for 6K reduction...
[ne, ns]=size(K);

% parameters
sig1s=St / ne / (1+snr);
sig2s=St / trace(K*K') * (snr/(1+snr));

sig1=sqrt(sig1s);
sig2=sqrt(sig2s);
% forward
mu=10*randn(1,ns);
j=mvnrnd(mu,sig2s*eye(ns),nt)';

V0=K*j;

Noise=sqrt(sig1s)*randn(size(V0));
V1 = V0 + Noise;

figure,imagesc(cov(V1')); colorbar;

V0=V1;
end