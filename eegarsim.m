function [V0,V1,J, M,G] = eegarsim(K,Source,xp,A,snr)
% EEGARSIM to simulate the EEG signal with autoregressive model
%
% inputs: 
%       K: leadfield
%       Source: structure with vertice, amplitud and spread of each source
%       xp: information related with the spectral analysis (fs, nfft, ...)
%       A: coeficients of auto-regressive model model
%       snr: signal to noise variance ratio
% outputs:
%       V0: Pure EEG signals
%       V1: Total EEG signals
%         J : current density 
%       M:  patch position
%       G:  patch source time signal

% Eduardo, Andy Hu, April 2017

ns=size(Source,2); % number of activated patchs
nv=xp.nv; % number of voxels
nt=xp.nt; % number of time samples
M = zeros(nv,ns);

for ii = 1:ns
    M(Source(ii).Vert,ii) = Source(ii).Amp;
end;

% Creating direct wiring process via AR
p = size(A,3); % the order of AR model
G = ones(ns,nt); 

for tt = p+1:nt
    e = randn(ns,1); 
    for ll = 1:p, G(:,tt) = A(:,:,ll)*G(:,tt-ll)+e; end
end

G(:,1:xp.nfft*2) = []; %plotG;  % burn in time unstable components...
% Creating Data, here should be noted that the operation V = (K*M)*G allow
% to go for higher dimentionality of time serie instead of J = M*G and V = K*J...

J = M*G*20;
V0 = K*J;

% Adding noise
[V1, V0] = adwn(V0, snr);

end