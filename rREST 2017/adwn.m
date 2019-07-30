function [V1, V0, noise] = adwn(V0, snr)
% ADWN to add white noise to simulated potentials or lead field
% Input
%      V0: simulated EEG potentials or lead field;
%     snr: 10log10(signal to noise variance ratio) in dB unit
% Output
%      V1: EEG signal added with gausian noise
% See also preary, MVNRND

% Andy Hu, 4/9/2017

Nc = size(V0,1);

nsr = 10^(-snr/10);

sig = nsr*mean(diag(cov(V0')));
noise = mvnrnd(zeros(1,size(V0,1)),sig*eye(Nc),size(V0,2))';

V1 = V0 + noise; 

end