function [snr, nsr]= mksnr(i)
% make  signal to noise ratio
% The ratio may be based on variance.
% SNR = 10*log10(Ps/Pn), NSR= Pn/Ps

% Andy Hu, Nov. 10, 2017

% optioned values

nsr  = [1e-2, 1e-1, 10^(-0.8), 10^(-0.4), 10^(-0.2), 10^(-0.1), 10^(-0.05)];

snr = 10*log10(1./nsr);

snr = snr(i); 
nsr = nsr(i);
end