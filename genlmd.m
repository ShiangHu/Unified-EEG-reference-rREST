function lmds = genlmd(np,mode,visible)
% GENLMD to generate the lambdas, regularization parameters
% The LMDs values should be increased monotonously with the pace
% NSR is defined as sensor noise to pure EEG signal variance ratio
% NjSR is defined as sensor noise to brain source signal variance ratio
% Variance may represent the energy of a signal. Since the volume conduction
% acts as a lowpass filter that lower the energy of the source signal. Thus,
%                              NjSR << NSR
% Input:
%         np ---- number of lmds
%         mode ---- if mode = 'AR', LMDs are from the lowest NSR to 10
%                   if model = 'REST',  LMDs are from the lowest NjSR to 1
%         visible ---- plot lmds if visible = 'fig'
% Output:
%         lmds ---- a vector of lmd
% See also eegarsim, logspace, exp, log, linspace

% Andy Hu, Nov 8, 2017
if nargin == 2,
    visible = 'no plot';
end

% 1e-10 is nearly zero for the traditional AR and REST

if strcmp(mode,'ar')
    lmds = [logspace(-3, -1.5, 500), logspace(-1.5+eps, 1, np-500), 1e-10];
else
    lmds = [logspace(-3.5, -2.5, 500), logspace(-2.5+eps, -1, np-500), 1e-10];
end

lmds = lmds';

if strcmp(visible,'fig')
    figure,
    subplot(1,3,1)
    plot(log10(lmds)); ylabel('LMDs in log10'); xlabel(strcat(num2str(np), ' LMDs'));
    subplot(1,3,2)
    plot(lmds(1:500)); ylabel('LMDs'); xlabel('First 500 LMDs');
    subplot(1,3,3)
    plot(lmds(501:np)); ylabel('LMDs'); xlabel(strcat(num2str(np-500), ' LMDs at last'));
end