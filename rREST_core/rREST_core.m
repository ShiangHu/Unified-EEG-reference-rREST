function [data1, H1, L, s, lmd] = rREST_core(data,K)
% rREST_core
% Input
%         data: EEG data (channels-time)
%             K: lead field (channels-dipole)
% Output
%         data1: redefined eeg data
%              H1: redefined reference transfer matrix
%                L: regularization matrix
%                s: singular vaules
%            lmd: optimal regularization paramter
% See also: 
%                   https://doi.org/10.3389/fnins.2018.00297
%

% Shiang Hu, 07/30/2019


[Nc,Nt] = size(data);
[data1,H1,L] = rREST_vrbtrans(data,K);

[~,S] = svd(H1'*H1);
s = diag(S);

lmds = rREST_genlmd(1000);
rss = zeros(size(lmds));
df   = zeros(size(lmds));
data1_scaled=data1./norm(data1,'fro');

for i=1:length(lmds)
    lmd = lmds(i);
    data2 = H1*diag(1./(s+lmd))*H1'*data1_scaled;
    rss(i) = norm(data2 - data1_scaled,'fro')^2;
    df(i)  = sum(s./(s+lmd));
end

gcv = rss./(Nc*Nt-df).^2;
lmd = lmds(gcv == min(gcv));
if lmd==min(lmds) || lmd==max(lmds)
    warning('LMD range needs to be shifted!')
end

% figure, 
% subplot(121), semilogx(lmds,gcv), xlabel('LMDs'), ylabel('GCV');
% subplot(122), plot(df,gcv), xlabel('DF'), ylabel('GCV');
end