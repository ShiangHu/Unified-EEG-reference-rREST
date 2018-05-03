function [bic, nu, LogP]=refbys(cp0,K,mode)
% REFBYS use bayesian approach to estimate the infinity reference and average reference, as well as output model comparison
% Methods: Bayesian Interpolation (Mackay 1992);
%                  Bayesian model averaging (Pedro, 2009);
%                  Wikipedia and Patrick Breheny BST 764: Applied Statistical Modeling
%                  BIC=n*log(RSS/n)+df*log(n);df equals to the trace of hat matrix;
% Input:
%          Nc ---- number of channels
%          cp0 ---- [Nc by 1] EEG potentials with signal channel reference
%          K  ---- lead field matrix for EEG inverse solution
% Output:
%          bic ---- bayesian information criteria; smaller bic is preferred
%          nu ---- EEG potentials with infinity or average reference
%          LogP ---- log of posterior estimate (Evidence)
%   alpha:   Precision of the spatial prior
%   beta:    Precission of the observation noise
% See also EVIDENCER

% Andy Hu, 07/10/2017
if nargin~=3
    error('Input Vr, K, mode')
end

[cp1,D,U,H,L1,L2] = vrbtrans(cp0,K);

if strcmp('ar',mode), L=L1; else L=L2; end

[Nd,Ns]=size(H);

%% Re-estimation process
diff=inf;
c=1;
LogP = [];
bic=[];

% Intialization
alpha = gamrnd(1000,0.001);
beta = gamrnd(1000,0.001);
A=beta*(H')*H+alpha*(L')*L;
nu = beta*pinv(A)*H'*cp1;
Enu = norm(nu).^2./2;
Ed = norm(cp1-nu).^2./2;

while diff>=eps && c<=1000,
    
    c=c+1;
    A=beta*(H')*H+alpha*(L')*L;
    gamma = Ns - alpha*trace(pinv(A));
    alpha=gamma/(2*Enu);
    beta=(Nd-gamma)/(2*Ed);
    nu = beta*pinv(A)*H'*cp1;
    Enu = norm(nu).^2./2;
    Ed = norm(cp1-nu).^2./2;
    
    if  gamma<0||alpha>10||beta>10, continue; end
    if isinf(det(A)), continue; end;
    Ai=pinv(A); t1=isnan(Ai); t2=isinf(Ai);
    if  ismember(1,t1), continue; end;
    if  ismember(1,t2), continue; end;
    
    M = -1/2*(beta*(cp1')*cp1-nu'*A*nu);
    constL = (-Nd/2)*log(2*pi*1/beta*1/alpha);
    LogP(c) = constL - 1/2*log(det(pinv(L'*L))) - log(det(A))  + M;
    
    diff = abs((LogP(c)-LogP(c-1))./LogP(c));
end;

LogP=LogP(end);

Hh=beta*pinv(A)*H'*D'*U';
df=trace(Hh);
H0=Hsc(Nd,Nd);
RSS=norm(cp0-H0*nu)^2;
bic = Nd*log(RSS/Nd)+df*log(Nd);
end