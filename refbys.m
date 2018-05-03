function [v_rcon,msi]=refbys(v_r,K,mode)
% REFBYS a Bayesian approach to simultaneously estimate EEG potentials at infinity and denosing, output model selections indexes as well
% Methods: Bayesian Interpolation (Mackay 1992);Bayesian model averaging (Pedro, 2009);
%          Wikipedia and Patrick Breheny BST 764: Applied Statistical Modeling
%          BIC=n*log(RSS/n)+df*log(n);df equals to the trace of hat matrix;
%          Search light: Thousand LMDs are generated to test. Optimal lmd is picked by ground truth or GCV, AIC, BIC 
% Input:
%          v_r ---- [Nc-Nt] EEG potentials referenced by the last channel
%          K  ---- lead field matrix
%          mode ---- 'ar' or 'rest'
% Output:
%          v_rcon ---- [Np-Nc-Nt] reconstructed EEG potentials at infinity
%          msi      ---- [Np-6] model selection indexes [lmds, df, res, gcv, aic, bic]
%          lmds ---- regularization parameter         
%          df ---- degree of freedom
%          res ---- goodness of fitting under the standard ridge regression    
%          gcv ---- generalized cross validation
%          aic ---- akaike information criteria
%          bic ---- bayesian information criteria; smaller bic is preferred
% See also EVIDENCER

% Andy Hu, 07/10/2017

if ~ischar(mode) || nargin<3
    error('Input Vr, K, mode')
end

% variable redefination
[V,H,L] = vrbtrans(v_r,K,mode);

% remove the zero rows 
H(end,:)=[]; V(end,:)=[];

% reference model comparison
[Ui,si,Vi] = svdrapid(H);s2=si.^2;
[Nc, Nt] = size(V);

np=1001;                              % Last lmd is extreme samll value for AR and REST
lmds=genlmd(np-1,mode); % LMDs has slight difference for AR and REST

df=zeros(np,1);  gcv=zeros(np,1); res=zeros(np,1);
aic=zeros(np,1); bic=zeros(np,1); 
v_rcon = zeros(np,Nc+1,Nt); 
% v_rcon = [];  % for real data

for i=1:np
    lmd=lmds(i);Nct = Nc * Nt;
    
    v_rcon (i,:,:) = pinv(L)*Vi*diag(si./(s2+lmd))*Ui'*V;
    
%     t=s2./(s2+lmd); df(i)=sum(t);
    
%     err=(eye(Nc)-diag(t))*Ui'*V;res(i)=norm(err,'fro')^2;
    
%     gcv(i)=res(i)/(Nct-Nt*df(i))^2;
%     aic(i)=Nct*log(res(i)/Nct)+Nt*2*df(i);
%     bic(i)=Nct*log(res(i)/Nct)+Nt*df(i)*log(Nct);
end

msi=[lmds, df, res, gcv, aic, bic];

end