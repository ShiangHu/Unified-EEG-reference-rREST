function [v_rcon,df, gcv,regs,res,bic,aic]=reforg(v_r, K, mode)
% without variable transformation
% REFBYS use bayesian approach to estimate the infinity reference and average reference, as well as output model comparison
% Methods: Bayesian Interpolation (Mackay 1992);Bayesian model averaging (Pedro, 2009);
%                  Wikipedia and Patrick Breheny BST 764: Applied Statistical Modeling
%                  BIC=n*log(RSS/n)+df*log(n);df equals to the trace of hat matrix;
% Input:
%          v_r ---- [Nc-Nt] EEG potentials with signal channel reference
%          K  ---- lead field matrix for EEG inverse solution
%          mode ---- 'ar' or 'rest'
% Output:
%          v_rcon ---- reconstructed EEG potentials at infinity
%          df ---- degree of freedom
%          gcv ---- generalized cross validation
%          regs ---- regularizers
%          res ---- goodness of fit under the transformed ridge regression
%          bic ---- bayesian information criteria; smaller bic is preferred
% See also EVIDENCER

% Andy Hu, 07/10/2017
if ~ischar(mode) || nargin<3
    error('Input Vr, K, mode')
end

% variable redefination
Nc=size(v_r,1); H1=Hsc(Nc,Nc);
[cp1,~,L1,L2] = vrbtrans(v_r,K);
H1(end,:)=[]; v_r(end,:)=[];

if strcmp('ar',mode), L=L1; else L=L2; end
H=H1*pinv(L);

% reference model comparison
[Ui,si,Vi] = svdrapid(H);s2=si.^2;
[Nc, Nt] =size(cp1);

% np=1000; regs=logspace(-6,3,np);
regs=0.1;np=1;
df=zeros(np,1); gcv=zeros(np,1);
res=zeros(np,1); bic=zeros(np,1);
aic=zeros(np,1);
v_rcon = zeros(Nc,Nt);
% v_rcon = []; 

for i=1:np
    lmd=regs(i);
    
    v_rcon (:,i) = pinv(L)*pinv(H'*pinv(H1*H1')*H+lmd*pinv(L'*L))*H'*v_r;
    
    df(i)=trace(H*pinv(H'*pinv(H1*H1')*H+lmd*pinv(L'*L))*H');
    
    res(i)=norm(v_r-H1*v_rcon,'fro')^2;
    
    Nct = Nc * Nt;
        
    gcv(i)=res(i)/(Nct-Nt*df(i))^2;
    
    aic(i)=Nct*log(res(i)/Nct)+Nt*2*df(i);
    
    bic(i)=Nct*log(res(i)/Nct)+Nt*df(i)*log(Nct);
end

% ind=find( (1<df) & (df<=0.3*Nc));
% regs=regs(ind);
% df=df(ind); gcv=gcv(ind);
% res=res(ind); bic=bic(ind);
v_rcon = v_rcon';
end