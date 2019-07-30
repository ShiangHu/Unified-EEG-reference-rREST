function [V,H,L] = rREST_vrbtrans(v_r,K)
%VRBTRANS to transform the variables in reference estimation into ridge
% regression form so as to eliminate the tricky covariance structure of the
% EEG measurement noise
%      General form         ||v_r-H1*phi||^2+lmd*||L*phi||^2
%      Standard form       ||v_r-H*phi'||^2+lmd*||phi'||^2, H=H1*pinv(L), phi'=L*phi
% Input:
%          v_r ----  real EEG data with Hsc
%          K  ---- EEG lead field
% Output:
%          V ---- redefined potentials
%          H ---- redefined ref trans matrix right divided by L
% See also:
%                   https://doi.org/10.3389/fnins.2018.00297

% Shiang Hu, 07/10/2017

% reference transformation
Nc=size(v_r,1); 
H0=rREST_Hsc(Nc);

% redefine H matrix to have the identity covariance of sensor noise
[D,U]=genD(H0); 
H1=D'*U'*H0; 
V=D'*U'*v_r;

% standard ridge regression form
[U_k,S_k]=svd(K*K');
L = U_k*pinv(sqrt(S_k))*U_k'; % full lead field
H=H1*pinv(L);
end

function [D,U] = genD(H)
% GEND is to eliminate the covariance structure sigma^(-2)*pinv(HH^T) of the scalp
% noise. The solution is to take the SVD, take the square root and plugin
% it into the terms in the misfit function ||v-H*phi||^2
% Input: H is the reference transforming matrix
% Output:
%            D is the inverse of the singular value matrix
%            U is the eigvectors

[U,S]=svd(H*H');
d=diag(S);
idx=logical(d~=0);
d(idx)=1./sqrt(d(idx));
D=diag(d);
end