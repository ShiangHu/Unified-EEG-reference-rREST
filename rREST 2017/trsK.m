function [K1, K2, K3, K4] = trsK( kd, malos, sbjid, J)
% TRSK to transform the lead field by normalizing, trace averaging
% Usage:
%              Simlation  [K1, K2, K3, K4] = trsK( kd, [], sbjid, J)
%              Validation [K1, K2, K3] = trsK( kd, malos, sbjid)
%              Individual LF [K1, K2] = trsK( kd, [], sbjid, J)
%              mean LF       [~, ~, K3, K4] = trsK( kd, malos)
% Inputs:
%        kd ---- the dir file of the raw BEM lead fields
%        malos ---- indices to remove the bad channels
%        subjid --- the indice for the current subject
%        J ---- brain source time course
% Outpus:
%        K ---- a suqare matrix. Here K = K*K';
% See also vrbtrans

% Andy Hu, Nov. 1, 2017
if nargin ==2,
    sbjid = []; J = [];
elseif nargin ==3,
    J=[];
end

K1=[]; K2=[]; K3=[]; K4=[];
idx = importdata('indices-red-6000.mat');

if ~isempty(sbjid)
    % K1 ---- individual lead field
    K = importdata(kd(sbjid).name);
    K=K(:,idx); K(malos,:)=[];
    K1 = (K*K')./trace(K*K');
    
    if ~isempty(J)
        % K2 ---- individual sparse lead field
        K = importdata(kd(sbjid).name);
        K=K(:,idx); K(malos,:)=[];
        Nc = size(K,1);
        K = ones(Nc,1)*(sum(abs(J'))~=0).*K;
        K2 = (K*K')./trace(K*K');
        return;
    end
end

if nargin <= 3
    % K3---- trace averaged lead field
    Nc = 58 - length(malos); ak=zeros(Nc,Nc);
    for i=1:length(kd)
        K = importdata(kd(i).name);
        K=K(:,idx);K(malos,:)=[];
        tmp=K*K' / trace(K*K'); ak=ak+tmp;
    end
    K3 = ak/length(kd);
    
    % K4 ---- Yao spherical lead field
    K = importdata('yaolf.mat'); K(malos,:)=[];
    K4 = (K*K')./trace(K*K');
end

end