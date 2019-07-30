function [re,rebm,oli,omi] = refpre(v_rcon,msi,v0,v_r,K1, K2, K3, K4)
% REFPRE to calculate the relative error of the estimated potentials against the ground tureth
% Reference estimates could be AR, rAR, REST, rREST.
% For REST and rREST, the lead field (LF) could be individual LF, Sparse individual LF, trace avg LF, Yao LF
% This only works for the theoretical simulation.
% Inputs:
%              v_rcon ---- [Np*Ne*Nt*5] estimated EEG potentials at infinity by AR(identity prior), REST (lead field prior)
%              msi  ---- [Np*6*5] msi=[lmds, df, res, gcv, aic, bic]; AR, individual LF, Sparse individual LF, trace avg LF, Yao LF
%              v0   ---- [Ne*Nt] simulated EEG pure signal taken as the ground truth
%              v_r  ---- [Ne*Nt] observed EEG potentials with last channel reference
%              K    ---- [Ne*Nv] lead field matrix

% outputs:
%              rebm ---- [1*10] 'best' pre by least model selection indexes
%              omi ---- [1*6] optimal lambda indices by the least gcv,aic,bic
%              re   ---- [1*10] real best pre by least re
%              oli ---- [1*2] optimal lambda indices if the ground true is known
% See also plotcre

% Andy Hu, Nov. 9, 2017

%--------------------RE by rAR, rREST--------------------
% Nr is the number of lambdas (regularization parameter)
[Np, Nc, Nt] = size(v_rcon(:,:,:,1));

v01 = reshape(v0,[1,Nc*Nt]);
v_rcon=reshape(v_rcon,[Np,Nc*Nt,5]);

RE = zeros(Np,5);
for i=1:5
    RE(:,i)=sqrt(sum((v_rcon(:,:,i)-ones(Np,1)*v01).^2,2)) ./ norm(v01,'fro');  % rAR, rREST[Y,I,A,S]
end

% last lmd (extermly samll) is for AR and REST
% re(end,1) is re of AR; re(end,2:5) is re of REST only if tol=[] in Hrt;

%-----------------RE of truncation REST-------------
for i=1:4
    if i==1, K=K1;          % ind lf
    elseif i==2, K=K2;   % Sparse inf lf
    elseif i==3, K=K3;   % trace averaged lf
    elseif i==4, K=K4;   % Yao lf
    end
    [U,S] = svd(K); K = U*(S^0.5)*U';
    v_rt = Hrt(Nc,K,0.05)*v_r;
    RE(end,i+1) = norm(v_rt-v0,'fro') / norm(v0, 'fro');
end

%--------------------Least RE & oli----------------------
[Y,oli] = min(RE(1:end-1,:));
re = [RE(end,1), Y(1), RE(end,2:5), Y(2:5)]; % AR, rAR, REST, rREST

rebm = [];oli = []; omi = [];
return;
%--------------------Least RE & omi----------------------
% return the LMDs indices with the least GCV, AIC, BIC
[~,omi] = min(msi(1:end-1,4:6,:)); omi = squeeze(omi);
[~,omi_ar] = min(msi(:,4:6,1)); % last row expected
%         AR-[gab], rAR-[gab],               REST(I,S,A,Y), rREST-[gab]-I, rREST-[gab]-S, rREST-[gab]-A, rREST-[gab]-Y
rebm = [RE(omi_ar,1)', RE(omi(:,1),1)', RE(end,2:5), RE(omi(:,2),2)', RE(omi(:,3),3)', RE(omi(:,4),4)', RE(omi(:,5),5)'];

omi(:,2:6)=omi; omi(:,1)=omi_ar;
end