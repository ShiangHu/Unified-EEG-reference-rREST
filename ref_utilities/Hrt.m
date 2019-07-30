function [H_rest,i] = Hrt(n,K,tol,indices)
%HRT REST transformation matrix, it adopts truncating the eigvalues for regularization
% tol: threshold to cut off the tail of eigvalues of HK, Refer to Zhai & Yao, 2004
% Usage: Hrt(n, K)
%              Hrt(n, K, [], indices)
% H is by default the last channel reference
% if add noise, real data or head model mismatched (special cases), tol for REST should be set!
% tol = 0.05 by default; also for real data.
% See also Har, Hsc

% Andy Hu, Apr. 2017

if nargin<4,
        H = Hsc(n);
else H = Hsc(n,indices);
end

i=1;

if nargin==2||isempty(tol), H_rest=K*pinv(H*K);
    
elseif length(tol)==1,  H_rest=K*pinv(H*K,tol); 
    
elseif length(tol)>1,
    
    H_rest = zeros(size(K,1), size(K,1), length(tol));
    for i=1:length(tol)
    H_rest(:,:,i)=K*pinvrt(H*K,tol(i));
    end
    tol = tol - 0.05; [~,i]=min(abs(tol));
 
end
end

