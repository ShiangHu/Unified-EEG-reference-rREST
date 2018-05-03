function H = Hlm(n, indices)
%HLM Linked mastoids 
% n: number of channels
% indices: number of left and right matoids

f=zeros(n,1);

f(indices)=0.5;

H=eye(n)-ones(n,1)*f';

end

