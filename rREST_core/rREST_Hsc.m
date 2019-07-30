function H = rREST_Hsc(n,indices)
%   Reference transformation matrix for single channel reference
%   Detailed explanation goes here
%
% See Har, Hrt
%   https://iopscience.iop.org/article/10.1088/1741-2552/aaa13f
%   https://link.springer.com/article/10.1007/s10548-019-00706-y

% Andy Hu, Apr 10, 2017

if nargin==1, indices=n; end

f=zeros(n,1);

f(indices)=1;

H=eye(n) - ones(n,1)*f';

end

