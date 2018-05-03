function H = Hsc(n,indices)
%   Reference transformation matrix for single channel reference
%   Detailed explanation goes here
% See Har, Hrt

% Andy Hu, Apr 10, 2017

if nargin==1, indices=n; end

f=zeros(n,1);

f(indices)=1;

H=eye(n) - ones(n,1)*f';

end

