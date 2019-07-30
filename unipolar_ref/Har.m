function H = Har( n )
% The reference transformation matrxi for average reference
% n: number of channels

H=eye(n)-ones(n)/n;

end

