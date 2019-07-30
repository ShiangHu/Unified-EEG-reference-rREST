% Esta funcion realiza una descomposicion rapida de una
% matriz en valores singulares.
%
% Sintaxis:
%           [U,d,V]=svdrapid(A);
% 
% Nota: Si alguno de los autovalores d son cero, lo que hay en las columnas
%       de V (si m<n) o de U (si m>n) es cascara. 
%       d sale ordenado de mayor a menor.
function [U,d,V]=svdrapid(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Codigo agregado por Nelson para contemplar el caso m > n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(A);
econ = 0;
if m>n, A = A';econ=1;end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[U,D] = eig(A*A');
d=diag(D);
ind=find(d <= eps*max(d));
d(ind)=zeros(size(ind));
d=sqrt(abs(d));
[d,i]=sort(d);
d=flipud(d(:));
i=flipud(i(:));
U=U(:,i);
ind=find(d > eps*max(d));
invd=zeros(size(d));
invd(ind)=1./d(ind);
V1=A'*U; % V es VD
%V=V*diag(invd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Codigo agregado por Nelson para contemplar el caso m > n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~econ,
    V=V1*diag(invd);
else
    V = U;
    U = V1*diag(invd);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
