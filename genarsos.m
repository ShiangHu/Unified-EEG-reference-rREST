function [Source, A] = genarsos(patchtype,vertices,faces)

% The AR model and source are generated.
%
% input:
%      patchtype: the type of manually defined patch sources
%                        1, Paz-Linares (2017); 2, Baccala and Sameshima (2001)/Pascual-Marqui (2014)
% output:
%       Source: structure with Vertices (Vert), amplitud (Amp) and patch spread distance (d)
%       A: tensor with AR model values
% tips: to localize a source find(abs(vertices(:,1) + repmat(xx.xx,size(vertices,1),1))<0.01)
% Deriel, Andy Hu, 2017

switch patchtype
    case 1
        a1 = 5; a2 = 5;   % amplitude
        d1 = 9; d2 = 9;   % patch spread distance
        idx1=fdvert(1689,vertices,faces,d1);
        idx2=fdvert(1520,vertices,faces,d2);
        
    case 2
        a1 = 5; a2 = 5;
        d1 = 9; d2 = 30;
        idx1=fdvert(1689,vertices,faces,d1);
        idx2=fdvert(4421,vertices,faces,d2);
        
    case 3
        a1 = 5; a2 = 5;
        d1 = 30; d2 = 9;
        idx1 =fdvert(4869,vertices,faces,d1);
        idx2 =fdvert(4421,vertices,faces,d2);
        
    case 4
        a1 = 5;   a2 = 5;
        d1 = 30; d2 = 30;
        idx1 =fdvert(4869,vertices,faces,d1);
        idx2 =fdvert(1520,vertices,faces,d2);
        
end

Source(1).Amp = a1;Source(2).Amp = a2;
Source(1).Patch = d1;Source(2).Patch = d2;
Source(1).Vert = idx1;Source(2).Vert = idx2;

A = zeros(2,2,4);
A(:,:,1) = [1 0; 0 1];       A(:,:,2) = [-0.25 0;0 -0.25];
A(:,:,3) = [0 0.2;-0.2 0]; A(:,:,4) = [0 0.1;0 0];

end