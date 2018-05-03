function [idx,fidx]=fdvert(csos,vertices,faces,d0)
% find the vertices index of patch source in the 6K space
% Inputs:
%            csos ---- central source
%            d0 ---- patch spread distance
% Outputs:
%            idx ---- the idx of the patch sources
% See also surfpatch

% Andy, 2017


% take the central source as the seed to build the first circle
[row,~] = find(faces==csos); % seed
fidxod = faces(row,:); fidxod = setdiff(fidxod(:),csos); %find old index
vect = vertices(fidxod,:) - repmat(vertices(csos,:),size(fidxod));
dold = sum(vect.^2,2).^(1/2); % distance
fidxod(dold>d0) = []; dold(dold>d0) = [];
idx = [csos;fidxod];

% search toward outside from the first circle
fidx = [];
for i=1:size(fidxod)
    [row,~] = find(faces==fidxod(i)); % seed
    idx1 = faces(row,:); idx1 = setdiff(idx1(:),idx);
    fidx = [fidx;idx1];
end

% exlarge the searching area
d = zeros(size(fidx)); fidxnw = [];  % find index new
for i=1:10
    ii=1;
    while ii<=size(fidx,1),
        if sum(fidx)==0, fidxnw=[]; break; end;
        [row,~] = find(faces==fidx(ii)); % seed
        idx1 = faces(row,:); [evalindex,ia,~] = intersect(fidxod,idx1(:));
        vect = vertices(evalindex,:) - repmat(vertices(fidx(ii),:),size(evalindex));
        d(ii) = min(dold(ia) + sum(vect.^2,2).^(1/2));
        if d(ii)>d0, fidx(ii) = []; d(ii) = []; continue; end;
        idx1 = setdiff(idx1,[idx;fidx]);
        fidxnw = [fidxnw;idx1];
        ii=ii+1;
    end
    
    idx = [idx;fidx];
    fidxod = fidx;
    fidx = fidxnw;  fidxnw = [];
    dold = d; d = zeros(size(fidx));
    if sum(fidx)==0, break; end
end

idx = clnidx(idx);
end

function newidx = clnidx(idx)
% clean index
idx = sort(idx);
for i = length(idx)-1:-1:1
    if idx(i)==idx(i+1), idx(i+1) = []; end;
end;
newidx = idx;
end