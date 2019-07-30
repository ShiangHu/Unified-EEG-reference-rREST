function [index,findex]=surfpatch(Source,vertices,faces,d0)
index = Source;
[row,col] = find(faces==Source);
findex_old = faces(row,:);
findex_old = findex_old(:);
findex_old = setdiff(findex_old,Source);
vect = vertices(findex_old,:) - repmat(vertices(Source,:),[size(findex_old),1]);
dold = sum(abs(vect).^2,2).^(1/2);
findex_old(dold>d0) = [];
dold(dold>d0) = [];
index = [index;findex_old];
findex = [];

for i=1:size(findex_old)
    [row,col] = find(faces==findex_old(i));
    index1 = faces(row,:);
    index1 = index1(:);
    index1 = setdiff(index1,index);
    findex = [findex;index1];
end
d = zeros(size(findex,1),1);
findex_new = [];

for i=1:10
    ii=1;
    while ii<=size(findex,1)
        if sum(findex)==0
            findex_new=[];
            break
        end
        %Search the neighbors in the prvious frontier and evaluate the distance
        [row,col] = find(faces==findex(ii));
        index1 = faces(row,:);
        index1 = index1(:);
        [evalindex,ia,ib] = intersect(findex_old,index1);
        vect = vertices(evalindex,:) - repmat(vertices(findex(ii),:),size(evalindex,1),1);
        d(ii) = min(dold(ia) + sum(abs(vect).^2,2).^(1/2));
        if d(ii)>d0
            findex(ii) = [];
            d(ii) = [];
            continue
        end
        index1 = setdiff(index1,[index;findex]);
        findex_new = [findex_new;index1];
        ii=ii+1;
    end
    dold = d;
    findex_old = findex;
    index = [index;findex];
    findex = findex_new;
    findex_new = [];
    d = zeros(size(findex,1),1);
    if sum(findex)==0
        break
    end
end
end