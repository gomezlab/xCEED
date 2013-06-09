function gidvec = define_neighborhood(neibound, dmat) 

% define_neighborhood.m
%    Divide nodes into clusters based on the size of neighborhood
%
% Input:
%    neibound: The size of neighborhood
%    dmat: Distance between nodes
%
% Output:
%    gidvec: Vector of cluster center node id
%

[dmat_sort dmat_rank] = sort(dmat, 2);
indx = find(dmat_sort >= neibound);
neimat = dmat_rank;
neimat(indx) = 0;

numpts = size(dmat, 1);
gidvec = zeros(1, numpts);
curnei = neimat(1, find(neimat(1, :) > 0));

gidvec(curnei) = 1;
for pid=2:numpts
    curnei = neimat(pid, find(neimat(pid, :) > 0));
    gids = unique(nonzeros(gidvec(curnei))');
    numgids = length(gids);
    if numgids
        curgid = gids(1);
        gidvec(curnei) = curgid;
        if numgids > 1 % Found previously unconnected group(s)
            for j = 2:numgids
                othergid = gids(j);
                gidvec(find(gidvec == othergid)) = curgid;
            end
        end
    else
        gidvec(curnei) = pid;
    end
end
