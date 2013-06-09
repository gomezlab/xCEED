function [dmat, dmat_sorted, rmat] = get_nearnei(coord1, coord2)

%
% get_nearnei.m
%    Computes the distances between two point patterns and sorts them in
%    the increasing order
%
% Input:
%    coord1 & 2: Coordinates of point pattern 1 & 2
%
% Output:
%    dmat: Distances between point pattern 1 & 2
%    dmat_sorted: 'dmat' sorted row-wise
%    rmat: Rankings in closeness (Each row vector lists the node id's of
%        the most close point to the farthest
%

[m p] = size(coord1);
n = size(coord2, 1);
Dsqsum = zeros(m, n);
for d = 1:p
    Dsqsum = Dsqsum + (repmat(coord1(:,d), 1, n) - repmat(coord2(:,d)', m, 1)) .^ 2;
end
dmat = sqrt(Dsqsum);
[dmat_sorted, rmat] = sort(Dsqsum, 2);
