function pairs = compute_closest_points(dmat, known_pairs1, known_pairs2)
% compute_closest_points.m
%   Assigns correspondence between two point patterns based on distances
%   between them. In case, we have known pairs, we consider their priority
%   in assignment is 
%
% Input:
%    dmat: Distance matrix (d_ij is distance between x_i and y_j)
%    known_pairs1 & 2: node id's of known pairs in point set 1 & 2
%    
% Output:
%    node id's of corresponding pairs
%
% Note:
%    known_pairs 1 & 2 are column vectors
%

[m n] = size(dmat);
rowbook = ones(m, 1);
colbook = ones(1, n);
num_pairs = 0;
total = min(m, n); % The max number of one-to-one correspondence
pairs = zeros(total, 2);

%
% Give top priority to the known pairs
%
num_known_pairs = length(known_pairs1);
if num_known_pairs
    pairs(1:num_known_pairs, :) = [known_pairs1 known_pairs2];
    rowbook(known_pairs1) = 0;
    colbook(known_pairs2) = 0;
    num_pairs = num_known_pairs;
end

%
% Traverse sorted distances in search of unassigned pairs
%
[dvec_sorted rankvec] = sort(dmat(:), 1);
trav = 1;
while num_pairs < total
%     [rowpoint colpoint] = get_index(rankvec(trav), m, n);
    [rowpoint colpoint] = ind2sub([m n], rankvec(trav));
    if rowbook(rowpoint) && colbook(colpoint) % A pair found
        rowbook(rowpoint) = 0;
        colbook(colpoint) = 0;
        num_pairs = num_pairs + 1;
        pairs(num_pairs, :) = [rowpoint, colpoint];
    end
    trav = trav + 1;
end
