function W = get_gauss_weight_mat(X, Y, scale)
%
% get_gauss_weight_mat.m
%    Converts distances into probabilities with Gaussian distribution, and
%    normalizes each row of the distance matrix the probabilities
%
% Input:
%    X & Y: Coordinates of the point pattern X & Y
%    scale: variance for gaussian distribution
%
% Output:
%    W: Normalized probability matrix (w_ij is the probability of
%    correspondence between x_i and y_j)
%
% Note:
%    get_sqdist.m is internally called to compute squared distances between
%    points
%

Dsq = get_sqdist(X, Y);
ny = size(Dsq, 2);
W = exp(-Dsq/sqrt(scale));
rowsum = sum(W, 2);
W = W ./ repmat(rowsum, 1, ny);
