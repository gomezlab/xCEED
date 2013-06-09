function X = embed_in_space(dmat, method, sizenei)
%
% embed_in_space.m
%    Projects points that best satisfy the given distance relationship
%
% Input:
%    dmat: Distances between nodes
%    method: Method for embedding (metric MDS, non-metric MDS, or SPE)
%    sizenei: The size of neighborhood if SPE (stochastic proximity
%        embedding method is chosen)
%
% Output:
%    X: The coordinates of points in the embedded structure
%

[W, e]  = cmdscale(dmat);
switch(method)
    case 'metric'
        X = W;
    case 'non-metric'
        opts = statset('MaxIter', 1000); %, 'TolFun', 0.001);
        [X, stress, disparities] = mdscale(dmat, [], 'Start', W, 'Options', opts);
    case 'spe'
        if nargin < 3 || isempty(sizenei)
            sizenei = 0.2;
        end
        [X, stress] = spembed(dmat, W, sizenei);
    otherwise
        error('Only ''metric'', ''non-metric'' or ''spe'' is allowed.');
end
