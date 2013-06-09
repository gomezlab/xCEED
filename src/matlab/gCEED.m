function [model_tr scene pairs_prerec pairs2] = gCEED(dmat1obs, dmat2obs, numrep, neibound, nodes1, nodes2, known_pairs1, known_pairs2, dbound, method)
%
% gCEED.m
% Kwangbom Choi
% 06/03/2009
%
% Aligns two point patterns embedded in Euclidean space and recovers
% correspondence
% 
% Input:
%   dmat1obs & dmat2obs: distance matrix 1 & 2
%   numrep: the number of iteration to perform
%   neighbound: the size of neighborhood (used to cluster embedded nodes/points)
%   nodes1 & 2: the initially assigned node id's
%   known_pairs1 & 2: 
%   dbound: distance constraint between known interacting nodes, if any
%   method: embedding method ('stochastic proximity embedment' is
%       supported)
%
% Output:
%   (1) Coordinates of embedded evolutionary structure (protein family 1) in
%       Euclidean space, superimposed onto that of protein family 2
%   (2) Coordinates of embedded evolutionary structure (protein family 2)
%   (3) node id's of predicted interacting pairs for protein family 1
%   (4) node id's of predicted interacting pairs for protein family 2
%

%
% Embed in 3-space
%
mydim = 3;
[Xobs, LambdaX] = cmdscale(dmat1obs);
[Yobs, LambdaY] = cmdscale(dmat2obs);
[nx px] = size(Xobs);
[ny py] = size(Yobs);
maxdim  = max(px, py);
mindim  = min(px, py);
if maxdim < mydim
    if(px < maxdim)
        model = [Xobs zeros(nx, maxdim-px)];
        scene = Yobs;
    elseif(py < maxdim)
        model = Xobs;
        scene = [Yobs zeros(ny, maxdim-py)];
    else
        model = Xobs;
        scene = Yobs;
    end
elseif mindim < mydim
    if(px < mydim)
        model = [Xobs zeros(nx, mydim-px)];
        scene = Yobs(:, 1:mydim);
    elseif(py < mydim)
        model = Xobs(:, 1:mydim);
        scene = [Yobs zeros(ny, mydim-py)];
    end
else
    model = Xobs(:, 1:mydim);
    scene = Yobs(:, 1:mydim);
end
if px == 1 && py == 1 % Make sure that we are playing at least in 2-dim space
    model = [model zeros(nx, 1)];
    scene = [scene zeros(ny, 1)];
    disp('Matching is in one dimensional space...')
end
if strcmp(method, 'spe')
    [model, stress_model] = spembed(dmat1obs, model, neibound, mydim);
    [scene, stress_scene] = spembed(dmat2obs, scene, neibound, mydim);
end
model = standardize_config(model);
scene = standardize_config(scene);

%
% Alignment by gaussian mixture model
%
[rotation, translation, model_tr] = gmmreg_highdim(model, scene, numrep, known_pairs1, known_pairs2, dbound);
[dmat dump dump] = get_nearnei(model_tr, scene);
pairs = compute_closest_points(dmat, known_pairs1', known_pairs2');
[dump neworder] = sort(pairs(:,1));
pairs_prerec = pairs(neworder, 2);
pairs2 = pairs_prerec;

%
% Define neighborhood, and realign sub-clusters by recursion
%
neighbors = define_neighborhood(neibound, dmat1obs); % Need to modify neibound?

for i=1:nx
    % 'cluster' is the locations of nodes in a particular neighborhood
    cluster = find(neighbors == i);
    if length(cluster) > 2
        % 'partners' is the locations of the cluster-nodes' partner
        partners = pairs2(cluster);
        % 'nodes' 1 & 2 are node id's
        disp(sprintf('A cluster of Node %d is being re-aligned...%s::%s', nodes1(cluster(1)), mat2str(nodes1(cluster)), mat2str(nodes2(partners))))
        dmat1nxt = dmat1obs(cluster, cluster);
        dmat2nxt = dmat2obs(partners, partners);
        [dump known_pairs_nxt dump] = intersect(cluster, known_pairs1);
        if ~isempty(find(dmat1nxt > 1e-8)) && ~isempty(find(dmat2nxt > 1e-8))
            [dump1 dump2 dump3 cluster_pair] = gCEED(dmat1nxt, dmat2nxt, ceil(numrep / 2), neibound - 0.1, ...
                nodes1(cluster), nodes2(partners), known_pairs_nxt, known_pairs_nxt, dbound, []);
            pairs2(cluster) = partners(cluster_pair);
        end
    end
end
