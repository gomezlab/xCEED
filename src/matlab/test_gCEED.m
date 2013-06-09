clear all;
close all;

[sp1 dmat1orig] = read_phylip_dst('../../data/lytt_regulators_4_tco.dst');
[sp2 dmat2orig] = read_phylip_dst('../../data/lytt_sensors_4_tco.dst');

%
% See the best alignment first
%
[Xorig, LambdaX] = cmdscale(dmat1orig);
[Yorig, LambdaY] = cmdscale(dmat2orig);
[nx px] = size(Xorig);
[ny py] = size(Yorig);
dim = 3;
X3 = Xorig(:, 1:dim);
Y3 = Yorig(:, 1:dim);
X3 = standardize_config(X3);
Y3 = standardize_config(Y3);
[stdd X3proc param] = procrustes(Y3, X3)
figure
draw_3d_snapshot(X3proc, Y3, 'k', 'k', [0.5, 0.5, 0.5], [], [], [], [])

nx = size(dmat1orig, 1);
ny = size(dmat2orig, 1);
known_pairs1orig = []; % known_pairs1orig should be in increasing order
known_pairs2orig = [];

%
% Shuffle matches
%
nummatch = min(nx, ny);
mperm = randperm(nx);
mpoints = mperm(1:nummatch);
dmat1obs = dmat1orig(mpoints, mpoints);
dmat2obs = dmat2orig;
[dump known_pairs1 dump] = intersect(mpoints, known_pairs1orig);
known_pairs2 = known_pairs2orig;
[known_pairs1' known_pairs2']

nodes1 = 1:nx;
nodes2 = 1:ny;
numrep = 128;
neibound = 0.4;
dbound = 0.05;
tic
[model scene pairs1 pairs2] = gCEED(dmat1obs, dmat2obs, numrep, neibound, ...
    nodes1, nodes2, known_pairs1, known_pairs2, dbound, []);
toc

xpart = nodes1;
ypart  = [];
ypart0 = [];
ypart(mpoints) = pairs2;
ypart0(mpoints) = pairs1;
[xpart' ypart' ypart0']
label1 = int2str((1:nx)');
label2 = int2str((1:ny)');
accuracy_strict = 100*length(find(xpart == ypart))/length(xpart)

model_tr(mpoints, :) = model;
figure
draw_3d_snapshot(model_tr, scene, 'k', 'k', [0.5, 0.5, 0.5], label1, label2, xpart, ypart)

Wscaled = get_gauss_weight_mat(model_tr, scene, 0.01);
figure
draw_weightmat(Wscaled);
