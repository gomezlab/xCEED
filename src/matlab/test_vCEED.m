%% test_vCEED.m
	
clear all;


% The input matrices need to be the same size
% Rows/columns need to be in the same order
%[sp1 dmat1] = read_phylip_dst('cog2255.probcons.dst');
%[sp2 dmat2] = read_phylip_dst('cog0020.probcons.dst');
[sp1 dmat1] = read_phylip_dst('../../data/DistanceMatrix1');
[sp2 dmat2] = read_phylip_dst('../../data/DistanceMatrix2');

rho = corr(squareform(dmat1)', squareform(dmat2)', 'type', 'Pearson')

% Superimpose
% Aligned pairs with error above the threshold will be regarded as an outlier.
threshold = 0.05;
[wRMSD12 Yhat Ystd param] = vCEED(dmat1, dmat2, 'huber', threshold);

% Output: [wRMSD12 Yhat Ystd param]
%    (1) Weighted rmsd between protein family 1 and 2
%    (2) Embedding of protein family 1
%    (3) Embedding of protein family 2 superimposed upon 1
%    (4) Transformation parameters

W12   = param.W;
err12 = diag((Ystd-Yhat)*(Ystd-Yhat)');

[errval errrank] = sort(err12)

% Output from the following figure can be rotated w/ the mouse to get a
% better view.
figure
draw_3d_snapshot(Yhat, Ystd, 'k', 'k', [0.5, 0.5, 0.5], [], [], [], [])

%
% !!! Labeling assumes that sp1 = sp2 (same species and order for both).
%
% Show errors associated with superimposition
figure
subplot(1, 2, 1);
labels = sp1(errrank);
yticks = 1:length(sp1);
barh(W12(errrank), 'FaceColor', [0.5 0.5 0.5])
set(gca, 'YTick', yticks);
set(gca,'YTickLabel', labels);
xlabel('Weighted Error for Pair');
% The following is similar to Figure 4 in the 2009 paper.
subplot(1, 2, 2);
barh(errval, 'FaceColor', [0.5 0.5 0.5])
set(gca, 'YTick', yticks);
set(gca,'YTickLabel', labels);
xlabel('Error Associated with Pair');


% List of species names and their error
% Comment out if you don't want this wriiten to a file.
fid = fopen('vceed_errors.txt','w');
fprintf(fid, 'Tree 1\tTree 2\tError\n');
for i = 1:length(sp1)
    fprintf(fid, '%s\t%s\t%f\n', sp1{i}, sp2{i}, err12(i));
end
fclose(fid);


% Show trees as based on the distance matrices
tree1 = seqneighjoin(dmat1, 'firstorder', sp1);
phytreetool(tree1);
tree2 = seqneighjoin(dmat2, 'firstorder', sp2);
phytreetool(tree2);
