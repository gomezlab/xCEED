function [stdrmsd1R stdrmsd2R stdrmsd12] = rCEED(dmat1, dmat2, dmatR)
%
% rCEED.m
% 
%   Computes the rmsd from indirect superimposition of two given evolutionary
%   structures (or point patterns)
%
% Input:
%    (1) distance vector of protein family 1
%    (2) distance vector of protein family 2
%    (3) distance vector of reference family (e.g., 16S rRNA orthologs)
%
% Output:
%    (1) Normalized rmsd between protein family 1 and reference
%    (2) Normalized rmsd between protein family 2 and reference
%    (3) Normalized rmsd between protein family 1 and 2
%        (which is the rmsd of indirect superimposition)
%

X1 = embed_in_space(dmat1);
X2 = embed_in_space(dmat2);
XR = embed_in_space(dmatR);

[x1n x1p] = size(X1);
[x2n x2p] = size(X2);
[xRn xRp] = size(XR);

maxdim = max(x1p, xRp);
X1 = [X1 zeros(x1n, maxdim-x1p)];
XR1 = [XR zeros(xRn, maxdim-xRp)];
[stdd1R, TOL1, param1] = procrustes(X1, XR1);

maxdim = max(x2p, xRp);
X2 = [X2 zeros(x2n, maxdim-x2p)];
XR2 = [XR zeros(xRn, maxdim-xRp)];
[stdd2R, TOL2, param2] = procrustes(X2, XR2);

[x1n x1p] = size(TOL1);
[x2n x2p] = size(TOL2);
comdim = max(x1p, x2p);
X1   = [X1 zeros(x1n, comdim-x1p)];
X2   = [X2 zeros(x2n, comdim-x2p)];
TOL1 = [TOL1 zeros(x1n, comdim-x1p)];
TOL2 = [TOL2 zeros(x2n, comdim-x2p)];

[stddcheck, TOL2rot, param] = procrustes(TOL1, TOL2);
X2rot = param.b * X2 * param.T + param.c;
sse12 = trace((X1-X2rot)*(X1-X2rot)');
sse1R = trace((X1-TOL1)*(X1-TOL1)');
sse2R = trace((X2rot-TOL1)*(X2rot-TOL1)');
varTOL1 = trace(TOL1'*TOL1);
stdrmsd12 = sqrt(sse12/varTOL1/x1n);
stdrmsd1R = sqrt(sse1R/varTOL1/x1n);
stdrmsd2R = sqrt(sse2R/varTOL1/x1n);
