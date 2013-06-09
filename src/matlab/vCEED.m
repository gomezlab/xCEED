function [wRMSD, Yhat, Ystd, param] = vCEED(dmatX, dmatY, method, threshold)
%
% vCEED.m
% 
%   Computes the weighted rmsd based on the robust superimposition of two
%   input distance matrices.
%
% Input:
%    (1) distance matrix of protein family 1
%    (2) distance matrix of protein family 2
%
% Output:
%    (1) Weighted rmsd between protein family 1 and 2
%    (2) Embedding of protein family 1
%    (3) Embedding of protein family 2 superimposed upon 1
%    (4) Transformation parameters
%

X = embed_in_space(dmatX, 'metric', []);
Y = embed_in_space(dmatY, 'metric', []);

if isempty(method)
    method = 'huber';
end

if isempty(threshold)
    threshold = 0.05;
end

[xn xp] = size(X);
[yn yp] = size(Y);
maxdim  = max(xp, yp);
X = [X zeros(xn, maxdim-xp)];
Y = [Y zeros(yn, maxdim-yp)];
Ystd = standardize_config(Y);
% Yhat is transformed X
[wRMSD Yhat param] = robust_procrustes(X, Ystd, method, threshold);
