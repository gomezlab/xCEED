function Z = standardize_config(X)
%
% standardize_config.m
%   Re-centers and rescales the point pattern. 
%
% Input:
%    X: Coordinates of points
%
% Output:
%    Z: A transformed coordinates, centered at the origin and rescaled with
%    the standard deviation
%

[n, m]   = size(X);
muX = mean(X,1);
XC = X - repmat(muX, n, 1);
ssqX = sum(XC.^2,1);
normX = sqrt(sum(ssqX)); % == sqrt(trace(XC*XC')) ~ the "centered" Frobenius norm
Z = XC / normX;
