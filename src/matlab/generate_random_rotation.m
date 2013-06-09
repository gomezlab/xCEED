function rotmat = generate_random_rotation(dim, projection)
%
% projection = 0 if projection is not allowed; 1 otherwise.
%
%
X = rand(dim, dim);
Y = rand(dim, dim);
[L, D, R] = svd(X'*Y);
rotmat = L * R';
if ~projection && det(rotmat) < 0
    projmat = eye(dim);
    projmat(dim, dim) = -1;
    rotmat = L * projmat * R';
end
