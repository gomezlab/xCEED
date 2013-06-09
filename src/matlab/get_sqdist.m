function Dsq = get_sqdist(X, Y)

[nx, px] = size(X);
[ny, py] = size(Y);

if px ~= py
    error('The dimensionality of two point sets should match.');
end
dim = px;

Dsq = zeros(nx, ny);
for d = 1:dim
    Dsq = Dsq + (repmat(X(:,d), 1, ny) - repmat(Y(:,d)', nx, 1)) .^ 2;
end