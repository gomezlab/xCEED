function [c, ceq] = gceedConstraint(x, d, A, B, dboundsq)
%
% Specifies that the constraint R'*R = R*R' = I, s > 0
%
dsq = d^2;
I = eye(d);
R = reshape(x(1:dsq), d, d);
c = [];
ceq = reshape(R' * R - I, dsq, 1);

%
% Specifies that the known pairs should be close to each other
%
m = size(A, 1);
if m
    t = x(dsq+1:dsq+d)';
    vecs = A * R' + repmat(t, m, 1) - B;
    dist = sum(vecs .* vecs, 2) - dboundsq * ones(m, 1);
    c = dist;
end
end
