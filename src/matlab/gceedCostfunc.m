function [f,g] = gceedCostfunc(x, model, scene, scale)
%=====================================================================
% gceedCostfunc.m
% Kwangbom Choi
% 06/28/2008
% Based intensively on Bing Jian's implementation
%=====================================================================
[n d] = size(model);
m = size(scene, 1);
rotation = reshape(x(1:d^2), d, d);
translation = x(d^2+1:end)';
transformed_model = model * rotation' + repmat(translation, n, 1);
mn = m * n;

%
% Compute Gauss transform
%
Dsqsum_ij = zeros(n, m);
grad = zeros(n, d);
for dim = 1:d
    Dsqsum_ij = Dsqsum_ij + (repmat(transformed_model(:,dim), 1, m) - repmat(scene(:,dim)', n, 1)).^2;
end
fij = exp(-Dsqsum_ij ./ scale^2);
f = sum(sum(fij)) ./ mn;

for i = 1:n
    grad(i,:) = sum((repmat(transformed_model(i,:), m, 1) - scene) .* repmat(-2 .* fij(i, :)', 1, d));
end
grad = grad ./ (mn * scale^2);

f = -f;
grad = -grad;
grad = grad';
if isnan(grad)
    grad(find(isnan(grad))) = 0;
end
gm = grad * model;
rotg = gm(:);
tg = sum(grad, 2);
g = cat(1, rotg, tg)';

% if isnan(gm)
%     disp('There is NaN')
% end
% if isinf(gm)
%     disp('There is Inf')
% end

% disp(sprintf('%6.4f %6.4f %6.4f', cond(rotation), cond(model * rotation'), cond(gm)))

