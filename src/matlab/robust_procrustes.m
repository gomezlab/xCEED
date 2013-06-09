function [wRMSD Yhat param] = robust_procrustes(X, Y, method, threshold)

% Superimpose X onto Y
% Methods implemented are 'huber' and 'biweight' function.
% Threshold is a critical input. Aligned pairs with erros above the
% threshold will be regarded as an outlier. 
% Returns the following parameters of the rigid transformation:
%    c: dilation
%    R: rotation
%    t: translation
%    W: weights

%
% Todo's: Add Lorenzian kernel
%

[nx px] = size(X);
[ny py] = size(Y);

if nx ~= ny
    error('Input point sets should have the same size.')
end

if px ~= py
    error('Input point sets should have the same dimensionality.')
end

% Initialization
ONE = ones(nx, 1);
W   = ones(nx, 1);
wSSE0 = 1;
improvement = 1;

while abs(improvement) > 1e-6
    
    % Compute Y*, X*, and 1*, and the centering matrix B
    sqrtW = diag(sqrt(W));
    Ystar = sqrtW * Y;
    Xstar = sqrtW * X;
    ONEstar = sqrtW * ONE;
    sqlenONEstar = ONEstar' * ONEstar;
    B = eye(nx) - ONEstar * ONEstar' / sqlenONEstar;

    % Solve weighted procrustes
    Xcenter = Xstar' * B;
    Xvarcov = Xcenter * Xstar;
    Yvarcov = Ystar' * B * Ystar;
    Covmat =  Xcenter * Ystar;
    [U Sigma V] = svd(Covmat);
    R = U * V';
    c = trace(R' * Covmat) / trace(Xvarcov);
    t = (Ystar - c * Xstar * R)' * ONEstar / sqlenONEstar;
%     Ystarhat = c * Xstar * R + ONEstar * t';
    Yhat = c * X * R + ONE * t';

    % Compute residues w.r.t. k(i), T(i), and t(i)
%     Ydiff = inv(sqrtW) * (Ystar - Ystarhat);
%     Ydiff = Ystar - Ystarhat;
    Ydiff = Y - Yhat;
    Residue = sqrt(diag(Ydiff * Ydiff'));

    % Compute the standardized distance between Y and Yhat
    wSSE1 = 1 - c^2 * trace(Xvarcov) / trace(Yvarcov);
    wRMSD = sqrt(wSSE1 / nx);
%     dispstr = sprintf('Weighted SSE: %.4f', wSSE1);
    dispstr = sprintf('Weighted wRMSD: %.4f', wRMSD);
    disp(dispstr);

    % Update weights
    switch(method)
        case 'huber'
            for i=1:nx
                if Residue(i) < threshold
                    W(i) = 1;
                else
                    W(i) = threshold / Residue(i);
                end
            end
        case 'biweight'
            for i=1:nx
                if Residue(i) > threshold
                    W(i) = 0;
                else
                    W(i) = (1 - (Residue(i)/threshold)^2)^2;
                end
            end
        otherwise
            error('Only ''huber'' and ''biweight'' are allowed')
    end

    % Update loop condition
    improvement = wSSE0 - wSSE1;
    wSSE0 = wSSE1;

end
param.c = c;
param.R = R;
param.t = t;
param.W = W;
