function [rotmat, translation, transformed_model, opt_fval] = gmmreg_highdim(model, scene, numreps, known_pairs1, known_pairs2, dbound)
%=====================================================================
% gmmreg_highdim.m
% Kwangbom Choi
% 06/28/2008
% Based intensively on Bing Jian's implementation  
%   'model' and 'scene' are two standardized point sets
%   'scale' is a free scalar parameter
%   'display': display the intermediate steps or not. 
%=====================================================================
if nargin < 2
    error('Usage: gmmreg_highdim(model, scene, numreps)');
end
[nx, px] = size(model);
[ny, py] = size(scene);
if px ~= py
    error('The model and scene should have the identical dimensionality');
end
scale = power(det(model'*model/nx), 1/(2^px));
% scale = 2; % Assuming both model & scene are standardized

% Set up options for optimization algorithm
max_iter =10000;
options = optimset('display',    'off', ...
                   'LargeScale', 'on', ...
                   'GradObj',    'on', ...
                   'TolFun',     1e-8, ...
                   'TolX',       1e-8, ...
                   'TolCon',     1e-8);
% options = optimset(options, 'outputfcn', @outfun);
options = optimset(options, 'Algorithm', 'active-set');
options = optimset(options, 'MaxFunEvals', max_iter);
opt_fval = 1e6;

A = model(known_pairs1,:);
B = scene(known_pairs2,:);
dboundsq = dbound^2;

% fprintf(1, 'Running GMM superimposition: ');
h = waitbar(0, 'Running GMM superimposition...');
% tic
for i=1:numreps
    %     fprintf(1, '%d/%d', i, numreps);
    waitbar(i/numreps, h, sprintf('Running GMM superimposition (%d/%d)\nCurent min: %6.4f', i, numreps, opt_fval))
    init_rotmat = generate_random_rotation(px, 1); init_rotmat = init_rotmat(:);
    init_translation = generate_random_translation(px, 1)';
    x0 = cat(1, init_rotmat, init_translation);
    [x fval] = fmincon(@(x) gceedCostfunc(x, model, scene, scale), x0, [], [], [], [], [], [], @(x) gceedConstraint(x, px, A, B, dboundsq), options);
    if fval < opt_fval
        opt_fval = fval;
        pxsq = px^2;
        rotmat = reshape(x(1:pxsq), px, px);
        translation = x(pxsq+1:end)';
        transformed_model = model * rotmat' + repmat(translation, nx, 1);
    end
end
% toc
% fprintf(1, '\n');
close(h);
end
