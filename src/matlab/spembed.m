function [Y, stress] = spembed(proximity, Y, rc, numdims)
%
% SPEMBED Perform the Stochastic Proximity Embedding algorithm
%
%   Y = spembed(proximity, Y, rc)
%
% Perform the Stochastic Proximity Embedding algorithm. proximity is the
% proximity information used in the iterative updation proceess. Y is the
% initial coordinates of data points on which SPE has to be applied. The
% dimensionality of the final map is same as that of Y. The size of the
% neighborhood is defined by the variable rc (default = 0.1).
%
%

    if nargin < 2
        error('At least proximity and initial coordinates should be given.');
    end

    if nargin < 3
        rc = 0.5;
    end
    
    if nargin < 4
        numdims = size(Y, 2);
    end

    lambda = 1;                                     % the learning rate
    numpts = size(Y, 1);                            % the size of the input space
    Y = [Y zeros(numpts, numdims-size(Y, 2))];
    numsteps = round((20000 + numpts^2) / numpts);  % the number of steps
    numcycles = 100;                                % the number of cycles
    epsilon = 1e-6;                                 % a small number to avoid /zero
    declambda = (lambda - 0.1) / numcycles;
    
    for i=1:numcycles
        for j=1:numsteps
            indices1 = randperm(numpts);
            indices2 = randperm(numpts);
            d = sqrt(sum((Y(indices1,:) - Y(indices2,:)) .^ 2, 2));
            r = proximity((indices1 - 1) * size(proximity, 1) + indices2)';
            % Get the point pairs to be updated
            qualified = find(indices1' ~= indices2');
            qualified = intersect(qualified, (union(find(d<r), find(r<rc))));
            if ~isempty(qualified)
                indices1 = indices1(qualified);
                indices2 = indices2(qualified);
                d = d(qualified);
                r = r(qualified);
                % Update locations of points
                Y(indices1,:) = Y(indices1,:) + lambda * (1/2) * repmat(((r - d) ./ (d + epsilon)), 1, numdims) .* (Y(indices1,:) - Y(indices2,:));
                Y(indices2,:) = Y(indices2,:) + lambda * (1/2) * repmat(((r - d) ./ (d + epsilon)), 1, numdims) .* (Y(indices2,:) - Y(indices1,:));
            end
        end
        
        % Update learning parameter
        lambda = lambda -  declambda;

    end
    Z = proximity(tril(true(size(proximity,2)),-1));
    Z = Z(:)';                 % force to a row vector, even if empty
    stress = get_stress(Z, Y, rc)
