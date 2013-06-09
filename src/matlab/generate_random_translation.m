function X = generate_random_translation(dim, numpts)

t = - ones(1, dim);
X = rand(numpts, dim) * 2 + repmat(t, numpts, 1);
% draw_graph(X, [], 'k', [0.5 0.5 0.5], 0);
