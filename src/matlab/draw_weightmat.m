function [] = draw_weightmat(W)

Z = [W zeros(size(W, 1), 1)];
Z = [Z; zeros(1, size(Z, 2))];
% caxis([0, 1])
pcolor(Z)
colormap(1-gray(20))
axis ij
% axis equal
axis normal
set(gca, 'XTickLabel', {})
set(gca, 'YTickLabel', {})
colorbar

end
