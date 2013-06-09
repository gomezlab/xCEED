function [] = draw_3d_snapshot(Z1, Z2, markercolor1, markercolor2, linecolor, labels1, labels2, pairs1, pairs2)
%
% draw_3d_snapshot.m
%    Displays point patterns and the correspondence between points
%
% Input:
%    Z1 & Z2: Coordinates of point pattern 1 & 2
%    markercolor1 & 2: Color for point pattern 1 & 2 
%    linecolor: Color for edge drawn between corresponding pairs
%    labels1 & 2: Text label for each point
%    pairs1 & 2: Node id's of corresponding pairs
%
% Output:
%    None, but a figure window will pop up with the 3D plot
%
% Note:
%    If labels1 = [], then no label will be printed.
%

hold on;
box on;
axis vis3d;
rotate3d on;

[n1 p1] = size(Z1);
[n2 p2] = size(Z1);

if isempty(labels1)
    labels1 = 1:size(Z1, 1);
    labels1 = num2str(labels1');
end
if isempty(labels2)
    labels2 = 1:size(Z2, 1);
    labels2 = num2str(labels2');
end

if p1 == 2 && p2 == 2
    objh1 = scatter(Z1(:,1), Z1(:,2));
    set(objh1, 'Marker', '.', 'SizeData', 96, 'MarkerEdgeColor', markercolor1)%, 'MarkerFaceColor', markercolor1)
    text(Z1(:,1), Z1(:,2), labels1, 'color', markercolor1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    objh2 = scatter(Z2(:,1), Z2(:,2));
    set(objh2, 'Marker', 'o', 'SizeData', 48, 'MarkerEdgeColor', markercolor2)%, 'MarkerFaceColor', markercolor2)
    text(Z2(:,1), Z2(:,2), labels2, 'color', markercolor2, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    if isempty(pairs1)
        linex = [Z1(:, 1), Z2(:, 1)]';
        liney = [Z1(:, 2), Z2(:, 2)]';
    else
        linex = [Z1(pairs1, 1), Z2(pairs2, 1)]';
        liney = [Z1(pairs1, 2), Z2(pairs2, 2)]';
    end
    line(linex, liney, 'color', linecolor, 'LineWidth', 1);
elseif p1 > 2 && p2 > 2
    objh1 = scatter3(Z1(:,1), Z1(:,2), Z1(:,3));
    set(objh1, 'Marker', '.', 'SizeData', 96, 'MarkerEdgeColor', markercolor1)%, 'MarkerFaceColor', markercolor1)
    text(Z1(:,1), Z1(:,2), Z1(:,3), labels1, 'color', markercolor1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    objh2 = scatter3(Z2(:,1), Z2(:,2), Z2(:,3));
    set(objh2, 'Marker', 'o', 'SizeData', 48, 'MarkerEdgeColor', markercolor2)%, 'MarkerFaceColor', markercolor2)
    text(Z2(:,1), Z2(:,2), Z2(:,3), labels2, 'color', markercolor2, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    if isempty(pairs1)
        linex = [Z1(:, 1), Z2(:, 1)]';
        liney = [Z1(:, 2), Z2(:, 2)]';
        linez = [Z1(:, 3), Z2(:, 3)]';
    else
        linex = [Z1(pairs1, 1), Z2(pairs2, 1)]';
        liney = [Z1(pairs1, 2), Z2(pairs2, 2)]';
        linez = [Z1(pairs1, 3), Z2(pairs2, 3)]';
    end
    line(linex, liney, linez, 'color', linecolor, 'LineWidth', 1);
end

