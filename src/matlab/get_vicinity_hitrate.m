function score = get_vicinity_hitrate(rmat, bound)
%
% get_vicinity_hitrate.m
%    Computes the percentage of points whose known interaction partner is
%        included in its top three (or the given threshold ranking) closest
%        points
%
% Input:
%    rmat: Rankings of nodes in terms of location (or vicinity)
%    bound: Threshold ranking
%
% Output:
%    score: The vicinity hitrate
%
n = size(rmat, 1);
cnt = 0;
for i=1:n
    if find(rmat(i,:) == i) <= bound
        cnt = cnt + 1;
    end
end
score = 100*cnt/n;
