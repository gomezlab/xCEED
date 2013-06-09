function [species dmat] = read_phylip_dst(file)
%
% read_phylip_dst.m
%    Reads in phylip distance file
%
% Input: distance file name
%
% Output:
%    species: a cell array of species names  
%    dmat: a matrix of distances (double & symmetric)
%

fid = fopen(file, 'r');
str = fgetl(fid);
numsp = sscanf(str, '%i');
species = cell(numsp, 1);
dmat = [];
for i=1:numsp
    numread = 0;
    dmatrow = [];
    str = fgetl(fid);
    [sp remstr] = strtok(str);
    vec = sscanf(remstr, '%g')';
    numread = numread + length(vec);
    dmatrow = horzcat(dmatrow, vec);
    species{i} = sp;
    while numread < numsp
        str = fgetl(fid);
        vec = sscanf(str, '%g')';
        numread = numread + length(vec);
        dmatrow = horzcat(dmatrow, vec);
    end
    dmat = vertcat(dmat, dmatrow);
end
fclose(fid);
