function [ sourcemesh, xstrmesh ] = setupFSP( source_list, xstr_list, matmesh, igroup )
%SETUPFSP Summary of this function goes here
%   source_list - multigroup source array
%   xstr_list   - multigroup transport cross-section array
%   matmesh     - materials mesh (FSR)
%   igroup      - group index

nfinecells = size(matmesh,1);
sourcemesh(1:nfinecells,1) = 0.0;
xstrmesh(1:nfinecells,1) = 0.0;
for i=1:nfinecells
    sourcemesh(i) = source_list(matmesh(i),igroup);
    xstrmesh(i) = xstr_list(matmesh(i),igroup);
end

end

