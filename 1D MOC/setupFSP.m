function [ mesh ] = setupFSP( source_list, xstr_list, mesh, igroup )
%SETUPFSP Summary of this function goes here
%   source_list - multigroup source array
%   xstr_list   - multigroup transport cross-section array
%   matmesh     - materials mesh (FSR)
%   igroup      - group index

nfinecells = size(mesh.materials,1);
mesh.source(1:nfinecells,1) = 0.0;
mesh.xs = xsClass(1);
for i=1:nfinecells
    mesh.source(i) = source_list(mesh.materials(i),igroup);
    mesh.xs.transport(i) = xstr_list(mesh.materials(i),igroup);
end

end

