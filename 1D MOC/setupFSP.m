function [ mesh ] = setupFSP( source_list, xsLib, mesh, igroup )
%SETUPFSP Summary of this function goes here
%   source_list - multigroup source array
%   xsLib       - Cross-section library object
%   matmesh     - materials mesh (FSR)
%   igroup      - group index

nfinecells = size(mesh.materials,1);
mesh.source(1:nfinecells,1) = 0.0;
for i=1:nfinecells
    mesh.source(i) = source_list(mesh.materials(i),47);
    mesh.xs.transport(i) = xsLib.xsSets(mesh.materials(i)).transport(igroup);
end

end

