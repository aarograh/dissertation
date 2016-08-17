function [ mesh ] = setupFSP( source_list, xsLib, mesh, igroup )
%SETUPFSP Sets up source and XS mesh for fixed source MOC problem
%   source_list - multigroup source array
%   xsLib       - Cross-section library object
%   mesh        - The mesh for this problem
%   igroup      - Group index

nfinecells = size(mesh.materials,1);
mesh.source(1:nfinecells,1) = 0.0;
for i=1:nfinecells
    mesh.source(i) = source_list(mesh.materials(i),47);
    mesh.xs.transport(i) = xsLib.xsSets(mesh.materials(i)).transport(igroup);
end

end

