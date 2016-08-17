function [ mesh ] = setupFSP( solution, source_list, xsLib, mesh, igroup )
%SETUPFSP Sets up source and XS mesh for fixed source MOC problem
%   solution    - The solution data to use to determine the source
%   source_list - multigroup source array
%   xsLib       - Cross-section library object
%   mesh        - The mesh for this problem
%   igroup      - Group index

mesh.source(1:mesh.nfsrcells,1) = 0.0;
for i=1:mesh.nfsrcells
    mesh.source(i) = source_list(mesh.materials(i),47);
    mesh.xs.transport(i) = xsLib.xsSets(mesh.materials(i)).transport(igroup);
end

end