function [ mesh ] = setupFSP( solution, xsLib, mesh, igroup )
%SETUPFSP Sets up source and XS mesh for fixed source MOC problem
%   solution    - The solution data to use to determine the source
%   xsLib       - Cross-section library object
%   mesh        - The mesh for this problem
%   igroup      - Group index

mesh.source(1:mesh.nfsrcells,1) = 0.0;
for i=1:mesh.nfsrcells
    % Use old scalar flux to do Jacobi style iteration
    mesh.source(i) = solution.fisssrc(i)*xsLib.xsSets(mesh.materials(i)).chi(igroup) + ...
        sum(solution.scalflux(i,1:igroup,2).*xsLib.xsSets(mesh.materials(i)).scatter(igroup,1:igroup));
    mesh.xs.transport(i) = xsLib.xsSets(mesh.materials(i)).transport(igroup);
end

end