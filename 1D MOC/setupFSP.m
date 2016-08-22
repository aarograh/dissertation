function [ mesh ] = setupFSP( solution, xsLib, mesh, igroup )
%SETUPFSP Sets up source and XS mesh for fixed source MOC problem
%   solution    - The solution data to use to determine the source
%   xsLib       - Cross-section library object
%   mesh        - The mesh for this problem
%   igroup      - Group index

mesh.source(1:mesh.nfsrcells,1) = 0.0;
for i=1:mesh.nfsrcells
    % Use old scalar flux to do Jacobi style iteration
    matID = mesh.materials(i);
    mesh.source(i,1) = (solution.fisssrc(i,1)*xsLib.xsSets(matID).chi(igroup)/solution.keff(1) + ...
        solution.scalflux(i,:,2)*xsLib.xsSets(matID).scatter(igroup,:)')*0.5;
    mesh.xstr(i,1) = xsLib.xsSets(matID).transport(igroup);
end

end