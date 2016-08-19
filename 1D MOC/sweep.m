function solution = sweep( igroup, solution, mesh, quad )
%SWEEP Performs 1D MOC sweep for a single ray with multiple polars
%   igroup   - The energy group being swept
%   solution - The solution object to store the solution in
%   mesh     - The mesh object to solve
%   quad     - Quadrature to use for the MOC sweep

for i=1:mesh.nfsrcells
    k = mesh.nfsrcells-i+1;
    for j=1:quad.npol
        % Forward Sweep
        dx = mesh.fsredges(i+1) - mesh.fsredges(i);
        exparg = exp(-mesh.xstr(i)*dx);
        solution.angflux(i+1,j,1,igroup) = solution.angflux(i,j,1,igroup)*exparg + ...
            mesh.source(i)/mesh.xstr(i)*(1 - exparg);
        solution.scalflux(i,igroup) = solution.scalflux(i,igroup) + ...
            0.5*sum(solution.angflux(i:i+1,j,1,igroup));
        
        % Backward Sweep
        dx = mesh.fsredges(k+1) - mesh.fsredges(k);
        exparg = exp(-mesh.xstr(k)*dx);
        solution.angflux(k,j,2,igroup) = solution.angflux(k+1,j,2,igroup)*exparg + ...
            mesh.source(k)/mesh.xstr(k)*(1 - exparg);
        solution.scalflux(k,igroup) = solution.scalflux(k,igroup) + ...
            0.5*sum(solution.angflux(k:k+1,j,1,igroup));
    end
end


end

