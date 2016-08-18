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
        dx = mesh.fsredges(i+1)-mesh.fsredges(i);
        tmp = exp(-mesh.xstr(i)*dx/quad.cosines(j));
        solution.angflux(i+1,j,1,igroup) = solution.angflux(i,j,1,igroup)*tmp + mesh.source(i)*(1.0-tmp);
        solution.angflux_cellavg(i,j,1,igroup) = ...
            0.5*(solution.angflux(i,j,1,igroup) + solution.angflux(i+1,j,1,igroup));
        solution.scalflux(i,igroup) = solution.scalflux(i,igroup) + ...
            solution.angflux_cellavg(i,j,1,igroup)*quad.weights(j);
        
        % Backward Sweep
        dx = mesh.fsredges(k+1)-mesh.fsredges(k);
        tmp = exp(-mesh.xstr(k)*dx/quad.cosines(j));
        solution.angflux(k,j,2,igroup) = solution.angflux(k+1,j,2,igroup)*tmp + mesh.source(k)*(1.0-tmp);
        solution.angflux_cellavg(k,j,2,igroup) = ...
            0.5*(solution.angflux(k,j,2,igroup) + solution.angflux(k+1,j,2,igroup));
        solution.scalflux(k,igroup) = solution.scalflux(k,igroup) + ...
            solution.angflux_cellavg(k,j,2,igroup)*quad.weights(j);
    end
end


end

