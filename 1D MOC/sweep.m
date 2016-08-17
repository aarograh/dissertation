function solution = sweep( solution, mesh, quad )
%SWEEP Performs 1D MOC sweep for a single ray with multiple polars
%   solution - The solution object to store the solution in
%   mesh     - The mesh object to solve
%   quad     - Quadrature to use for the MOC sweep

nfinecells = size(mesh.fsredges,1)-1;
solution = solution.update( );
solution.angflux_cellavg(1:nfinecells,1:quad.npol,1:2) = 0.0;
solution.scalflux(:,:,2) = solution.scalflux(:,:,1)
solution.scalflux(1:nfinecells) = 0.0;

for i=1:nfinecells
    k = nfinecells-i+1;
    for j=1:quad.npol
        % Forward Sweep
        dx = mesh.fsredges(i+1)-mesh.fsredges(i);
        tmp = exp(-mesh.xs.transport(i)*dx/quad.cosines(j));
        solution.angflux(i+1,j,1) = solution.angflux(i,j,1)*tmp + mesh.source(i)*(1.0-tmp);
        solution.angflux_cellavg(i,j,1) = 0.5*(solution.angflux(i,j,1) + solution.angflux(i+1,j,1));
        solution.scalflux(i) = solution.scalflux(i) + solution.angflux_cellavg(i,j,1)*quad.weights(j);
        
        % Backward Sweep
        dx = mesh.fsredges(k+1)-mesh.fsredges(k);
        tmp = exp(-mesh.xs.transport(k)*dx/quad.cosines(j));
        solution.angflux(k,j,2) = solution.angflux(k+1,j,2)*tmp + mesh.source(k)*(1.0-tmp);
        solution.angflux_cellavg(k,j,2) = 0.5*(solution.angflux(k,j,2) + solution.angflux(k+1,j,2));
        solution.scalflux(k) = solution.scalflux(k) + solution.angflux_cellavg(k,j,2)*quad.weights(j);
    end
end


end

