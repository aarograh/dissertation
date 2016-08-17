function [ angflux_edge,angflux_avg,scalflux ] = sweep( mesh, BC, quad )
%SWEEP Performs 1D MOC sweep for a single ray with multiple polars
%   mesh       - The mesh object to solve
%   BC         - Boundary condition for the flux
%   quad       - Quadrature to use for the MOC sweep

nfinecells = size(mesh.fsredges,1)-1;
angflux_edge(1:nfinecells+1,1:quad.npol,1:2) = BC;
angflux_avg(1:nfinecells,1:quad.npol,1:2) = 0.0;
scalflux(1:nfinecells) = 0.0;

for i=1:nfinecells
    k = nfinecells-i+1;
    for j=1:quad.npol
        % Forward Sweep
        dx = mesh.fsredges(i+1)-mesh.fsredges(i);
        tmp = exp(-mesh.xs.transport(i)*dx/quad.cosines(j));
        angflux_edge(i+1,j,1) = angflux_edge(i,j,1)*tmp + mesh.source(i)*(1.0-tmp);
        angflux_avg(i,j,1) = 0.5*(angflux_edge(i,j,1) + angflux_edge(i+1,j,1));
        scalflux(i) = scalflux(i) + angflux_avg(i,j,1)*quad.weights(j);
        
        % Backward Sweep
        dx = mesh.fsredges(k+1)-mesh.fsredges(k);
        tmp = exp(-mesh.xs.transport(k)*dx/quad.cosines(j));
        angflux_edge(k,j,2) = angflux_edge(k+1,j,2)*tmp + mesh.source(k)*(1.0-tmp);
        angflux_avg(k,j,2) = 0.5*(angflux_edge(k,j,2) + angflux_edge(k+1,j,2));
        scalflux(k) = scalflux(k) + angflux_avg(k,j,2)*quad.weights(j);
    end
end


end

