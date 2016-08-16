function [ angflux_edge,angflux_avg,scalflux ] = sweep( finemesh, xstrmesh, sourcemesh, BC, quad )
%SWEEP Performs 1D MOC sweep for a single ray with multiple polars
%   finemesh  - Flat Source Region (FSR) cell-edge positions
%   xstrmesh  - FSR transport cross-sections (1-group)
%   sourcemesh - FSR sources (1-group)
%   BC         - Boundary condition for the flux
%   quad       - Quadrature to use for the MOC sweep

nfinecells = size(finemesh,1)-1;
angflux_edge(1:nfinecells+1,1:quad.npol,1:2) = BC;
angflux_avg(1:nfinecells,1:quad.npol,1:2) = 0.0;
scalflux(1:nfinecells) = 0.0;

for i=1:nfinecells
    k = nfinecells-i+1;
    for j=1:quad.npol
        % Forward Sweep
        dx = finemesh(i+1)-finemesh(i);
        tmp = exp(-xstrmesh(i)*dx/quad.cosines(j));
        angflux_edge(i+1,j,1) = angflux_edge(i,j,1)*tmp + sourcemesh(i)*(1.0-tmp);
        angflux_avg(i,j,1) = 0.5*(angflux_edge(i,j,1) + angflux_edge(i+1,j,1));
        scalflux(i) = scalflux(i) + angflux_avg(i,j,1)*quad.weights(j);
        
        % Backward Sweep
        dx = finemesh(k+1)-finemesh(k);
        tmp = exp(-xstrmesh(k)*dx/quad.cosines(j));
        angflux_edge(k,j,2) = angflux_edge(k+1,j,2)*tmp + sourcemesh(k)*(1.0-tmp);
        angflux_avg(k,j,2) = 0.5*(angflux_edge(k,j,2) + angflux_edge(k+1,j,2));
        scalflux(k) = scalflux(k) + angflux_avg(k,j,2)*quad.weights(j);
    end
end


end

