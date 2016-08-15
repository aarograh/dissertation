close all, clear all, clc;

%% Material IDs
nmats = 8;
id_mod = 1;
id_fuel = 2;
id_control = 3;
id_clad = 4;
id_gap = 5;
id_gtube = 6;
id_controlmod = 7;
id_controlgap = 8;

%% User Options
% Core Information
npins = 3;
pinmap = ['fuel   '; 'control'; 'fuel   '];
pintypes = ['fuel   ';'control'];
cells = [id_mod,id_clad,id_gap,id_fuel,id_gap,id_clad,id_mod;
    id_mod,id_clad,id_gap,id_control,id_gap,id_clad,id_mod];

% Pin information
pitch = 1.26;
fuelrad = 0.4096;
controlrad = 0.4;
cladrad = [0.418, 0.475];
gtuberad = cladrad;
modsize = 1; % 1: narrow water; 2: wide water

% Quadrature
npol = 1;

% Mesh
nmesh_mod = 2;
nmesh_fuel = 5;
nmesh_control = 9;
nmesh_gap = 1;
nmesh_clad = 1;
nmesh_gtube = nmesh_clad;

% Cross-sections
igroup = 8;
xstr_list = ones(nmats,8); %(ngroups,nmats)
source_list = ones(nmats,8); %(ngroups,nmats)

%% Setup Problem
% Geometry
mesh_list = [nmesh_mod; nmesh_fuel; nmesh_control; nmesh_clad; nmesh_gap; nmesh_gtube; nmesh_mod; nmesh_gap];
tmp = pitch/2.0;
if (modsize == 2)
    tmp = sqrt(2)*tmp;
end
widths(id_mod,1) = tmp - cladrad(2);
widths(id_fuel,1) = fuelrad*2.0;
widths(id_control,1) = controlrad*2.0;
widths(id_clad,1) = cladrad(2) - cladrad(1);
widths(id_gap,1) = cladrad(1) - fuelrad;
widths(id_gtube,1) = gtuberad(2) - gtuberad(1);
widths(id_controlmod,1) = tmp - gtuberad(2);
widths(id_controlgap,1) = gtuberad(1) - controlrad;
widths = widths./mesh_list;

% Ray Parameters
if (npol == 1)
    mu = pi/4;
    mucos = cos(mu);
    weights = 2*pi;
elseif (npol == 2)
elseif (npol == 4)
elseif (npol == 8)
elseif (npol == 16)
end

% Mesh
ncells = 0;
for i=1:npins
    for j=1:size(pintypes,1)
        if strcmp(pinmap(i,:),pintypes(j,:))
            for k=size(cells,2):-1:1
                if (cells(j,k) ~= 0)
                    break
                end
            end
            k = k-1;
            break
        end
    end
    regions = cells(j,1:k)';
    for j=1:size(regions)
        for k=1:mesh_list(regions(j))
            ncells = ncells+1;
            matmesh(ncells) = regions(j);
            mesh(ncells) = widths(matmesh(ncells));
        end
    end
end

% Setup source and cross-sections
sourcemesh(1:ncells) = 0.0;
xstrmesh(1:ncells) = 0.0;
for i=1:ncells
    sourcemesh(i) = source_list(matmesh(i),igroup);
    xstrmesh(i) = xstr_list(matmesh(i),igroup);
end

% Set up solution variables
% Angular flux at cell boundaries (ncells+1, npol, 2); third index is forward and backward directions
angflux_edge(1:ncells+1,1:npol,1:2) = 0.0;
% Average angular flux at cell center (ncells, npol, 2); third index is forward and backward directions
angflux_avg(1:ncells,1:npol,1:2) = 0.0;
% Scalar flux at cell center
scalflux(1:ncells) = 0.0;

%% Solve Problem

for i=1:ncells
    k = ncells-i+1;
    for j=1:npol
        % Forward Sweep
        tmp = exp(-xstrmesh(i)*mesh(i)/mucos(j));
        angflux_edge(i+1,j,1) = angflux_edge(i,j,1)*tmp + sourcemesh(i)*(1.0-tmp);
        angflux_avg(i,j,1) = 0.5*(angflux_edge(i,j,1) + angflux_edge(i+1,j,1));
        scalflux(i) = scalflux(i) + angflux_avg(i,j,1)*weights(j);
        
        % Backward Sweep
        tmp = exp(-xstrmesh(k)*mesh(k)/mucos(j));
        angflux_edge(k,j,2) = angflux_edge(k+1,j,2)*tmp + sourcemesh(k)*(1.0-tmp);
        angflux_avg(k,j,2) = 0.5*(angflux_edge(k,j,2) + angflux_edge(k+1,j,2));
        scalflux(k) = scalflux(k) + angflux_avg(k,j,2)*weights(j);
    end
end

%% Generate Plots
% Setup
position_edge(1:ncells+1) = 0.0;
for i=2:ncells+1
    position_edge(i)=position_edge(i-1)+mesh(i-1);
end
position_center = mesh(1)/2.0;
for i=2:ncells
    position_center(i) = position_center(i-1) + 0.5*(mesh(i) + mesh(i-1));
end
position_all(1:2:2*ncells+1,:,:) = position_edge;
position_all(2:2:2*ncells,:,:) = position_center;

% Plot Angular Flux v. position
angflux_plot(1:2:2*ncells+1,:,:) = angflux_edge;
angflux_plot(2:2:2*ncells,:,:) = angflux_avg;

figure(1);
ipol=1;
plot(position_all,angflux_plot(:,ipol,1),'b',position_all,angflux_plot(:,ipol,2),'r')
xlabel('Position (cm)')
ylabel('Angular Flux')

% Scalar Flux v. position
figure(2);
plot(position_center,scalflux)
xlabel('Position (cm)')
ylabel('Scalar Flux (cm^-2)')