function [ angflux, scalflux, mesh ] = ...
    MOC_1D( pinmap, pitch, diag, pinmats, radii, pinmesh, npol, igroup )
%MOC_1D Solves 1D MOC given geometry and quadrature inputs
%   pinmap  - Map of pins in problem (vector, integer)
%   pitch   - Pitch for each pin (scalar, double)
%   diag    - Flag to indicate of pins are trace horizontally or vertically (0: horizontal, 1: diagonal)
%   pinmats - List of materials in each pin (2D array, integer, one row per pin)
%   radii   - List of radii for each pin (2D array, double, one row per pin, one less column than material)
%   pinmesh - Number of subregions to use for each region of each pin 
%             (2D array, integer, one row per pin, one column per material)
%   npol    - Number of polar angles to use for the ray (scalar, integer, range[1,1])
%   igroup  - Index of the energy group that is being solved (TODO: remove need for this)

%% Material IDs
nmats = 8;
id_mod = 1;
id_clad = 2;
id_gap = 3;
id_fuel = 4;
id_gtube = 5;
id_controlgap = 6;
id_controlmod = 7;
id_control = 8;

%% Cross-sections
ngroups = 47;
% TODO: Get rid of this, add real cross-sections
xstr_list = genXSData(nmats, ngroups, id_mod, id_clad, id_gap, id_fuel, id_gtube, id_controlgap, ...
    id_controlmod, id_control);
source_list = ones(nmats,ngroups);
source_list(id_controlmod,:) = source_list(id_mod,:);
source_list(id_gtube:id_control,47) = [7.58E-04; 1.50E-06; 2.44E-02; 8.44399E-05];
source_list(id_controlgap,:) = source_list(id_gap,:);
source_list(id_mod:id_fuel,47) = [0.136409169; 4.03E-03; 1.50E-06; 6.74E-03];

%% Quadrature
quad = quadratureClass(npol);

%% Mesh
mesh = meshClass(pinmap, pinmats, radii, pinmesh, pitch, diag);

%% Perform Sweep
mesh = setupFSP(source_list, xstr_list, mesh, igroup);
[angflux, ~, scalflux] = sweep(mesh, 0.0, quad);

end