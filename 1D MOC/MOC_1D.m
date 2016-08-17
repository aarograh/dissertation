function [ solution, mesh ] = ...
    MOC_1D( pinmap, pitch, diag, pinmats, radii, pinmesh, npol, filename, scattype, BCond )
%MOC_1D Solves 1D MOC given geometry and quadrature inputs
%   pinmap   - Map of pins in problem (vector, integer)
%   pitch    - Pitch for each pin (scalar, double)
%   diag     - Flag to indicate of pins are trace horizontally or vertically (0: horizontal, 1: diagonal)
%   pinmats  - List of materials in each pin (2D array, integer, one row per pin)
%   radii    - List of radii for each pin (2D array, double, one row per pin, one less column than material)
%   pinmesh  - Number of subregions to use for each region of each pin 
%              (2D array, integer, one row per pin, one column per material)
%   npol     - Number of polar angles to use for the ray (scalar, integer, range[1,1])
%   filename - Name of the XS Library file
%   scattype - Transport Scattering option.  Currently accepted values are P0.
%   BCond    - Boundary condition for MOC sweeps

%% Material IDs
nmats = 5;
id_mod = 1;
id_clad = 2;
id_gap = 3;
id_fuel = 4;
id_control = 5;

%% Cross-sections
ngroups = 47;
xsLib = xsLibraryClass( filename, scattype );
source_list = ones(nmats,ngroups);
source_list(id_mod:id_control,47) = [0.136409169; 4.03E-03; 1.50E-06; 6.74E-03; 8.44399E-05];

%% Quadrature
quad = quadratureClass(npol);

%% Mesh
mesh = meshClass(pinmap, pinmats, radii, pinmesh, pitch, diag);

%% Initialize Solution and Perform Sweeps
% Initialize solution
solution = solutionClass(mesh.nfsrcells,npol,xsLib.ngroups,BCond);

% Solve
solution = solution.update();
for igroup=1:xsLib.ngroups
    mesh = setupFSP(solution, source_list, xsLib, mesh, igroup);
    solution = sweep(solution, mesh, quad);
end

end