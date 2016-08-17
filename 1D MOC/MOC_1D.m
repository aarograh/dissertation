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

%% Cross-sections
display('Setting up XS Library...')
xsLib = xsLibraryClass( filename, scattype );

%% Quadrature
display('Setting up Quadrature...')
quad = quadratureClass(npol);

%% Mesh
display('Setting up Mesh...')
mesh = meshClass(pinmap, pinmats, radii, pinmesh, pitch, diag);

%% Initialize Solution and Perform Sweeps
% Initialize solution
display('Initializing Solution...')
solution = solutionClass(mesh.nfsrcells,npol,xsLib.ngroups,BCond);

% Solve
for igroup=1:xsLib.ngroups
    display(sprintf('Sweeping Energy Group %i',igroup))
    mesh = setupFSP(solution, xsLib, mesh, igroup);
    solution = sweep(igroup, solution, mesh, quad);
    solution = solution.update();
end

end