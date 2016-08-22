function [ solution, mesh ] = ...
    MOC_1D( pinmap, pitch, diag, pinmats, radii, pinmesh, npol, filename, scattype, BCond, nouters )
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
%   nouters  - Maximum number of outer iterations for an eigenvalue calculation

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
solution = solution.calcFissSrc( mesh, xsLib );
solution.fisssrc(:,1) = 1.0;

% Solve
for iouter=1:nouters
    display(sprintf('Eigenvalue iteration %i',iouter));
    solution = solution.update( iouter );
    for igroup=1:xsLib.ngroups
        mesh = setupFSP(solution, xsLib, mesh, igroup);
        solution = sweep(igroup, solution, mesh, quad);
    end
    solution = solution.calcFissSrc( mesh, xsLib );
    solution = solution.updateEig( );
    [conv_flux, conv_keff] = solution.calcResidual( mesh, xsLib );
    display(sprintf('Flux norm : %g',conv_flux));
    display(sprintf('k-eff norm: %g',conv_keff));
    display(sprintf('k-eff     : %g',solution.keff(1)));
    if abs(conv_flux) < 1.0e-7 && abs(conv_keff) < 1.0e-8
        display(sprintf('Converged after %i iterations...',iouter));
        break
    end
    display('  ');
end

end