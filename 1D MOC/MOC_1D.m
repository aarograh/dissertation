function solver = MOC_1D( input )
%MOC_1D Solves 1D MOC given geometry and quadrature inputs
%   input - The inputClass object to use while setting up the problem

%% Cross-sections
display('Setting up XS Library...')
xsLib = xsLibraryClass( input.xsfilename, input.scattype );

%% Quadrature
display('Setting up Quadrature...')
quad = quadratureClass(input.npol);

%% Mesh
display('Setting up Mesh...')
mesh = meshClass(input.pinmap, input.pinmats, input.radii, input.pinmesh, input.pitch, input.diag);

%% Initialize Solution and Perform Sweeps
% Initialize solution
display('Initializing Solution...')
solution = solutionClass(mesh.nfsrcells,quad.npol,xsLib.ngroups,input.BCond);

%% Initialize Eigensolver
solver = eigensolverClass(xsLib, mesh, quad, solution, input.nouters); %, critera);

%% Solve
solver = solver.setup();
solver = solver.solve();

end