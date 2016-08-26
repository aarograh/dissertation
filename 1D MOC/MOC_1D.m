function solver = MOC_1D( input )
%MOC_1D Solves 1D MOC given geometry and quadrature inputs
%   input - The inputClass object to use while setting up the problem

%% Initialize Eigensolver
solver = eigensolverClass(input); %, critera);

%% Solve
solver = solver.setup();
solver = solver.solve();

end