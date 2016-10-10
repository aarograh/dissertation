classdef eigensolverClass < handle
    %EIGENSOLVERCLASS Solves an eigenvalue problem
    %   This class contains the mesh, cross-section library,
    %   and quadrature required to solve a 1D MOC problem.
    %   It also contains methods to setup, solve, and checkConv
    %   the solution, as well as check for convergence.
    
    properties
        xsLib
        mesh
        quad
        fss
        accel
        cmfd
        conv_crit=[1.0e-8, 1.0e-7];
        input
        nouters=1000
        converged
        verbose=true
    end
    
    methods
        function obj = eigensolverClass( input )
            %EIGENSOLVERCLASS Sets up the eigensolver
            %   input - The inputClass container from which to initialize
            
            obj.xsLib = xsLibraryClass(input);
            obj.mesh = meshClass(input);
            obj.quad = quadratureClass(input);
            obj.fss = FixedSourceSolverClass(obj.xsLib, obj.quad, input);
            if ~isempty(input.cmfd)
                if input.cmfd
                    obj.cmfd = cmfdClass(input, obj.mesh, obj.xsLib);
                    obj.accel = true;
                else
                    obj.accel = false;
                end
            else
                obj.accel = false;
            end
            if input.conv_crit(1) > 0
                obj.conv_crit(1) = input.conv_crit(1);
            end
            if input.conv_crit(2) > 0
                obj.conv_crit(2) = input.conv_crit(2);
            end
            obj.nouters = input.nouters;
            obj.converged = false;
            if ~isempty(input.verbose)
                obj.verbose = input.verbose;
            end
            
        end
        
        function obj = solve(obj)
            %SOLVE Solves the eigenvalue problem
            %   obj - The eigensolver object to solve
            
            obj.mesh.source(1:20,1) = 0.0;
            for iouter=1:obj.nouters
                if obj.verbose
                    display(sprintf('Eigenvalue iteration %i',iouter));
                end
                if obj.accel && iouter < 180
                    obj.cmfd.solve(obj.fss.solution, obj.mesh);
                end
                obj.step();
                if obj.converged
                    if obj.verbose
                        display(sprintf('Converged after %i iterations...',iouter));
                    end
                    break
                elseif iouter == obj.nouters
                    if obj.verbose
                        display(sprintf('Reached maximum number of iterations...'));
                    end
                end
            end
            
        end
        
        function obj = step(obj)
            %STEP Performs a single iteration of the eigenvalue solve
            %   obj       - The eigensolver object to solve
            
            if obj.accel
                obj.fss.solve(0, 1);
            else
                obj.fss.solution.updateEig( );
                obj.fss.solve(1, 1);
            end
            
            obj.checkConv();
            
        end
        
        function obj = checkConv(obj)
            %CHECKCONV Updates eigenvalue object after each iteration
            %   obj - The eigensolver object whose convergence should be checked
            
            [conv_flux, conv_keff] = obj.fss.solution.calcResidual( obj.mesh, obj.xsLib);
            if obj.verbose
                display(sprintf('Flux norm : %0.7f',conv_flux));
                display(sprintf('k-eff norm: %0.7f',conv_keff));
                display(sprintf('k-eff     : %0.7f\n',obj.fss.solution.keff(1)));
            end
            if abs(conv_flux) < obj.conv_crit(2) && abs(conv_keff) < obj.conv_crit(1)
                obj.converged = true;
            end
            
        end
    end
    
end

