classdef eigensolverClass
    %EIGENSOLVERCLASS Solves an eigenvalue problem
    %   This class contains the mesh, cross-section library,
    %   and quadrature required to solve a 1D MOC problem.
    %   It also contains methods to setup, solve, and update
    %   the solution, as well as check for convergence.
    
    properties
        xsLib
        mesh
        quad
        solution
        criteria
        input
        nouters
        converged
    end
    
    methods
        function obj = eigensolverClass(xsLib, mesh, quad, solution, nouters) %, criteria)
            %EIGENSOLVERCLASS Sets up the eigensolver
            %   xsLib    - The XS Library to use for the calculation
            %   mesh     - The mesh to solve the problem on
            %   quad     - The quadrature to use for the calculation
            %   solution - The solution object that contains solution data
            %   criteria - The convergence criteria object
            
            obj.xsLib = xsLib;
            obj.mesh = mesh;
            obj.quad = quad;
            obj.solution = solution;
%             obj.criteria = criteria;
            obj.nouters = nouters;
            obj.converged = false;
            
        end
        
        function obj = setup(obj)
            %SETUP Prepares the solver for the next iteration
            %   obj - The eigensolver object to set up
            
            obj.solution = obj.solution.calcFissSrc( obj.mesh, obj.xsLib );
            
        end
        
        function obj = solve(obj)
            %SOLVE Solves the eigenvalue problem
            %   obj - The eigensolver object to solve
            
            for iouter=1:obj.nouters
                display(sprintf('Eigenvalue iteration %i',iouter));
                obj.solution = obj.solution.update(iouter);
                obj = obj.step();
                obj = obj.update();
               if obj.converged
                   display(sprintf('Converged after %i iterations...',iouter));
                   break
               elseif iouter == obj.nouters
                   display(sprintf('Reached maximum number of iterations...'));
               end
            end
            
        end
        
        function obj = step(obj)
            %STEP Performs a single iteration of the eigenvalue solve
            %   obj - The eigensolver object to solve
            
            for igroup=1:obj.xsLib.ngroups
                obj.mesh = setupFSP(obj.solution, obj.xsLib, obj.mesh, igroup);
                obj.solution = sweep(igroup, obj.solution, obj.mesh, obj.quad);
            end
            
        end
        
        function obj = update(obj)
            %UPDATE Updates eigenvalue object after each iteration
            %   obj - The eigensolver object to update
            
            obj.solution = obj.solution.calcFissSrc( obj.mesh, obj.xsLib );
            obj.solution = obj.solution.updateEig( );
            [conv_flux, conv_keff] = obj.solution.calcResidual( obj.mesh, obj.xsLib);
            display(sprintf('Flux norm : %0.7f',conv_flux));
            display(sprintf('k-eff norm: %0.7f',conv_keff));
            display(sprintf('k-eff     : %0.7f\n',obj.solution.keff(1)));
            if abs(conv_flux) < 1.0e-7 && abs(conv_keff) < 1.0e-8
                obj.converged = true;
            end
            
        end
    end
    
end

