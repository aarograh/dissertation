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
        conv_crit=[1.0e-8, 1.0e-7];
        input
        nouters=1000
        converged
    end
    
    methods
        function obj = eigensolverClass( input )
            %EIGENSOLVERCLASS Sets up the eigensolver
            %   input - The inputClass container from which to initialize
            
            obj.xsLib = xsLibraryClass(input);
            obj.mesh = meshClass(input);
            obj.quad = quadratureClass(input);
            obj.solution = solutionClass(obj.mesh.nfsrcells,obj.xsLib.ngroups,input);
            if input.conv_crit(1) > 0
              obj.conv_crit(1) = input.conv_crit(1);
            end
            if input.conv_crit(2) > 0
              obj.conv_crit(2) = input.conv_crit(2);
            end
            obj.nouters = input.nouters;
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
                obj = obj.setupFSP(igroup);
                obj = obj.sweep(igroup);
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
            if abs(conv_flux) < obj.conv_crit(2) && abs(conv_keff) < obj.conv_crit(1)
                obj.converged = true;
            end
            
        end
        
        function obj = setupFSP( obj, igroup )
            %SETUPFSP Sets up source and XS mesh for fixed source MOC problem
            %   solution    - The solution data to use to determine the source
            %   xsLib       - Cross-section library object
            %   mesh        - The mesh for this problem
            %   igroup      - Group index
            
            obj.mesh.source(1:obj.mesh.nfsrcells,1) = 0.0;
            for i=1:obj.mesh.nfsrcells
                % Use old scalar flux to do Jacobi style iteration
                matID = obj.mesh.materials(i);
                obj.mesh.source(i,1) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).chi(igroup)/obj.solution.keff(1) + ...
                    obj.solution.scalflux(i,:,2)*obj.xsLib.xsSets(matID).scatter(igroup,:)')*0.5;
                obj.mesh.xstr(i,1) = obj.xsLib.xsSets(matID).transport(igroup);
            end
            
        end
        
        function obj = sweep( obj, igroup )
            %SWEEP Performs 1D MOC sweep for a single ray with multiple polars
            %   obj    - The eigensolver object to sweep
            %   igroup - The energy group being swept
            
            for i=1:obj.mesh.nfsrcells
                k = obj.mesh.nfsrcells-i+1;
                for j=1:obj.quad.npol
                    % Forward Sweep
                    dx = (obj.mesh.fsredges(i+1) - obj.mesh.fsredges(i))/obj.quad.cosines(j);
                    exparg = exp(-obj.mesh.xstr(i)*dx);
                    obj.solution.angflux(i+1,j,1,igroup) = obj.solution.angflux(i,j,1,igroup)*exparg + ...
                        obj.mesh.source(i)/obj.mesh.xstr(i)*(1 - exparg);
                    obj.solution.scalflux(i,igroup,1) = obj.solution.scalflux(i,igroup,1) + ...
                        0.5*sum(obj.solution.angflux(i:i+1,j,1,igroup))*obj.quad.weights(j);
                    
                    % Backward Sweep
                    dx = (obj.mesh.fsredges(k+1) - obj.mesh.fsredges(k))/obj.quad.cosines(j);
                    exparg = exp(-obj.mesh.xstr(k)*dx);
                    obj.solution.angflux(k,j,2,igroup) = obj.solution.angflux(k+1,j,2,igroup)*exparg + ...
                        obj.mesh.source(k)/obj.mesh.xstr(k)*(1 - exparg);
                    obj.solution.scalflux(k,igroup,1) = obj.solution.scalflux(k,igroup,1) + ...
                        0.5*sum(obj.solution.angflux(k:k+1,j,2,igroup))*obj.quad.weights(j);
                end
            end            
        end
    end
    
end

