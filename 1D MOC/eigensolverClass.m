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
            obj.solution = solutionClass(obj.mesh.nfsrcells,obj.xsLib.ngroups,input);
            if ~isempty(input.cmfd)
                if input.cmfd
                    obj.cmfd = cmfdClass(input, obj.xsLib.ngroups);
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
            
            obj.solution = obj.solution.calcFissSrc(obj.mesh, obj.xsLib);
            for iouter=1:obj.nouters
                if obj.verbose
                    display(sprintf('Eigenvalue iteration %i',iouter));
                end
                if obj.accel
                    obj.solution = obj.cmfd.solve(obj.solution, obj.mesh);
                end
                obj = obj.step();
                obj = obj.update();
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
        
        function obj = step(obj, source_in)
            %STEP Performs a single iteration of the eigenvalue solve
            %   obj       - The eigensolver object to solve
            %   source_in - Flag to indicate if a user-specified source is present (true) or not (false).
            %               This is an optional argument
            
            obj.solution = obj.solution.update();
            if ~exist('source_in','var')
                source_in=false;
            end
            for igroup=1:obj.xsLib.ngroups
                if ~source_in && ~obj.accel
                    obj = obj.setupFSP(igroup);
                end
                obj = obj.sweep(igroup);
            end
            
        end
        
        function obj = update(obj)
            %UPDATE Updates eigenvalue object after each iteration
            %   obj - The eigensolver object to update
            
            obj.solution = obj.solution.calcFissSrc( obj.mesh, obj.xsLib );
            if ~obj.accel
                obj.solution = obj.solution.updateEig( );
            end
            [conv_flux, conv_keff] = obj.solution.calcResidual( obj.mesh, obj.xsLib);
            if obj.verbose
                display(sprintf('Flux norm : %0.7f',conv_flux));
                display(sprintf('k-eff norm: %0.7f',conv_keff));
                display(sprintf('k-eff     : %0.7f\n',obj.solution.keff(1)));
            end
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
            
            obj.mesh.source(1:obj.mesh.nfsrcells,igroup) = 0.0;
            for i=1:obj.mesh.nfsrcells
                % Use old scalar flux to do Jacobi style iteration
                matID = obj.mesh.materials(i);
                obj.mesh.source(i,igroup) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).chi(igroup)/obj.solution.keff(1) + ...
                    obj.solution.scalflux(i,:,2)*obj.xsLib.xsSets(matID).scatter(igroup,:)')*0.5;
                obj.mesh.xstr(i,igroup) = obj.xsLib.xsSets(matID).transport(igroup);
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
                    exparg = exp(-obj.mesh.xstr(i,igroup)*dx);
                    obj.solution.angflux(i+1,j,1,igroup) = obj.solution.angflux(i,j,1,igroup)*exparg + ...
                        obj.mesh.source(i,igroup)/obj.mesh.xstr(i,igroup)*(1 - exparg);
                    obj.solution.scalflux(i,igroup,1) = obj.solution.scalflux(i,igroup,1) + ...
                        0.5*sum(obj.solution.angflux(i:i+1,j,1,igroup))*obj.quad.weights(j);
                    
                    % Backward Sweep
                    dx = (obj.mesh.fsredges(k+1) - obj.mesh.fsredges(k))/obj.quad.cosines(j);
                    exparg = exp(-obj.mesh.xstr(k,igroup)*dx);
                    obj.solution.angflux(k,j,2,igroup) = obj.solution.angflux(k+1,j,2,igroup)*exparg + ...
                        obj.mesh.source(k,igroup)/obj.mesh.xstr(k,igroup)*(1 - exparg);
                    obj.solution.scalflux(k,igroup,1) = obj.solution.scalflux(k,igroup,1) + ...
                        0.5*sum(obj.solution.angflux(k:k+1,j,2,igroup))*obj.quad.weights(j);
                end
            end
        end
    end
    
end

