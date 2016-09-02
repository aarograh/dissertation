classdef eigensolverClass < handle
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
                if iouter == 1
                    tmp = obj.mesh.source(:,1);
                end
%                 display([obj.solution.fisssrc(:,1),obj.solution.fisssrc(:,2), ...
%                     obj.mesh.source(:,1),obj.solution.scalflux(:,1,1),obj.solution.scalflux(:,1,2),...
%                     tmp]);
                if obj.accel && iouter < 180
                    obj.cmfd.solve(obj.solution, obj.mesh);
                else
                    obj.solution.calcFissSrc( obj.mesh, obj.xsLib );
                    obj.solution.updateEig( );
                end
                obj.solution.calcFissSrc(obj.mesh, obj.xsLib);
                obj.solution.update();
                obj.step();
                display([obj.solution.fisssrc(:,1),obj.solution.fisssrc(:,2), ...
                    obj.mesh.source(:,1),obj.solution.scalflux(:,1,1),obj.solution.scalflux(:,1,2),...
                    tmp]);
                tmp = obj.mesh.source(:,1);
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
            
            if ~exist('source_in','var')
                source_in=false;
            end
            for inners=1:1
                for igroup=1:obj.xsLib.ngroups
                    if ~source_in
                        obj.setupFSP(igroup);
                    end
                    if obj.accel
                        obj.sweep_wCur(igroup);
                    else
                        obj.sweep(igroup);
                    end
                end
                obj.update();
            end
            
        end
        
        function obj = update(obj)
            %UPDATE Updates eigenvalue object after each iteration
            %   obj - The eigensolver object to update
            
            if ~obj.accel
                obj.solution.calcFissSrc( obj.mesh, obj.xsLib );
                obj.solution.updateEig( );
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
            obj.solution.update();
            
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
                    obj.solution.scalflux(i,:,1)*obj.xsLib.xsSets(matID).scatter(igroup,:)')*0.5;
                obj.mesh.xstr(i,igroup) = obj.xsLib.xsSets(matID).transport(igroup);
            end
            
        end
        
        function obj = sweep( obj, igroup )
            %SWEEP Performs 1D MOC sweep for a single ray with multiple polars
            %   obj    - The eigensolver object to sweep
            %   igroup - The energy group being swept
            
            obj.solution.scalflux(:,igroup,2) = obj.solution.scalflux(:,igroup,1);
            obj.solution.scalflux(:,igroup,1) = 0.0;
            for i=1:obj.mesh.nfsrcells
                k = obj.mesh.nfsrcells-i+1;
                for j=1:obj.quad.npol
                    % Forward Sweep
                    dx = (obj.mesh.fsredges(i+1) - obj.mesh.fsredges(i))/obj.quad.cosines(j);
                    exparg = exp(-obj.mesh.xstr(i,igroup)*dx);
                    obj.solution.angflux(i+1,j,1,igroup) = obj.solution.angflux(i,j,1,igroup)*exparg + ...
                        obj.mesh.source(i,igroup)/obj.mesh.xstr(i,igroup)*(1 - exparg);
                    
                    psibar = 0.5*sum(obj.solution.angflux(i:i+1,j,1,igroup));
                    obj.solution.scalflux(i,igroup,1) = obj.solution.scalflux(i,igroup,1) + ...
                        psibar*obj.quad.weights(j);
                    
                    % Backward Sweep
                    dx = (obj.mesh.fsredges(k+1) - obj.mesh.fsredges(k))/obj.quad.cosines(j);
                    exparg = exp(-obj.mesh.xstr(k,igroup)*dx);
                    obj.solution.angflux(k,j,2,igroup) = obj.solution.angflux(k+1,j,2,igroup)*exparg + ...
                        obj.mesh.source(k,igroup)/obj.mesh.xstr(k,igroup)*(1 - exparg);
                    
                    psibar = 0.5*sum(obj.solution.angflux(k:k+1,j,2,igroup));
                    obj.solution.scalflux(k,igroup,1) = obj.solution.scalflux(k,igroup,1) + ...
                        psibar*obj.quad.weights(j);
                end
            end
        end
        
        function obj = sweep_wCur( obj, igroup )
            %SWEEP_wCur Performs 1D MOC sweep for a single ray with multiple polars while
            %           tallying currents for cmfd acceleration
            %   obj    - The eigensolver object to sweep
            %   igroup - The energy group being swept
            
            obj.solution.scalflux(:,igroup,2) = obj.solution.scalflux(:,igroup,1);
            obj.solution.current(:,igroup,2) = obj.solution.current(:,igroup,1);
            obj.solution.scalflux(:,igroup,1) = 0.0;
            obj.solution.current(:,igroup,1) = 0.0;
            for j = 1:obj.quad.npol
                obj.solution.current(1,igroup,1) = obj.solution.current(1,igroup,1) + ...
                    obj.solution.angflux(1,j,1,igroup)*obj.quad.cosines(j)*obj.quad.weights(j);
                obj.solution.current(end,igroup,1) = obj.solution.current(end,igroup,1) - ...
                    obj.solution.angflux(end,j,2,igroup)*obj.quad.cosines(j)*obj.quad.weights(j);
            end
            
            for i=1:obj.mesh.nfsrcells
                k = obj.mesh.nfsrcells-i+1;
                for j=1:obj.quad.npol
                    % Forward Sweep
                    dx = (obj.mesh.fsredges(i+1) - obj.mesh.fsredges(i))/obj.quad.cosines(j);
                    exparg = exp(-obj.mesh.xstr(i,igroup)*dx);
                    obj.solution.angflux(i+1,j,1,igroup) = obj.solution.angflux(i,j,1,igroup)*exparg + ...
                        obj.mesh.source(i,igroup)/obj.mesh.xstr(i,igroup)*(1 - exparg);
                    
                    psibar = 0.5*sum(obj.solution.angflux(i:i+1,j,1,igroup));
                    obj.solution.scalflux(i,igroup,1) = obj.solution.scalflux(i,igroup,1) + ...
                        psibar*obj.quad.weights(j);
                    obj.solution.current(i+1,igroup,1) = obj.solution.current(i+1,igroup,1) + ...
                        obj.solution.angflux(i+1,j,1,igroup)*obj.quad.cosines(j)*obj.quad.weights(j);
                    
                    % Backward Sweep
                    dx = (obj.mesh.fsredges(k+1) - obj.mesh.fsredges(k))/obj.quad.cosines(j);
                    exparg = exp(-obj.mesh.xstr(k,igroup)*dx);
                    obj.solution.angflux(k,j,2,igroup) = obj.solution.angflux(k+1,j,2,igroup)*exparg + ...
                        obj.mesh.source(k,igroup)/obj.mesh.xstr(k,igroup)*(1 - exparg);
                    
                    psibar = 0.5*sum(obj.solution.angflux(k:k+1,j,2,igroup));
                    obj.solution.scalflux(k,igroup,1) = obj.solution.scalflux(k,igroup,1) + ...
                        psibar*obj.quad.weights(j);
                    obj.solution.current(k,igroup,1) = obj.solution.current(k,igroup,1) - ...
                        obj.solution.angflux(k,j,2,igroup)*obj.quad.cosines(j)*obj.quad.weights(j);
                end
            end
        end
    end
    
end

