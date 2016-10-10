classdef FixedSourceSolverClass < handle
    %FIXEDSOURCESOLVERCLASS Performs a fixed source calculation
    %   Uses a fixed fission source to perform a calculation for
    %   the scalar flux.
    
    properties
        xsLib
        mesh
        quad
        solution
    end
    
    methods
        function obj = FixedSourceSolverClass( xsLib, quad, input )
            %FIXEDSOURCESOLVERCLASS Initializes a fixedSourceSolverClass object
            %   xsLib    - XSLibraryClass object
            %   quad     - QuadratureClass object
            
            obj.xsLib = xsLib;
            obj.mesh = meshClass(input);
            obj.quad = quad;
            obj.solution = solutionClass(obj.mesh.nfsrcells, obj.xsLib.ngroups, input);
        end
        
        function obj = initFromEigenSolver( eig )
            %INITFROMEIGENSOLVER Initializes a fixedSourceSolverClass object from
            %an eigensolverClass object
            %   eig - EigensolverClass object to use for initialization
            
            
            
        end
        
        function obj = solve( obj, wCur, ninners )
            %SOLVE Solves a fixed source problem
            %   obj     - The FixedSourceSolverClass object to solve
            %   wCur    - Logical to tally currents (1) or not (0)
            %   ninners - The number of inner iterations to perform
            
            obj.solution.calcFissSrc(obj.mesh, obj.xsLib);
            obj.solution.updateBC();
            for inner=1:ninners
                for igroup=1:obj.xsLib.ngroups
                    obj.setup(igroup);
                    if wCur
                        obj.sweep_wCur(igroup);
                    else
                        obj.sweep(igroup);
                    end
                end
            end
            
            obj.solution.calcFissSrc(obj.mesh, obj.xsLib);
        end
        
        function obj = setup( obj, igroup )
            %SETUPFSP Sets up source and XS mesh for fixed source MOC problem
            %   obj    - The FixedSourceSolverClass object to set up
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

