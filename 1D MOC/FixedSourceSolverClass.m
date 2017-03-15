classdef FixedSourceSolverClass < handle
    %FIXEDSOURCESOLVERCLASS Performs a fixed source calculation
    %   Uses a fixed fission source to perform a calculation for
    %   the scalar flux.
    
    properties
        xsLib
        mesh
        quad
        solution
        converged
        verbose=false
        relax=0.25
        subray=false
        submesh_vol
        homog=0 %0 mixes subray scalar fluxes, 1 mixes cross-sections and re-solves
    end
    
    methods
        function obj = FixedSourceSolverClass( varargin )
            %FIXEDSOURCESOLVERCLASS Initializes a fixedSourceSolverClass object
            %   varargin - Input arguments.  These can be "xsLib, quad, input" or
            %              "input, eig".  Anything else results in an error.
            %   xsLib - XSLibraryClass object
            %   quad  - QuadratureClass object
            %   input - inputClass object to initialize from
            %   eig   - eigensolverClass object to initialize from (optional; other arguments
            %           ignored if this one is present)
            
            if nargin == 2
                obj.mesh = meshClass(varargin{1});
                obj.xsLib = varargin{2}.xsLib;
                obj.quad = varargin{2}.quad;
                obj.solution = solutionClass(obj.mesh, obj.xsLib, varargin{1});
                obj.initFromEigenSolver( varargin{2});
                obj.verbose = varargin{1}.verbose;
                obj.subray = varargin{1}.subray;
                obj.homog = varargin{1}.homog;
            elseif nargin == 3
                obj.xsLib = varargin{1};
                obj.quad = varargin{2};
                obj.mesh = meshClass(varargin{3});
                obj.solution = solutionClass(obj.mesh, obj.xsLib, varargin{3});
                obj.verbose = varargin{3}.verbose;
                obj.subray = varargin{3}.subray;
                obj.homog = varargin{3}.homog;
            end
            if obj.subray
                for i=1:obj.mesh.nfsrcells
                    if obj.xsLib.xsSets(obj.mesh.materials(i)).nsubxs > 0
                        break
                    end
                end
                obj.submesh_vol = obj.xsLib.xsSets(obj.mesh.materials(i)).subfracs;
                if obj.homog == 1
                    i = length(obj.solution.angflux(1,1,1,1,:))+1;
                    obj.solution.angflux(:,:,:,:,i) = 0.0;
                    obj.solution.current(:,:,i,:) = 0.0;
                end
            end
        end
        
        function obj = initFromEigenSolver( obj, eig )
            %INITFROMEIGENSOLVER Initializes a fixedSourceSolverClass object from
            %an eigensolverClass object
            %   obj - fixedSourceSolverClass object to initialize
            %   eig - EigensolverClass object to use for initialization
            
            obj.solution.keff(:) = eig.fss.solution.keff(:);
            obj.solution.scalflux(:) = eig.fss.solution.scalflux(:);
            if (length(obj.solution.angflux(1,1,1,1,:)) == 1)
                obj.solution.angflux(:,:,:,:,1)=eig.fss.solution.angflux(:,:,:,:,1);
            elseif (length(eig.fss.solution.angflux(1,1,1,1,:)) == 1)
                for i=1:length(obj.solution.angflux(1,1,1,1,:))
                    obj.solution.angflux(:,:,:,:,i) = eig.fss.solution.angflux(:,:,:,:,1);
                end
            else
                obj.solution.angflux(:) = eig.fss.solution.angflux(:);
            end
            obj.solution.fisssrc(:) = eig.fss.solution.fisssrc(:);
            
        end
        
        function obj = solve( obj, wCur, ninners, convcrit )
            %SOLVE Solves a fixed source problem
            %   obj      - The FixedSourceSolverClass object to solve
            %   wCur     - Logical to tally currents (1) or not (0)
            %   ninners  - The number of inner iterations to perform
            %   convcrit - Convergence criteria for scattering source calculation
            
            obj.solution.updateBC();
            obj.converged = false;
            inner=0;
            if ~exist('convcrit','var')
                convcrit=1.0e-3;
            end
            
            while ~obj.converged
                inner = inner + 1;
                obj.solution.scalflux(:,:,2) = obj.solution.scalflux(:,:,1);
                if wCur && inner == ninners
                    obj.setup( 0 );
                    obj.sweep_wCur( );
                    obj.solution.scalflux(:,:,1) = obj.relax*obj.solution.scalflux(:,:,1) + ...
                        (1.0-obj.relax)*obj.solution.scalflux(:,:,2);
                    obj.solution.current(:,:,:,1) = obj.relax*obj.solution.current(:,:,:,1) + ...
                        (1.0-obj.relax)*obj.solution.current(:,:,:,2);
                else
                    if obj.subray
                        for i=1:length(obj.submesh_vol)
                            obj.setup( i );
                            obj.sweep( i );
                        end
                        if obj.homog == 1
                            obj.setup_subrayHomog( );
                            nsubmesh = length(obj.submesh_vol)+1;
                            obj.sweep( nsubmesh );
                            obj.postprocess_unified( nsubmesh );
                        else
                            obj.postprocess_subray( );
                        end
                    else
                        obj.setup( 0 );
                        obj.sweep( 0 );
                        obj.postprocess_unified( );
                    end
                end
                scatconv = obj.checkConv( );
                if obj.verbose
                    display(sprintf('Fixed Source Iteration %i - %g',inner,scatconv));
                end
                if inner == ninners || scatconv < convcrit
                    obj.converged = true;
                    if obj.verbose
                        display('Fixed Source Solve Converged.')
                    end
                end
            end
            
        end
        
        function obj = setup( obj, level )
            %SETUP Sets up source and XS mesh for fixed source MOC problem
            %   obj    - The FixedSourceSolverClass object to set up
            %   level  - The submesh level being set up
            
            obj.mesh.source(:) = 0.0;
            if level == 0
                for i=1:obj.mesh.nfsrcells
                    matID = obj.mesh.materials(i);
                    for j=1:obj.xsLib.ngroups
                        % Use old scalar flux to do Jacobi style iteration
                        obj.mesh.source(j,i) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).chi(j)/obj.solution.keff(1) + ...
                            obj.xsLib.xsSets(matID).scatter(j,:)*obj.solution.scalflux(:,i,2))*0.5;
                        obj.mesh.xstr(j,i) = obj.xsLib.xsSets(matID).transport(j);
                    end
                end
            else
                for i=1:obj.mesh.nfsrcells
                    matID = obj.mesh.materials(i);
                    for j=1:obj.xsLib.ngroups
                        % Use old scalar flux to do Jacobi style iteration
                        if obj.xsLib.xsSets(matID).nsubxs > 0
                            obj.mesh.source(j,i) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).subxs(level).chi(j)/obj.solution.keff(1) + ...
                                obj.xsLib.xsSets(matID).subxs(level).scatter(j,:)*obj.solution.scalflux(:,i,2))*0.5;
                            obj.mesh.xstr(j,i) = obj.xsLib.xsSets(matID).subxs(level).transport(j);
                        else
                            obj.mesh.source(j,i) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).chi(j)/obj.solution.keff(1) + ...
                                obj.xsLib.xsSets(matID).scatter(j,:)*obj.solution.scalflux(:,i,2))*0.5;
                            obj.mesh.xstr(j,i) = obj.xsLib.xsSets(matID).transport(j);
                        end
                    end
                end
            end
            
        end
        
        function obj = setup_subrayHomog( obj )
            %SETUP_SUBRAYHOMOG Sets up source and XS mesh for final fixed source MOC problem after subray
            %sweeps
            %   obj - The FixedSourceSolverClass object to set up
            
            obj.mesh.source(:) = 0.0;
            nsubmesh = length(obj.submesh_vol);
            
            for k=1:nsubmesh
                obj.postprocess_unified( 1 );
                scalflux = obj.solution.scalflux(:,:,1)*obj.submesh_vol(k);
                fisssrc = obj.solution.fisssrc(:,1)*obj.submesh_vol(k);
                for i=1:obj.mesh.nfsrcells
                    matID = obj.mesh.materials(i);
                    for j=1:obj.xsLib.ngroups
                        if obj.xsLib.xsSets(matID).nsubxs > 0
                            obj.mesh.source(j,i) = obj.mesh.source(j,i) + ...
                                (fisssrc(i)*obj.xsLib.xsSets(matID).subxs(k).chi(j)/obj.solution.keff(1) + ...
                                obj.xsLib.xsSets(matID).subxs(k).scatter(j,:)*scalflux(:,i))*0.5;
                            obj.mesh.xstr(j,i) = obj.mesh.xstr(j,i) + ...
                                fisssrc(i)*obj.xsLib.xsSets(matID).subxs(k).transport(j)*scalflux(j,i);
                        else
                            obj.mesh.source(j,i) = obj.mesh.source(j,i) + ...
                                (fisssrc(i)*obj.xsLib.xsSets(matID).chi(j)/obj.solution.keff(1) + ...
                                obj.xsLib.xsSets(matID).scatter(j,:)*scalflux(:,i))*0.5;
                            obj.mesh.xstr(j,i) = obj.mesh.xstr(j,i) + ...
                                fisssrc(i)*obj.xsLib.xsSets(matID).transport(j)*scalflux(j,i);
                        end
                    end
                end
            end
            
            for i=1:obj.mesh.nfsrcells
                obj.mesh.xstr(:,i) = obj.mesh.xstr(:,i)./scalflux(:,i);
            end
            
        end
        
        function maxdiff = checkConv( obj )
            %CHECKCONV Checks for convergence of the scattering source
            %   obj - FixedSourceSolverClass object whose convergence is being checked
            
            maxdiff = 0.0;
            for igroup=1:obj.xsLib.ngroups
                oldsource = zeros(obj.mesh.nfsrcells,1);
                newsource = zeros(obj.mesh.nfsrcells,1);
                for i=1:obj.mesh.nfsrcells
                    matID = obj.mesh.materials(i);
                    oldsource(i) = obj.xsLib.xsSets(matID).scatter(igroup,:)*obj.solution.scalflux(:,i,2);
                    newsource(i) = obj.xsLib.xsSets(matID).scatter(igroup,:)*obj.solution.scalflux(:,i,1);
                end
                maxdiff = max(max(abs((oldsource - newsource)./oldsource)),maxdiff);
            end
            
        end
        
        function obj = postprocess_unified( obj, level )
            %POSTPROCESS_UNIFIED Post-rpocess the sweep result for regular MOC
            %   obj    - The fixedsourcesolver object to post-process
            %   level  - The submesh level to use
            
            if ~exist('level','var')
                level = 1;
            end
            obj.solution.scalflux(:,:,1) = 0.0;
            for j=1:obj.mesh.nfsrcells
                for k=1:obj.quad.npol
                    for g=1:obj.xsLib.ngroups
                        psibar = sum(obj.solution.angflux(1,g,k,j:j+1,level));
                        psibar = 0.5*(psibar + sum(obj.solution.angflux(2,g,k,j:j+1,level)));
                        obj.solution.scalflux(g,j,1) = obj.solution.scalflux(g,j,1) + ...
                            psibar*obj.quad.weights(k);
                    end
                end
            end
            
        end
        
        function obj = postprocess_subray( obj )
            %POSTPROCESS_SUBRAY Post-processes the sweep result for sub-ray MOC
            %   obj    - The fixedsourcesolver object to post-process
            
            obj.solution.scalflux(:,:,1) = 0.0;
            for j=1:obj.mesh.nfsrcells
                for k=1:obj.quad.npol
                    for g=1:obj.xsLib.ngroups
                        psibar = 0.0;
                        for i=1:length(obj.submesh_vol)
                            psibar = psibar+obj.submesh_vol(i)*...
                                sum(obj.solution.angflux(1,g,k,j:j+1,i));
                            psibar = psibar+obj.submesh_vol(i)*...
                                sum(obj.solution.angflux(2,g,k,j:j+1,i));
                        end
                        psibar = psibar*0.5;
                        obj.solution.scalflux(g,j,1) = obj.solution.scalflux(g,j,1) + ...
                            psibar*obj.quad.weights(k);
                    end
                end
           end
            
        end
        
        function obj = sweep( obj, isubmesh )
            %SWEEP Performs 1D MOC sweep for a single ray with multiple polars
            %   obj      - The fixedsourcesolver object to sweep
            %   isubmesh - The submesh level to sweep
            
            if isubmesh == 0
                isubmesh = 1;
            end
            for i=1:obj.mesh.nfsrcells
                k = obj.mesh.nfsrcells-i+1;
                for j=1:obj.quad.npol
                    dx1 = (obj.mesh.fsredges(i+1) - obj.mesh.fsredges(i))/obj.quad.cosines(j);
                    dx2 = (obj.mesh.fsredges(k+1) - obj.mesh.fsredges(k))/obj.quad.cosines(j);
                    for igroup=1:obj.xsLib.ngroups
                        % Forward Sweep
                        exparg = exp(-obj.mesh.xstr(igroup,i)*dx1);
                        obj.solution.angflux(1,igroup,j,i+1,isubmesh) = ...
                            obj.solution.angflux(1,igroup,j,i,isubmesh)*exparg + ...
                            obj.mesh.source(igroup,i)/obj.mesh.xstr(igroup,i)*(1 - exparg);

                        % Backward Sweep
                        exparg = exp(-obj.mesh.xstr(igroup,k)*dx2);
                        obj.solution.angflux(2,igroup,j,k,isubmesh) = ...
                            obj.solution.angflux(2,igroup,j,k+1,isubmesh)*exparg + ...
                            obj.mesh.source(igroup,k)/obj.mesh.xstr(igroup,k)*(1 - exparg);
                    end
                end
            end
        end
        
        function obj = sweep_wCur( obj )
            %SWEEP_wCur Performs 1D MOC sweep for a single ray with multiple polars while
            %           tallying currents for cmfd acceleration
            %   obj    - The eigensolver object to sweep
            
            isubmesh = 1;
            obj.solution.scalflux(:,:,2) = obj.solution.scalflux(:,:,1);
            obj.solution.current(:,:,2) = obj.solution.current(:,:,1);
            obj.solution.scalflux(:,:,1) = 0.0;
            obj.solution.current(:,:,1) = 0.0;
            for j = 1:obj.quad.npol
                for igroup=1:obj.xsLib.ngroups
                    obj.solution.current(igroup,1,isubmesh,1) = obj.solution.current(igroup,1,isubmesh,1) + ...
                        obj.solution.angflux(1,igroup,j,1,isubmesh)*obj.quad.cosines(j)*obj.quad.weights(j);
                    obj.solution.current(igroup,end,isubmesh,1) = obj.solution.current(igroup,end,isubmesh,1) - ...
                        obj.solution.angflux(2,igroup,j,end,isubmesh)*obj.quad.cosines(j)*obj.quad.weights(j);
                end
            end

            for i=1:obj.mesh.nfsrcells
                k = obj.mesh.nfsrcells-i+1;
                for j=1:obj.quad.npol
                    dx1 = (obj.mesh.fsredges(i+1) - obj.mesh.fsredges(i))/obj.quad.cosines(j);
                    dx2 = (obj.mesh.fsredges(k+1) - obj.mesh.fsredges(k))/obj.quad.cosines(j);
                    for igroup=1:obj.xsLib.ngroups
                        % Forward Sweep
                        exparg = exp(-obj.mesh.xstr(igroup,i)*dx1);
                        obj.solution.angflux(1,igroup,j,i+1,isubmesh) = ...
                            obj.solution.angflux(1,igroup,j,i,isubmesh)*exparg + ...
                            obj.mesh.source(igroup,i)/obj.mesh.xstr(igroup,i)*(1 - exparg);

                        psibar = 0.5*sum(obj.solution.angflux(1,igroup,j,i:i+1,isubmesh));
                        obj.solution.scalflux(igroup,i,1) = obj.solution.scalflux(igroup,i,1) + ...
                            psibar*obj.quad.weights(j);
                        obj.solution.current(igroup,i+1,isubmesh,1) = obj.solution.current(igroup,i+1,isubmesh,1) + ...
                            obj.solution.angflux(1,igroup,j,i+1,isubmesh)*obj.quad.cosines(j)*obj.quad.weights(j);

                        % Backward Sweep
                        exparg = exp(-obj.mesh.xstr(igroup,k)*dx2);
                        obj.solution.angflux(2,igroup,j,k,isubmesh) = ...
                            obj.solution.angflux(2,igroup,j,k+1,isubmesh)*exparg + ...
                            obj.mesh.source(igroup,k)/obj.mesh.xstr(igroup,k)*(1 - exparg);

                        psibar = 0.5*sum(obj.solution.angflux(2,igroup,j,k:k+1,isubmesh));
                        obj.solution.scalflux(igroup,k,1) = obj.solution.scalflux(igroup,k,1) + ...
                            psibar*obj.quad.weights(j);
                        obj.solution.current(igroup,k,isubmesh,1) = obj.solution.current(igroup,k,isubmesh,1) - ...
                            obj.solution.angflux(2,igroup,j,k,isubmesh)*obj.quad.cosines(j)*obj.quad.weights(j);
                    end
                end
            end
        end
    end
    
end

