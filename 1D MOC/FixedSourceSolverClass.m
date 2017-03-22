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
        nsubmesh=0
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
                % Get number of submeshes.  Assume the same number or 0 everywhere
                for i=1:obj.mesh.nfsrcells
                    if obj.xsLib.xsSets(obj.mesh.materials(i)).nsubxs > 0
                        break
                    end
                end
                % Set submesh volume fractions
                obj.nsubmesh = length(obj.xsLib.xsSets(obj.mesh.materials(i)).subfracs);
                obj.submesh_vol = obj.xsLib.xsSets(obj.mesh.materials(i)).subfracs;
                % Split mesh into multiple meshes
                for i=1:obj.mesh.nfsrcells
                    matID = obj.mesh.materials(i);
                    nsubmesh = obj.xsLib.xsSets(matID).nsubxs;
                    if nsubmesh > 0
                        for j=1:nsubmesh
                            obj.mesh.materials(i,j) = obj.xsLib.xsSets(matID).subxs(j).ID;
                        end
                    else
                        obj.mesh.materials(i,1:obj.nsubmesh) = matID;
                    end
                end
                % Initialize submesh flux to the scalar flux
                for j=1:obj.nsubmesh
                    obj.solution.submesh_scalflux(:,:,:,j) = obj.solution.scalflux(:,:,:);
                end
                obj.mesh.xstr(:,:,1:obj.nsubmesh) = 0.0;
                obj.mesh.source(:,:,1:obj.nsubmesh) = 0.0;
                if obj.homog == 1
                    i = length(obj.solution.angflux(1,1,1,1,:))+1;
                    obj.solution.angflux(:,:,:,:,i) = 0.0;
                    obj.solution.current(:,:,i,:) = 0.0;
                end
            % Initialize xstr since the mesh doesn't know how many groups there are
            else
                obj.mesh.xstr(1:obj.xsLib.ngroups,1:obj.mesh.nfsrcells,1) = 0.0;
                obj.mesh.source(1:obj.xsLib.ngroups,1:obj.mesh.nfsrcells,1) = 0.0;
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
                convcrit=1.0e-4;
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
                        for i=1:obj.nsubmesh
                            obj.setup_subray( i );
                            obj.sweep( i );
                            obj.postprocess_subray( i );
                        end
                        if obj.homog == 1
                            obj.setup_subrayHomog( );
                            obj.sweep( obj.nsubmesh+1 );
                            obj.postprocess_unified( obj.nsubmesh+1 );
                        else
                            obj.postprocess_subray( );
                        end
                    else
                        obj.setup( );
                        obj.sweep( );
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
        
        function obj = setup( obj, isubmesh )
            %SETUP Sets up source and XS mesh for fixed source MOC problem
            %   obj    - The FixedSourceSolverClass object to set up
            %   isubmesh  - The submesh level being set up
            
            obj.mesh.source(:) = 0.0;
            for i=1:obj.mesh.nfsrcells
                matID = obj.mesh.materials(i,1);
                for j=1:obj.xsLib.ngroups
                    % Use old scalar flux to do Jacobi style iteration
                    obj.mesh.source(j,i,1) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).chi(j)/obj.solution.keff(1) + ...
                        obj.xsLib.xsSets(matID).scatter(j,:)*obj.solution.scalflux(:,i,2))*0.5;
                    obj.mesh.xstr(j,i,1) = obj.xsLib.xsSets(matID).transport(j);
                end
            end
            
        end

        function obj = setup_subray( obj, isubmesh )
            %SETUP_SUBRAY Sets up source and XS mesh for fixed source sub-ray MOC problem
            %   obj   - The FixedSourceSolverClass object to set up
            %   isubmesh - The submesh level being set up

            obj.mesh.source(:,:,isubmesh) = 0.0;
            for i=1:obj.mesh.nfsrcells
                matID = obj.mesh.materials(i,isubmesh);
                for j=1:obj.xsLib.ngroups
                    if obj.xsLib.xsSets(matID).nsubxs > 0
                        obj.mesh.source(j,i,isubmesh) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).subxs(isubmesh).chi(j)/obj.solution.keff(1) + ...
                            obj.xsLib.xsSets(matID).subxs(isubmesh).scatter(j,:)*obj.solution.submesh_scalflux(:,i,isubmesh))*0.5;
                        % Commented out bit uses source from homogenized scalar flux, gives terrible results
                        % obj.mesh.source(j,i,isubmesh) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).subxs(isubmesh).chi(j)/obj.solution.keff(1) + ...
                        %     obj.xsLib.xsSets(matID).subxs(isubmesh).scatter(j,:)*obj.solution.scalflux(:,i,2))*0.5;
                        obj.mesh.xstr(j,i,isubmesh) = obj.xsLib.xsSets(matID).subxs(isubmesh).transport(j);
                    else
                        obj.mesh.source(j,i,isubmesh) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).chi(j)/obj.solution.keff(1) + ...
                            obj.xsLib.xsSets(matID).scatter(j,:)*obj.solution.submesh_scalflux(:,i,isubmesh))*0.5;
                        % Commented out bit uses source from homogenized scalar flux, gives terrible results
                        % obj.mesh.source(j,i,isubmesh) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).chi(j)/obj.solution.keff(1) + ...
                        %     obj.xsLib.xsSets(matID).scatter(j,:)*obj.solution.scalflux(:,i,2))*0.5;
                        obj.mesh.xstr(j,i,isubmesh) = obj.xsLib.xsSets(matID).transport(j);
                    end
                end
            end

        end
        
        function obj = setup_subrayHomog( obj )
            %SETUP_SUBRAYHOMOG Sets up source and XS mesh for final fixed source MOC problem after subray
            %sweeps
            %   obj - The FixedSourceSolverClass object to set up
            
            obj.mesh.source(:) = 0.0;
            
            for k=1:obj.nsubmesh
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
        
        function obj = postprocess_unified( obj, isubmesh )
            %POSTPROCESS_UNIFIED Post-rpocess the sweep result for regular MOC
            %   obj    - The fixedsourcesolver object to post-process
            %   isubmesh  - The submesh level to use
            
            if ~exist('isubmesh','var')
                isubmesh = 1;
            end
            obj.solution.scalflux(:,:,1) = 0.0;
            for j=1:obj.mesh.nfsrcells
                for k=1:obj.quad.npol
                    for g=1:obj.xsLib.ngroups
                        psibar = sum(obj.solution.angflux(1,g,k,j:j+1,isubmesh));
                        psibar = 0.5*(psibar + sum(obj.solution.angflux(2,g,k,j:j+1,isubmesh)));
                        obj.solution.scalflux(g,j,1) = obj.solution.scalflux(g,j,1) + ...
                            psibar*obj.quad.weights(k);
                    end
                end
            end
            
        end
        
        function obj = postprocess_subray( obj, isubmesh )
            %POSTPROCESS_SUBRAY Post-processes the sweep result for sub-ray MOC
            %   obj    - The fixedsourcesolver object to post-process
            
            if ~exist('isubmesh','var')
                isubmesh = 0;
            end

            if isubmesh == 0
                obj.solution.scalflux(:,:,1) = 0.0;
                for i=1:obj.nsubmesh
                    obj.solution.scalflux(:,:,1) = obj.solution.scalflux(:,:,1) + ...
                    obj.submesh_vol(i)*obj.solution.submesh_scalflux(:,:,i);
                end
            else
                obj.solution.submesh_scalflux(:,:,isubmesh) = 0.0;
                for j=1:obj.mesh.nfsrcells
                    for k=1:obj.quad.npol
                        for g=1:obj.xsLib.ngroups
                            psibar = sum(sum(obj.solution.angflux(1:2,g,k,j:j+1,isubmesh),4),1)*0.5;
                            obj.solution.submesh_scalflux(g,j,isubmesh) = obj.solution.submesh_scalflux(g,j,isubmesh) + ...
                                psibar*obj.quad.weights(k);
                        end
                    end
                end
            end
            
        end
        
        function obj = sweep( obj, isubmesh )
            %SWEEP Performs 1D MOC sweep for a single ray with multiple polars
            %   obj      - The fixedsourcesolver object to sweep
            %   isubmesh - The submesh level to sweep
            
            if ~exist('isubmesh','var')
                isubmesh = 0;
            end
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
                        exparg = exp(-obj.mesh.xstr(igroup,i,isubmesh)*dx1);
                        obj.solution.angflux(1,igroup,j,i+1,isubmesh) = ...
                            obj.solution.angflux(1,igroup,j,i,isubmesh)*exparg + ...
                            obj.mesh.source(igroup,i,isubmesh)/obj.mesh.xstr(igroup,i,isubmesh)*(1 - exparg);

                        % Backward Sweep
                        exparg = exp(-obj.mesh.xstr(igroup,k,isubmesh)*dx2);
                        obj.solution.angflux(2,igroup,j,k,isubmesh) = ...
                            obj.solution.angflux(2,igroup,j,k+1,isubmesh)*exparg + ...
                            obj.mesh.source(igroup,k,isubmesh)/obj.mesh.xstr(igroup,k,isubmesh)*(1 - exparg);
                    end
                end
            end
        end
        
        function obj = sweep_wCur( obj )
            %SWEEP_wCur Performs 1D MOC sweep for a single ray with multiple polars while
            %           tallying currents for cmfd acceleration
            %   obj    - The eigensolver object to sweep
            
            %TODO: update to handle sub-ray
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
                % [oldsource,newsource]
                maxdiff = max(max(abs((oldsource - newsource)./oldsource)),maxdiff);
            end
            
        end
    end
    
end

