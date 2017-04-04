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
        npinSubTrack
        accel=false
        fipin
        bipin
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
                obj.npinSubTrack = varargin{1}.npinSubTrack;
            elseif nargin == 3
                obj.xsLib = varargin{1};
                obj.quad = varargin{2};
                obj.mesh = meshClass(varargin{3});
                obj.solution = solutionClass(obj.mesh, obj.xsLib, varargin{3});
                obj.verbose = varargin{3}.verbose;
                obj.subray = varargin{3}.subray;
                obj.npinSubTrack = varargin{3}.npinSubTrack;
            end
            if obj.subray > 0
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
                        obj.mesh.ipin(i) = -obj.mesh.ipin(i);
                    else
                        obj.mesh.materials(i,1:obj.nsubmesh) = matID;
                    end
                end
                % Set pin index flag to turn subray on and off
                if obj.nsubmesh > 1
                    pinids = zeros(max(obj.mesh.ipin),1);
                    for i=1:obj.mesh.nfsrcells
                        if obj.mesh.ipin(i) < 0
                            pinids(abs(obj.mesh.ipin(i))) = 1;
                        end
                    end
                    for i=1:obj.mesh.nfsrcells
                        if pinids(abs(obj.mesh.ipin(i)))
                            obj.mesh.ipin(i) = -abs(obj.mesh.ipin(i));
                        end
                    end
                    obj.fipin(1:obj.mesh.nfsrcells) = obj.mesh.ipin(:);
                    obj.bipin(1:obj.mesh.nfsrcells) = obj.mesh.ipin(:);
                    k = obj.mesh.nfsrcells+1;
                    fthispin = 0;
                    bthispin = max(abs(obj.bipin))+1;
                    fpinspast = obj.npinSubTrack+1;
                    bpinspast = obj.npinSubTrack+1;
                    for i=1:obj.mesh.nfsrcells
                        k = k-1;
                        flastpin = fthispin;
                        blastpin = bthispin;
                        fthispin = obj.mesh.ipin(i);
                        bthispin = obj.mesh.ipin(k);
                        if fthispin < 0
                            fpinspast = 0;
                        elseif fthispin ~= flastpin
                            fpinspast = fpinspast + 1;
                        end
                        if fpinspast <= obj.npinSubTrack
                            obj.fipin(i) = -abs(obj.fipin(i));
                        end
                        if bthispin < 0
                            bpinspast = 0;
                        elseif bthispin ~= blastpin
                            bpinspast = bpinspast + 1;
                        end
                        if bpinspast <= obj.npinSubTrack
                            obj.bipin(k) = -abs(obj.bipin(k));
                        end
                    end
                end
                % Initialize submesh flux to the scalar flux
                for j=1:obj.nsubmesh
                    obj.solution.submesh_scalflux(:,:,j) = obj.solution.scalflux(:,:,1);
                end
                obj.mesh.xstr(:,:,1:obj.nsubmesh) = 0.0;
                obj.mesh.source(:,:,1:obj.nsubmesh) = 0.0;
            % Initialize xstr since the mesh doesn't know how many groups there are
            else
                obj.nsubmesh = 1;
                obj.submesh_vol = 1.0;
                for j=1:obj.nsubmesh
                    obj.solution.submesh_scalflux(:,:,j) = obj.solution.scalflux(:,:,1);
                end
                obj.mesh.xstr(1:obj.xsLib.ngroups,1:obj.mesh.nfsrcells,1) = 0.0;
                obj.mesh.source(1:obj.xsLib.ngroups,1:obj.mesh.nfsrcells,1) = 0.0;
                obj.fipin(1:obj.mesh.nfsrcells) = 0;
                obj.bipin(1:obj.mesh.nfsrcells) = 0;
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
                if wCur && ninners == inner
                    obj.accel = true;
                else
                    obj.accel = false;
                end
                obj.solution.scalflux(:,:,2) = obj.solution.scalflux(:,:,1);
                if obj.accel
                    obj.solution.current(:,:,2) = obj.solution.current(:,:,1);
                end
                if obj.subray
                    for i=1:obj.nsubmesh
                        obj.setup( i );
                    end
                    if obj.subray == 1
                        for i=1:obj.nsubmesh
                            obj.sweep( i );
                            obj.postprocess( i );
                        end
                        obj.postprocess( 0 );
                    elseif obj.subray == 2
                        obj.sweep_subray( );
                    end
                else
                    obj.setup( );
                    obj.sweep( );
                    obj.postprocess( );
                    obj.postprocess( 0 );
                end
                if obj.accel
                    obj.solution.scalflux(:,:,1) = obj.relax*obj.solution.scalflux(:,:,1) + ...
                        (1.0-obj.relax)*obj.solution.scalflux(:,:,2);
                    obj.solution.current(:,:,:,1) = obj.relax*obj.solution.current(:,:,:,1) + ...
                        (1.0-obj.relax)*obj.solution.current(:,:,:,2);
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
            %SETUP_SUBRAY Sets up source and XS mesh for fixed source sub-ray MOC problem
            %   obj   - The FixedSourceSolverClass object to set up
            %   isubmesh - The submesh level being set up

            if ~exist('isubmesh','var')
                isubmesh = 1;
            end
            obj.mesh.source(:,:,isubmesh) = 0.0;
            for i=1:obj.mesh.nfsrcells
                matID = obj.mesh.materials(i,isubmesh);
                for j=1:obj.xsLib.ngroups
                    if obj.fipin(i) < 0 || obj.bipin(i) < 0
                        obj.mesh.source(j,i,isubmesh) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).chi(j)/obj.solution.keff(1) + ...
                            obj.xsLib.xsSets(matID).scatter(j,:)*obj.solution.submesh_scalflux(:,i,isubmesh))*0.5;
                    else
                        obj.mesh.source(j,i,isubmesh) = (obj.solution.fisssrc(i,1)*obj.xsLib.xsSets(matID).chi(j)/obj.solution.keff(1) + ...
                            obj.xsLib.xsSets(matID).scatter(j,:)*obj.solution.scalflux(:,i,2))*0.5;
                    end
                    obj.mesh.xstr(j,i,isubmesh) = obj.xsLib.xsSets(matID).transport(j);
                end
            end

        end
        
        function obj = postprocess( obj, isubmesh )
            %POSTPROCESS_SUBRAY Post-processes the sweep result for sub-ray MOC
            %   obj    - The fixedsourcesolver object to post-process
            
            if ~exist('isubmesh','var')
                isubmesh = 1;
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
                if obj.accel
                    obj.solution.current(:,:,:,1) = 0.0;
                    for j=1:obj.mesh.nfsrcells+1;
                        for k=1:obj.quad.npol
                            for g=1:obj.xsLib.ngroups
                                psibar = (obj.solution.angflux(1,g,k,j,isubmesh)- ...
                                    obj.solution.angflux(2,g,k,j,isubmesh))*0.5;
                                obj.solution.current(g,j,isubmesh,1) = obj.solution.current(g,j,isubmesh,1) + ...
                                    psibar*obj.quad.cosines(k)*obj.quad.weights(k);
                            end
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
        
        function obj = sweep_subray( obj )
            %SWEEP_SUBRAY Performs 1D MOC sweep using subray
            %   obj - The fixedsourcesolver object to sweep
            
            % Some initialization
            obj.solution.scalflux(:,:,1) = 0.0;
            obj.solution.submesh_scalflux(:) = 0.0;
            psi_in(1:obj.nsubmesh,1:2,1:obj.xsLib.ngroups,1:obj.quad.npol) = 0.0;
            source(1:obj.nsubmesh,1:2,1:obj.xsLib.ngroups) = 0.0;
            
            % Loop over all regions
            for i=1:obj.mesh.nfsrcells
                k = obj.mesh.nfsrcells-i+1;
                % If fipin is negative, we do subray in the forward direction
                if obj.fipin(i) < 0
                    lforwardSub = true;
                else
                    lforwardSub = false;
                end
                % If bipin is negative, we do subray in the backward direction
                if obj.bipin(k) < 0
                    lbackwardSub = true;
                else
                    lbackwardSub = false;
                end
                % Set boundary conditions and source
                for j=1:obj.quad.npol
                    for igroup=1:obj.xsLib.ngroups
                        for isubmesh=1:obj.nsubmesh
                            psi_in(isubmesh,1,igroup,j) = obj.solution.angflux(1,igroup,j,i,isubmesh);
                            psi_in(isubmesh,2,igroup,j) = obj.solution.angflux(2,igroup,j,k+1,isubmesh);
                        end
                    end
                end
                for igroup=1:obj.xsLib.ngroups
                    for isubmesh=1:obj.nsubmesh
                        source(isubmesh,1,igroup) = obj.mesh.source(igroup,i,isubmesh);
                        source(isubmesh,2,igroup) = obj.mesh.source(igroup,k,isubmesh);
                    end
                end
                % Modify source and boundary condition if subray shouldn't be used
                if ~lforwardSub
                    for j=1:obj.quad.npol
                        for igroup=1:obj.xsLib.ngroups
                            psi_in(:,1,igroup,j) = obj.submesh_vol*psi_in(:,1,igroup,j);
                        end
                    end
                    for igroup=1:obj.xsLib.ngroups
                        source(:,1,igroup) = obj.submesh_vol*source(:,1,igroup);
                    end
                end
                if ~lbackwardSub
                    for j=1:obj.quad.npol
                        for igroup=1:obj.xsLib.ngroups
                            psi_in(:,2,igroup,j) = obj.submesh_vol*psi_in(:,2,igroup,j);
                        end
                    end
                    for igroup=1:obj.xsLib.ngroups
                        source(:,2,igroup) = obj.submesh_vol*source(:,2,igroup);
                    end
                end
                for j=1:obj.quad.npol
                    dx1 = (obj.mesh.fsredges(i+1) - obj.mesh.fsredges(i))/obj.quad.cosines(j);
                    dx2 = (obj.mesh.fsredges(k+1) - obj.mesh.fsredges(k))/obj.quad.cosines(j);
                    for igroup=1:obj.xsLib.ngroups
                        % Forward sweep with submeshes
                        for isubmesh=1:obj.nsubmesh
                            % Forward Sweep
                            exparg = exp(-obj.mesh.xstr(igroup,i,isubmesh)*dx1);
                            obj.solution.angflux(1,igroup,j,i+1,isubmesh) = ...
                                psi_in(isubmesh,1,igroup,j)*exparg + ...
                                source(isubmesh,1,igroup)/obj.mesh.xstr(igroup,i,isubmesh)*(1.0-exparg);

                            psibar = sum(obj.solution.angflux(1,igroup,j,i:i+1,isubmesh),4)*0.5;
                            contribution = psibar*obj.quad.weights(j);
                            obj.solution.submesh_scalflux(igroup,i,isubmesh) = ...
                                obj.solution.submesh_scalflux(igroup,i,isubmesh) + contribution;
                            obj.solution.scalflux(igroup,i,1) = obj.solution.scalflux(igroup,i,1) + ...
                                contribution*obj.submesh_vol(isubmesh);
                            
                            % Backward Sweep
                            exparg = exp(-obj.mesh.xstr(igroup,k,isubmesh)*dx2);
                            obj.solution.angflux(2,igroup,j,k,isubmesh) = ...
                                psi_in(isubmesh,2,igroup,j)*exparg + ...
                                source(isubmesh,2,igroup)/obj.mesh.xstr(igroup,k,isubmesh)*(1.0 - exparg);

                            psibar = sum(obj.solution.angflux(2,igroup,j,k:k+1,isubmesh),4)*0.5;
                            contribution = psibar*obj.quad.weights(j);
                            obj.solution.submesh_scalflux(igroup,k,isubmesh) = ...
                                obj.solution.submesh_scalflux(igroup,k,isubmesh) + contribution;
                            obj.solution.scalflux(igroup,k,1) = obj.solution.scalflux(igroup,k,1) + ...
                                contribution*obj.submesh_vol(isubmesh);
                        end
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

