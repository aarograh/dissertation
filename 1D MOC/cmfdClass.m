classdef cmfdClass < handle
    %CMFDCLASS CMFD Accelerator class
    %   This class contains the data and methods required to perform
    %   1D CMFD acceleration
    
    properties
        flux
        fisssrc
        xstr
        xst
        xsrm
        xssc
        xsnf
        xsch
        dtils
        dhats
        A
        x
        b
        keff
        xsLib
        mesh
        ngroups
        ncells
        regPerCell
        cellwidths
        firstIteration = true
    end
    
    methods
        function obj = cmfdClass( input, mesh, xsLib)
            %CMFDCLASS Constructor for cmfdClass object
            %   input - The inputClass object from which to initialize
            %   mesh  - The meshClass object to solve
            %   xsLib - The xsLibraryClass object with XS data
            
            obj.ncells = length(input.pinmap);
            obj.ngroups = xsLib.ngroups;
            obj.xsLib = xsLib;
            obj.mesh = mesh;
            for ipin=1:obj.ncells
                % Multiply by 2 since only half the pin is being described
                obj.regPerCell(ipin) = 2*sum(input.pinmesh(input.pinmap(ipin),:),2);
            end
%             for ipin=1:npins
%                 for imat=nmats:-1:1
%                     if input.pinmats(input.pinmap(ipin),imat)
%                         break
%                     end
%                 end
%                 obj.regPerCell(obj.ncells+1:obj.ncells+imat) = input.pinmesh(input.pinmap(ipin),imat:-1:1);
%                 obj.ncells = obj.ncells + imat;
%                 obj.regPerCell(obj.ncells+1:obj.ncells+imat) = input.pinmesh(input.pinmap(ipin),1:imat);
%                 obj.ncells = obj.ncells + imat;
%             end
            
            obj.flux(1:obj.ncells,1:obj.ngroups,1:2) = 0.0;
            obj.fisssrc(1:obj.ncells,1:2) = 0.0;
            obj.xstr(1:obj.ncells,1:obj.ngroups) = 0.0;
            obj.xst(1:obj.ncells,1:obj.ngroups) = 0.0;
            obj.xsrm(1:obj.ncells,1:obj.ngroups) = 0.0;
            obj.xssc(1:obj.ncells,1:obj.ngroups,1:obj.ngroups) = 0.0;
            obj.xsnf(1:obj.ncells,1:obj.ngroups) = 0.0;
            obj.xsch(1:obj.ncells,1:obj.ngroups) = 0.0;
            obj.dtils(1:obj.ncells+1,1:obj.ngroups) = 0.0;
            obj.dhats(1:obj.ncells+1,1:obj.ngroups) = 0.0;
            obj.A(1:obj.ncells*obj.ngroups,1:obj.ncells*obj.ngroups) = 0.0;
            obj.x(1:obj.ncells*obj.ngroups,1) = 0.0;
            obj.b(1:obj.ncells*obj.ngroups,1) = 0.0;
            
        end
        
        function solution = solve( obj, solution, mesh )
            %SOLVE Performs power iteration to solve the CMFD system
            %   obj      - The cmfdClass object to set up
            %   solution - The solutionClass object to use for the setup
            %   mesh     - The meshClass object on which to solve
            
            obj.keff(1:2) = solution.keff(1);
            converged = 0;
            
            obj.homogenize(solution, mesh);
            obj.setupMatrix();
            iters = 0;
            while ~converged
                iters = iters + 1;
                obj.setupSource();
                obj.step();
                [conv_flux, conv_keff] = obj.calcResidual( );
                if conv_flux < 1.0e-8 && conv_keff < 1.0e-8
                    converged = 1;
                elseif iters == 20
                    display(sprintf('Too many iterations.  Giving up.'));
                    break
                else
                    display(sprintf('  CMFD keff, residuals (k-eff, flux): %0.8f, %0.8f, %0.8f',...
                        obj.keff(1),conv_keff,conv_flux));
                end
            end
            obj.project(solution, mesh);
            
            solution.keff(2) = solution.keff(1);
            solution.keff(1) = obj.keff(1);
            
        end
        
        function obj = step( obj )
            %STEP Performs a single CMFD iteration
            %   obj      - The cmfdClass object to set up
            
            % Solve linear system
            obj.x = obj.A\obj.b;
            
            % Update flux and fission source
            obj.fisssrc(:,2) = obj.fisssrc(:,1);
            obj.fisssrc(:,1) = 0.0;
            irow = 0;
            for i=1:obj.ncells
                for g=1:obj.ngroups
                    irow = irow + 1;
                    obj.flux(i,g,1) = obj.x(irow);
                    obj.fisssrc(i,1) = obj.fisssrc(i,1) + ...
                        obj.flux(i,g,1)*obj.xsnf(i,g,1)*obj.cellwidths(i);
                end
            end
            
            % Update eigenvalue
            obj.keff(2) = obj.keff(1);
            obj.keff(1) = obj.keff(2)*sum(obj.fisssrc(:,1))/sum(obj.fisssrc(:,2));
            
        end
        
        function [conv_flux, conv_keff] = calcResidual(obj)
            %CALCRESIDUAL Calculates the scalar flux and k-eff residuals
            %   obj - cmfdClass object whose convergence needs to be calculated
            
            conv_flux = 0.0;
            for i=1:obj.ncells
                conv_flux = max(conv_flux,abs((obj.fisssrc(i,1) - obj.fisssrc(i,2))/...
                    obj.fisssrc(i,2)));
            end
            conv_keff = obj.keff(1) - obj.keff(2);
            
        end
        
        function obj = homogenize( obj, solution, mesh )
            %HOMOGENIZE Homogenizes the cmfd mesh
            %   obj      - The cmfdClass object to set up
            %   solution - The solutionClass object to use for the setup
            %   mesh     - The meshClass object on which to solve
            
            obj.xstr(:) = 0.0;
            obj.xst(:) = 0.0;
            obj.xsrm(:) = 0.0;
            obj.xssc(:) = 0.0;
            obj.xsnf(:) = 0.0;
            obj.xsch(:) = 0.0;
            display(solution.scalflux)
            
            for g=1:obj.ngroups
                icell = 0;
                for i=1:obj.ncells
                    volsum = 0.0;
                    flxvolsum = 0.0;
                    fisssrcsum = 0.0;
                    for j=1:obj.regPerCell(i)
                        icell = icell + 1;
                        imat = mesh.materials(icell);
                        dx = (mesh.fsredges(icell+1)-mesh.fsredges(icell));
                        volsum = volsum + dx;
                        flxvol = solution.scalflux(icell,g,1)*dx;
                        flxvolpsi = solution.fisssrc(icell,1)*dx;
                        flxvolsum = flxvolsum + flxvol;
                        fisssrcsum = fisssrcsum + flxvolpsi;
                        obj.xstr(i,g) = obj.xstr(i,g) + obj.xsLib.xsSets(imat).transport(g)*flxvol;
                        obj.xst(i,g) = obj.xst(i,g) + obj.xsLib.xsSets(imat).transport(g)*flxvol;
                        for g2=1:obj.ngroups
                            obj.xssc(i,g2,g) = obj.xssc(i,g2,g) + ...
                                obj.xsLib.xsSets(imat).scatter(g2,g)*flxvol;
                        end
                        obj.xsnf(i,g) = obj.xsnf(i,g) + obj.xsLib.xsSets(imat).nufission(g)*flxvol;
                        obj.xsch(i,g) = obj.xsch(i,g) + ...
                            flxvolpsi*obj.xsLib.xsSets(imat).chi(g);
                    end
                    
                    % Calculate Homogenized cross-sections and flux
                    obj.flux(i,g,1:2) = flxvolsum/volsum;
                    obj.xstr(i,g) = obj.xstr(i,g)/flxvolsum;
                    obj.xst(i,g) = obj.xst(i,g)/flxvolsum;
                    obj.xsnf(i,g) = obj.xsnf(i,g)/flxvolsum;
                    obj.xsch(i,g) = obj.xsch(i,g)/fisssrcsum;
                    if g == 1
                        obj.fisssrc(i,1:2) = obj.flux(i,g,1)*obj.xsnf(i,g);
                    end
                    for g2=1:obj.ngroups
                        obj.xssc(i,g2,g) = obj.xssc(i,g2,g)/flxvolsum;
                    end
                    obj.xsrm(i,g) = obj.xst(i,g) - obj.xssc(i,g,g);
                    
                    % Calculate coupling coefficients
                    %   Interior surface
                    if i > 1
                        obj.dtils(i-1,g) = 2.0/(3.0*(volsum*obj.xstr(i,g) + oldvolsum*obj.xstr(i-1,g)));
                        obj.dhats(i-1,g) = (solution.current(icell-obj.regPerCell(i),g,1) + ...
                            obj.dtils(i-1,g)*(obj.flux(i,g,1) - obj.flux(i-1,g,1)))/ ...
                            (obj.flux(i,g,1) + obj.flux(i-1,g,1));
                    %   Left boundary
                    else
                        if strcmp(strtrim(solution.BCond(1,:)),'vacuum')
                            alpha = 0.5;
                        else
                            alpha = 0.0;
                        end
                        obj.dtils(1,g) = alpha/(1.0 + 3.0*(volsum/2.0)*alpha*obj.xstr(1,g));
                        if strcmp(strtrim(solution.BCond(1,:)),'reflecting')
                            obj.dhats(1,g) = 0.0;
                        else
                            obj.dhats(1,g) = (solution.current(1,g,1) + obj.dtils(1,g)*obj.flux(1,g,1))/ ...
                                obj.flux(1,g,1);
                        end
                    end
                    %   Right boundary
                    if i == obj.ncells
                        if strcmp(strtrim(solution.BCond(2,:)),'vacuum')
                            alpha = 0.5;
                        else
                            alpha = 0.0;
                        end
                        obj.dtils(end,g) = alpha/(1.0 + 3.0*(volsum/2.0)*alpha*obj.xstr(end,g));
                        if strcmp(strtrim(solution.BCond(2,:)),'reflecting')
                            obj.dhats(end,g) = 0.0;
                        else
                            obj.dhats(end,g) = (solution.current(end,g,1) - obj.dtils(end,g)*obj.flux(end,g,1))/ ...
                                obj.flux(end,g,1);
                        end
                    end
                    obj.cellwidths(i) = volsum;
                    oldvolsum = volsum;
                end
            end
            display(obj.flux)
        end
        
        function obj = setupMatrix( obj )
            %SETUPMATRIX Sets up the matrix for the CMFD linear system
            %   obj - cmfdClass to set up
            
            irow = 0;
            obj.A(:) = 0.0;
            obj.x(:) = 0.0;
            for i=1:obj.ncells
                for g=1:obj.ngroups
                    irow = irow + 1;
                    % Sum of coupling coefficients goes on diagonal
                    obj.A(irow,irow) = obj.dtils(i+1,g) + obj.dhats(i+1,g) + ...
                        obj.dtils(i,g) - obj.dhats(i,g);
                    % Add coupling coefficients for east neighbor
                    if i < obj.ncells
                        obj.A(irow,irow+obj.ngroups) = -obj.dtils(i+1,g) + obj.dhats(i+1,g);
                    end
                    % Add coupling coefficients for west neighbor
                    if i > 1
                        obj.A(irow,irow-obj.ngroups) = -obj.dtils(i,g) - obj.dhats(i,g);
                    end
                    
                    % Add total reaction rate
                    obj.A(irow,irow) = obj.A(irow,irow) + obj.xsrm(i,g)*obj.cellwidths(i);
                    
                    % Add scattering terms
                    for g2=1:obj.ngroups
                        if g ~= g2
                            obj.A(irow,irow-g+g2) = obj.A(irow,irow-g+g2) + ...
                                obj.xssc(i,g,g2)*obj.cellwidths(i);
                        end
                    end
                end
            end
            
        end
        
        function obj = setupSource( obj )
            %SETUPSOURCE Sets up the source for the CMFD linear system
            %   obj - cmfdClass to set up
            
            irow = 0;
            obj.b(:) = 0.0;
            for i=1:obj.ncells
                for g=1:obj.ngroups
                    irow = irow + 1;
                    % Add source
                    obj.b(irow,1) = obj.fisssrc(i,1)/obj.keff(1);
                end
            end
        end
        
        function obj = project( obj, solution, mesh )
            %PROJECT Projects the CMFD solution onto the MOC mesh
            %   obj      - The cmfdClass object to set up
            %   solution - The solutionClass object to use for the setup
            %   mesh     - The meshClass object on which to solve
            
            icell=0;
            for i=1:obj.ncells
                regcells = obj.regPerCell(i);
                for g=1:obj.ngroups
                    scale = obj.flux(i,g,1)/obj.flux(i,g,2);
                    solution.scalflux(icell+1:icell+regcells,g,1) = ...
                        solution.scalflux(icell+1:icell+regcells,g,1)*scale;
                end
                icell = icell + regcells;
            end
            
        end
    end
    
end

