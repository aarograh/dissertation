classdef solutionClass
    %SOLUTIONCLASS Stores angular and scalar flux solution variables
    %   This is mainly to make it easier to pass information around.
    %   It will also be used to check convergence for eigenvalue
    %   problems.
    
    properties
        BCond
        keff % First index is current, second is previous iteration
        angflux
        angflux_cellavg
        scalflux % First index is current, second is previous iteration
        fisssrc
        fluxnorm
    end
    
    methods
        function obj = solutionClass( ncells,npol,ngroups,BCond )
            %SOLUTIONCLASS Constructor for solutionClass
            %   ncells  - The number of cells in the problem
            %   npol    - The number of polar angles in the problem
            %   ngroups - The number of energy groups in the problem
            %   BCond   - The boundary condition for the angular flux
            
            obj.keff(1:2) = 1.0;
            obj.angflux(1:ncells+1,1:npol,1:2,1:ngroups) = 0.0;
            obj.scalflux(1:ncells,1:ngroups,1:2) = 1.0;
            obj.BCond = BCond;
            obj.fisssrc(1:ncells) = 0.0;
            obj.fluxnorm = 0.0;
        end
        
        function obj = update( obj )
            %UPDATEBC Updates the solution to prepare for the next iteration
            %   Resets the rest of the angular flux array to 0.0
            %   Applies the appropriate boundary conditions to the angular flux
            %   Saves the scalar flux before zeroing it
            %   TODO: update eigenvalue
            
            obj.angflux(2:end,:,1,:) = 0.0;
            obj.angflux(1:end-1,:,2,:) = 0.0;
            obj.angflux_cellavg = 0.0;
            obj.scalflux(:,:,2) = obj.scalflux(:,:,1);
            obj.scalflux(:,:,1) = 0.0;
            obj.keff(2) = obj.keff(1);
            
            % Set angular flux BC
            if ischar(obj.BCond(1))
                if strcmp(strtrim(obj.BCond(1,:)),'vacuum')
                    obj.angflux(1,:,1,:) = 0.0;
                end
            end
            if ischar(obj.BCond(2))
                if strcmp(strtrim(obj.BCond(2,:)),'vacuum')
                    obj.angflux(end,:,2,:) = 0.0;
                end
            end
            
            if abs(obj.fluxnorm < eps)
                obj.fluxnorm = sum(obj.fisssrc);
            end
            
        end
        
        function obj = updateEig( obj, quad, mesh, xsLib )
            %UPDATEEIG Calculates the new k-eff eigenvalue
            %   obj   - The solutionClass object to update
            %   quad  - The quadrature to integrate the solution
            %   mesh  - The mesh that the solution is on
            %   xsLib - The XS Library used by the problem
            
            numerator = 0.0;
            denominator = 0.0;
            % Accumulate fission source and absorption rate
            for j=1:xsLib.ngroups
                if j == xsLib.ngroups
                    dE = xsLib.groupBounds(xsLib.ngroups);
                else
                    dE = xsLib.groupBounds(j) - xsLib.groupBounds(j+1);
                end
                for i=1:quad.npol
                    denominator = denominator + obj.angflux(1,i,2,j)*quad.weights(i)*dE;
                    denominator = denominator + obj.angflux(end,i,1,j)*quad.weights(i)*dE;
                end
                for i=1:mesh.nfsrcells
                    dxdE = dE*(mesh.fsredges(i+1) - mesh.fsredges(i));
                    flux = obj.scalflux(i,j,1);
                    matID = mesh.materials(i);
                    numerator = numerator + flux*xsLib.xsSets(matID).nufission(j)*dxdE;
                    denominator = denominator + flux*xsLib.xsSets(matID).absorption(j)*dxdE;
                end
            end
            
            obj.keff(1) = numerator/denominator;
        end
        
        function [conv_flux, conv_keff] = calcResidual( obj, mesh, xsLib)
            %CALCRESIDUAL Calculates the scalar flux and k-eff residuals
            %   obj   - The solutionClass object to check
            %   mesh  - The mesh that the solution is on
            %   xsLib - The XS Library used by the calculation
            
            conv_flux = 0.0;
            for j=1:xsLib.ngroups
                for i=1:mesh.nfsrcells
                    conv_flux = max(conv_flux,abs(obj.scalflux(i,j,1) - obj.scalflux(i,j,2)));
%                     display(sprintf('Engery group %i: %g %g',j,obj.scalflux(i,j,1),obj.scalflux(i,j,2)));
                end
            end
            conv_keff = obj.keff(1) - obj.keff(2);
            
        end
        
        function [ obj ] = calcFissSrc( obj, mesh, xsLib )
            %CALCFISSSRC Calculates the fission source in each cell
            %   obj   - The solution object to use for the FS calculation
            %   mesh  - The mesh to calculate the FS on
            %   xsLib - The XS Library to use for the calculation
            
            obj.fisssrc = 0.0;
            for i=1:mesh.nfsrcells
                matID = mesh.materials(i);
                obj.fisssrc(i) = sum(obj.scalflux(i,:,1).*xsLib.xsSets(matID).nufission)/...
                    obj.keff(1);
            end
            obj.fisssrc
        end
        
        function [ obj ] = normalize( obj )
            %NORMALIZE Normalizes the flux to the original fission source,
            %          or some specified constant
            %   obj - The solutionClass object to normalize
            
            obj.scalflux(:,:,1) = obj.scalflux(:,:,1)*obj.fluxnorm/sum(obj.fisssrc);
            
        end
    end
    
end

