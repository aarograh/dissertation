classdef solutionClass
    %SOLUTIONCLASS Stores angular and scalar flux solution variables
    %   This is mainly to make it easier to pass information around.
    %   It will also be used to check convergence for eigenvalue
    %   problems.
    
    properties
        BCond
        keff % First index is current, second is previous iteration
        angflux
        scalflux % First index is current, second is previous iteration
        fisssrc
        fluxnorm
    end
    
    methods
        function obj = solutionClass( ncells,ngroups,input )
            %SOLUTIONCLASS Constructor for solutionClass
            %   ncells  - The number of cells in the problem
            %   ngroups - The number of energy groups in the problem
            %   input   - The iput Class container from which to initialize
            
            obj.keff(1:2) = 1.0;
            obj.angflux(1:ncells+1,1:input.npol,1:2,1:ngroups) = 0.0;
            obj.scalflux(1:ncells,1:ngroups,1:2) = 1.0;
            obj.BCond = input.BCond;
            obj.fisssrc(1:ncells,1:2) = 0.0;
            obj.fluxnorm = 0.0;
        end
        
        function obj = update( obj, iouter )
            %UPDATEBC Updates the solution to prepare for the next iteration
            %   Resets the rest of the angular flux array to 0.0
            %   Applies the appropriate boundary conditions to the angular flux
            %   Saves the scalar flux before zeroing it
            %   TODO: update eigenvalue
            
            obj.scalflux(:,:,2) = obj.scalflux(:,:,1);
            obj.scalflux(:,:,1) = 0.0;
            obj.keff(2) = obj.keff(1);
            obj.fisssrc(:,2) = obj.fisssrc(:,1);
            
            % Set angular flux BC
            if ischar(obj.BCond(1))
                if strcmp(strtrim(obj.BCond(1,:)),'vacuum')
                    obj.angflux(1,:,1,:) = 0.0;
                elseif strcmp(strtrim(obj.BCond(1,:)),'reflecting')
                    if iouter == 1
                        obj.angflux(1,:,1,:) = 1.0;
                    else
                        obj.angflux(1,:,1,:) = obj.angflux(1,:,2,:);
                    end
                end
            end
            if ischar(obj.BCond(2))
                if strcmp(strtrim(obj.BCond(2,:)),'vacuum')
                    obj.angflux(end,:,2,:) = 0.0;
                elseif strcmp(strtrim(obj.BCond(2,:)),'reflecting')
                    if iouter == 1
                        obj.angflux(end,:,2,:) = 1.0;
                    else
                        obj.angflux(end,:,2,:) = obj.angflux(end,:,1,:);
                    end
                end
            end
            
        end
        
        function obj = updateEig( obj )
            %UPDATEEIG Calculates the new k-eff eigenvalue
            %   obj   - The solutionClass object to update
            
            obj.keff(1) = obj.keff(2)*sum(obj.fisssrc(:,1))/sum(obj.fisssrc(:,2));
            
        end
        
        function [conv_flux, conv_keff] = calcResidual( obj, mesh, xsLib)
            %CALCRESIDUAL Calculates the scalar flux and k-eff residuals
            %   obj   - The solutionClass object to check
            %   mesh  - The mesh that the solution is on
            %   xsLib - The XS Library used by the calculation
            
            conv_flux = 0.0;
            for j=1:xsLib.ngroups
                for i=1:mesh.nfsrcells
                    conv_flux = max(conv_flux,abs((obj.scalflux(i,j,1) - obj.scalflux(i,j,2))/...
                        obj.scalflux(i,j,2)));
                end
            end
            conv_keff = obj.keff(1) - obj.keff(2);
            
        end
        
        function [ obj ] = calcFissSrc( obj, mesh, xsLib )
            %CALCFISSSRC Calculates the fission source in each cell
            %   obj   - The solution object to use for the FS calculation
            %   mesh  - The mesh to calculate the FS on
            %   xsLib - The XS Library to use for the calculation
            
            obj.fisssrc(:,1) = 0.0;
            for i=1:mesh.nfsrcells
                matID = mesh.materials(i);
                obj.fisssrc(i,1) = sum(obj.scalflux(i,:,1).*xsLib.xsSets(matID).nufission);
            end
        end
    end
    
end

