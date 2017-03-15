classdef solutionClass < handle
    %SOLUTIONCLASS Stores angular and scalar flux solution variables
    %   This is mainly to make it easier to pass information around.
    %   It will also be used to check convergence for eigenvalue
    %   problems.
    
    properties
        BCond
        keff % First index is current, second is previous iteration
        angflux
        current
        scalflux % First index is current, second is previous iteration
        fisssrc
        submesh_scalflux
    end
    
    methods
        function obj = solutionClass( mesh,xsLib,input )
            %SOLUTIONCLASS Constructor for solutionClass
            %   mesh  - The mesh
            %   xsLib - The XS library
            %   input - The iput Class container from which to initialize
            
            nsubmesh = 1;
            for i=1:mesh.nfsrcells
                if xsLib.xsSets(mesh.materials(i)).nsubxs > 0
                    nsubmesh = xsLib.xsSets(mesh.materials(i)).nsubxs;
                    break
                end
            end
            obj.keff(1:2) = 1.0;
            obj.angflux(1:2,1:xsLib.ngroups,1:input.npol,1:mesh.nfsrcells+1,nsubmesh) = 1.0;
            obj.current(1:xsLib.ngroups,1:mesh.nfsrcells+1,nsubmesh,1:2) = 0.0;
            obj.scalflux(1:xsLib.ngroups,1:mesh.nfsrcells,1:2) = 1.0;
            obj.BCond = input.BCond;
            obj.fisssrc(1:mesh.nfsrcells,1:2) = 1.0;
            if nsubmesh > 1
                obj.submesh_scalflux(1:xsLib.ngroups,1:mesh.nfsrcells,nsubmesh) = 0.0;
            end
        end
        
        function obj = updateBC( obj )
            %UPDATEBC Updates the solution to prepare for the next iteration
            %   obj - The solutionClass object to update
            
            % Set angular flux BC
            if ischar(obj.BCond(1))
                if strcmp(strtrim(obj.BCond(1,:)),'vacuum')
                    obj.angflux(1,:,:,1,:) = 0.0;
                elseif strcmp(strtrim(obj.BCond(1,:)),'reflecting')
                    obj.angflux(1,:,:,1,:) = obj.angflux(2,:,:,1,:);
                end
            end
            if ischar(obj.BCond(2))
                if strcmp(strtrim(obj.BCond(2,:)),'vacuum')
                    obj.angflux(2,:,:,end,:) = 0.0;
                elseif strcmp(strtrim(obj.BCond(2,:)),'reflecting')
                    obj.angflux(2,:,:,end,:) = obj.angflux(1,:,:,end,:);
                end
            end
            
        end
        
        function obj = updateEig( obj )
            %UPDATEEIG Calculates the new k-eff eigenvalue
            %   obj   - The solutionClass object to update
            
            obj.keff(2) = obj.keff(1);
            obj.keff(1) = obj.keff(2)*sum(obj.fisssrc(:,1))/sum(obj.fisssrc(:,2));
            
        end
        
        function [conv_flux, conv_keff] = calcResidual( obj )
            %CALCRESIDUAL Calculates the scalar flux and k-eff residuals
            %   obj   - The solutionClass object to check
            
            conv_flux = 0.0;
            for i=1:length(obj.fisssrc(:,1))
                conv_flux = max(conv_flux,abs((obj.fisssrc(i,1) - obj.fisssrc(i,2))));
            end
            conv_keff = obj.keff(1) - obj.keff(2);
            
        end
        
        function [ obj ] = calcFissSrc( obj, mesh, xsLib )
            %CALCFISSSRC Calculates the fission source in each cell
            %   obj   - The solution object to use for the FS calculation
            %   mesh  - The mesh to calculate the FS on
            %   xsLib - The XS Library to use for the calculation
            
            obj.fisssrc(:,2) = obj.fisssrc(:,1);
            obj.fisssrc(:,1) = 0.0;
            for i=1:mesh.nfsrcells
                matID = mesh.materials(i);
                obj.fisssrc(i,1) = xsLib.xsSets(matID).nufission*obj.scalflux(:,i,1);
            end
        end
    end
    
end

