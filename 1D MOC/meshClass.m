classdef meshClass < handle
    %MESHCLASS Object to store FSR and XS meshes and materials
    %   This object store three different types of meshes.
    %   The first is matmesh.  This contains the materials
    %   for each region.  Currently these materials are
    %   stored on the FSR mesh, but this will be switched
    %   and a method added to determine the material on
    %   the XS mesh.
    %
    %   The FSR mesh stores the position of cell edges for
    %   the mesh that will be ray traced.  The XS mesh
    %   stores cell edges for unique material regions, and
    %   is a subset of the XS mesh.
    
    properties
        nxscells
        nfsrcells
        materials
        xsedges
        fsredges
        xstr
        source
    end
    
    methods
        function obj = meshClass( input )
            %MESHCLASS Generates mesh object given geometry information
            %   input - The inputClass container from which to initialize
            
            % Initialize information
            npins = size(input.pinmap,2);
            nmats = size(input.pinmats,2);
            nfinecells = 0;
            ncoarsecells = 0;
            obj.fsredges(nfinecells+1) = 0.0;
            obj.xsedges(ncoarsecells+1) = 0.0;
            hpitch = input.pitch/2.0;
            if input.diag
                hpitch = sqrt(2)*hpitch;
            end
            
            % Loop over pins
            for i=1:npins
                % Get index of last material for this cell
                for j=nmats:-1:1
                    if input.pinmats(input.pinmap(i),j) ~=0
                        nreg = j;
                        break
                    end
                end
                
                % Loop out -> in over cell descriptions
                for j=nreg:-1:1
                    ncoarsecells = ncoarsecells+1;
                    % Outermost region needs to use pin hpitch
                    if isempty(input.radii)
                        width = hpitch;
                    elseif j == nreg
                        width = (hpitch - input.radii(input.pinmap(i),j-1));
                    elseif j == 1
                        width = input.radii(input.pinmap(i),j);
                    else
                        width = (input.radii(input.pinmap(i),j) - input.radii(input.pinmap(i),j-1));
                    end
                    % Set coarse mesh
                    obj.xsedges(ncoarsecells+1,1) = obj.xsedges(ncoarsecells) + width;
                    % Set material
                    oldnfinecells = nfinecells;
                    nfinecells = nfinecells + input.pinmesh(input.pinmap(i),j);
                    obj.materials(oldnfinecells+1:nfinecells,1) = input.pinmats(input.pinmap(i),j);
                    % Set fine mesh
                    obj.fsredges(oldnfinecells+1:nfinecells+1,1) = obj.xsedges(ncoarsecells):width/...
                        input.pinmesh(input.pinmap(i),j):obj.xsedges(ncoarsecells+1);
                end
                % Loop in -> out over cell descriptions
                for j=1:nreg
                    ncoarsecells = ncoarsecells+1;
                    % Outermost region needs to use pin hpitch
                    if isempty(input.radii)
                        width = hpitch;
                    elseif j == nreg
                        width = (hpitch - input.radii(input.pinmap(i),j-1));
                    elseif j == 1
                        width = input.radii(input.pinmap(i),j);
                    else
                        width = (input.radii(input.pinmap(i),j) - input.radii(input.pinmap(i),j-1));
                    end
                    % Set coarse mesh
                    obj.xsedges(ncoarsecells+1,1) = obj.xsedges(ncoarsecells) + width;
                    % Set material
                    oldnfinecells = nfinecells;
                    nfinecells = nfinecells + input.pinmesh(input.pinmap(i),j);
                    obj.materials(oldnfinecells+1:nfinecells,1) = input.pinmats(input.pinmap(i),j);
                    % Set fine mesh
                    obj.fsredges(oldnfinecells+1:nfinecells+1,1) = obj.xsedges(ncoarsecells):width/...
                        input.pinmesh(input.pinmap(i),j):obj.xsedges(ncoarsecells+1);
                end
            end
            obj.nxscells = ncoarsecells;
            obj.nfsrcells = nfinecells;
            
        end
    end
    
end

