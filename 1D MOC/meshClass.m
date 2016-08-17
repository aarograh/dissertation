classdef meshClass
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
        materials
        xsedges
        fsredges
        xs
        source
    end
    
    methods
        function obj = meshClass(pinmap, cells, rad, cellmesh, pitch, diag)
            %MESHCLASS Generates mesh object given geometry information
            %   pinmap   - Map of pintypes used in the problem
            %   cells    - List of materials in each region of each pin
            %   rad      - List of radii for each region of each pin
            %   cellmesh - List of FSR per region of each pin
            %   pitch    - The pitch of the problem
            %   diag     - Flag to indicate if pin is being traced horizontally or diagonally
            
            % Initialize information
            npins = size(pinmap,2);
            nmats = size(cells,2);
            nfinecells = 0;
            ncoarsecells = 0;
            obj.fsredges(nfinecells+1) = 0.0;
            obj.xsedges(ncoarsecells+1) = 0.0;
            hpitch = pitch/2.0;
            if diag
                hpitch = sqrt(2)*hpitch;
            end
            
            % Loop over pins
            for i=1:npins
                % Get index of last material for this cell
                for j=nmats:-1:1
                    if cells(pinmap(i),j) ~=0
                        nreg = j;
                        break
                    end
                end
                
                % Loop out -> in over cell descriptions
                for j=nreg:-1:1
                    ncoarsecells = ncoarsecells+1;
                    % Outermost region needs to use pin hpitch
                    if j == nreg
                        width = (hpitch - rad(pinmap(i),j-1));
                    elseif j == 1
                        width = rad(pinmap(i),j);
                    else
                        width = (rad(pinmap(i),j) - rad(pinmap(i),j-1));
                    end
                    % Set coarse mesh
                    obj.xsedges(ncoarsecells+1,1) = obj.xsedges(ncoarsecells) + width;
                    % Set material
                    oldnfinecells = nfinecells;
                    nfinecells = nfinecells + cellmesh(pinmap(i),j);
                    obj.materials(oldnfinecells+1:nfinecells,1) = cells(pinmap(i),j);
                    % Set fine mesh
                    obj.fsredges(oldnfinecells+1:nfinecells+1,1) = ...
                        obj.xsedges(ncoarsecells):width/cellmesh(pinmap(i),j):obj.xsedges(ncoarsecells+1);
                end
                % Loop in -> out over cell descriptions
                for j=1:nreg
                    ncoarsecells = ncoarsecells+1;
                    % Outermost region needs to use pin hpitch
                    if j == nreg
                        width = (hpitch - rad(pinmap(i),j-1));
                    elseif j == 1
                        width = rad(pinmap(i),j);
                    else
                        width = (rad(pinmap(i),j) - rad(pinmap(i),j-1));
                    end
                    % Set coarse mesh
                    obj.xsedges(ncoarsecells+1,1) = obj.xsedges(ncoarsecells) + width;
                    % Set material
                    oldnfinecells = nfinecells;
                    nfinecells = nfinecells + cellmesh(pinmap(i),j);
                    obj.materials(oldnfinecells+1:nfinecells,1) = cells(pinmap(i),j);
                    % Set fine mesh
                    obj.fsredges(oldnfinecells+1:nfinecells+1,1) = ...
                        obj.xsedges(ncoarsecells):width/cellmesh(pinmap(i),j):obj.xsedges(ncoarsecells+1);
                end
            end
            
        end
    end
    
end

