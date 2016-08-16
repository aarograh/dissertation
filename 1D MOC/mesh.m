function [ matmesh, finemesh, coarsemesh ] = mesh(pinmap, cells, rad, cellmesh, hpitch)
%function [ matmesh, finemesh, coarsemesh ] = mesh(pintypes, pinmap, cells, mesh_list, widths)
%MESH Summary of this function goes here
%   pinmap   - map of pintypes used in the problem
%   cells    - list of materials in each region of each pin
%   rad      - list of radii for each region of each pin
%   cellmesh - List of FSR per region of each pin
%   hpitch    - The halfpitch of the problem

% Initialize information
npins = size(pinmap,2);
nmats = size(cells,2);
nfinecells = 0;
ncoarsecells = 0;
finemesh(nfinecells+1) = 0.0;
coarsemesh(ncoarsecells+1) = 0.0;

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
        coarsemesh(ncoarsecells+1,1) = coarsemesh(ncoarsecells) + width;
        % Set material
        oldnfinecells = nfinecells;
        nfinecells = nfinecells + cellmesh(pinmap(i),j);
        matmesh(oldnfinecells+1:nfinecells,1) = cells(pinmap(i),j);
        % Set fine mesh
        finemesh(oldnfinecells+1:nfinecells+1,1) = coarsemesh(ncoarsecells):width/cellmesh(pinmap(i),j):...
            coarsemesh(ncoarsecells+1);
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
        coarsemesh(ncoarsecells+1,1) = coarsemesh(ncoarsecells) + width;
        % Set material
        oldnfinecells = nfinecells;
        nfinecells = nfinecells + cellmesh(pinmap(i),j);
        matmesh(oldnfinecells+1:nfinecells,1) = cells(pinmap(i),j);
        % Set fine mesh
        finemesh(oldnfinecells+1:nfinecells+1,1) = coarsemesh(ncoarsecells):width/cellmesh(pinmap(i),j):...
            coarsemesh(ncoarsecells+1);
    end
end

end

