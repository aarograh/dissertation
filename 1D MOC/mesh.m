function [ matmesh, finemesh, coarsemesh ] = mesh(pintypes, pinmap, cells, mesh_list, widths)
%MESH Summary of this function goes here
%   pintypes - list of pintypes that could be used in the problem
%   pinmap   - map of pintypes used in the problem
%   cells    - list of materials in each region of each pin

% Get sizes
npins = size(pinmap,1);

ncoarsecells = 0;
nfinecells = 0;
coarsemesh(ncoarsecells+1) = 0.0;
finemesh(nfinecells+1) = 0.0;
nfinecells = nfinecells+1;
ncoarsecells = ncoarsecells+1;
for i=1:npins
    for j=1:size(pintypes,1)
        if strcmp(pinmap(i,:),pintypes(j,:))
            for k=size(cells,2):-1:1
                if (cells(j,k) ~= 0)
                    break
                end
            end
            break
        end
    end
    regions = cells(j,1:k)';
    for j=1:size(regions) 
        for k=1:mesh_list(regions(j))
            matmesh(nfinecells) = regions(j);
            finemesh(nfinecells+1,1) = finemesh(nfinecells) + widths(matmesh(nfinecells));
            nfinecells = nfinecells+1;
        end
         coarsemesh(ncoarsecells+1,1) = coarsemesh(ncoarsecells) + widths(regions(j))*mesh_list(regions(j));
         ncoarsecells = ncoarsecells+1;
    end
end

end

