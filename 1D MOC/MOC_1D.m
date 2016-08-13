close all, clear all, clc;

%% User Options
npins = 3;
pinmap = ['fuel   '; 'control'; 'fuel   '];
npol = 1;
nmesh_mod = 2;
nmesh_fuel = 5;
nmesh_control = 9;
nmesh_gap = 1;
nmesh_clad = 1;

%% Setup Problem
% Materials
id_mod = 1;
id_fuel = 2;
id_control = 3;
id_clad = 4;
id_gap = 5;
meshes = [nmesh_mod, nmesh_fuel, nmesh_control, nmesh_clad, nmesh_gap];
xstr = [1.0; 1.0; 1.0; 1.0; 1.0];

% Geometry
fuelrad = 0.4096;
crrad = 0.4;
pitch = 1.26;
halfpitch = pitch/2.0;

% Ray Parameters
if (npol == 1)
    quad = 1.0;
    mu = cos(pi/4);
end

% Mesh
n=0;
for i=1:npins
    if strcmp(pinmap(i,:),'fuel   ')
        %set up fuel pin
        fprintf('Setting up fuel pin\n')
        regions=[id_mod,id_clad,id_gap,id_fuel,id_gap,id_clad,id_mod];
    elseif strcmp(pinmap(i,:),'control')
        %set up control rod
        fprintf('Setting up control pin\n')
        regions = [id_mod,id_clad,id_gap,id_control,id_gap,id_clad,id_mod];
    end
    for j=1:size(regions)
        for k=1:meshes(regions(j))
            n = n+1;
            imats(n) = regions(j);
            dx(n) = widths(imats(n));
        end
    end
end

%% Solve Problem

%% Generate Plots