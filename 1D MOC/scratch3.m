close all
clear all
clc

% OnePin7group_mix_fuel_rodded
% ThreePin7group_mix_rodded_unrodded
TenPin7group_mix_rodded_unrodded

eSolver(1) = eigensolverClass(input);
eSolver(1).solve();
input.subray = 1;
eSolver(2) = eigensolverClass(input);
eSolver(2).solve();
input.subray = 2;
eSolver(3) = eigensolverClass(input);
eSolver(3).solve();

eSolver(1).solution.keff(1)
eSolver(2).solution.keff(1)
eSolver(3).solution.keff(1)