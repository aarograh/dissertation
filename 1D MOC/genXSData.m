function [ xstr_list ] = genXSData( nmats, ngroups, id_mod, id_clad, id_gap, id_fuel, id_gtube, ...
    id_controlgap, id_controlmod, id_control )
%GENXSDATA Summary of this function goes here
%   nmats       - Number of materials
%   ngroups     - Number of energy groups
%     Remaining variables are material IDs, and will be removed eventually

xstr_list = ones(nmats,ngroups);
xstr_list(id_mod:id_fuel,47) = [4.394922836; 0.299908964; 3.44E-05; 1.489061075];
xstr_list(id_gtube:id_control,47) = [0.299908964; 3.44E-05; 4.39E+00; 19.13440383];
xstr_list(id_controlmod,:) = xstr_list(id_mod,:);
xstr_list(id_controlgap,:) = xstr_list(id_controlgap,:);

end

