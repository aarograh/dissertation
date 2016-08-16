function [ mu, mucos, weights ] = setupQuad( npol )
%SETUPQUAD Summary of this function goes here
%   npol - number of polar angles requested by user

if (npol == 1)
    mu = pi/4;
    mucos = cos(mu);
    weights = 2*pi;
elseif (npol == 2)
elseif (npol == 4)
elseif (npol == 8)
elseif (npol == 16)
end

end

