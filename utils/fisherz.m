function z=fisherz(r);

% Z=fisherz(R)
%
% this function performs a fisher Z-transform of vector R and output the
% transformed vector Z. This function is used to modify for example
% correlations in [-1;+1] into a near gaussian population. It could be used
% to transform any uniform population in the interval [-1;+1] into a
% gaussian population. The inverse operation is done by function inverse_fisherz.m
%
% Vincent MORON
% feb. 2005

z=0.5*log((1+r)./(1-r));