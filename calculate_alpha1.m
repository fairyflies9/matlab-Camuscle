function [ alpha1 ] = calculate_alpha1( Ca, Caf, P, V, L )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global LO P0 Lis mu0 mu1 lambda2 alpham alphap alphamax k1 k2 k30 k40 k5 km1 km2 C S 

%if(calculate_vc(Ca, Caf, P, V, L)>0),
    alpha1=alpham;  %assume constant slope
%else
%    alpha1=alphap;
%end
end

