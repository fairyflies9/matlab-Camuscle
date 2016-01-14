function [ alpha ] = calculate_alpha(Ca, Caf, P, V, L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to calculate alpha, which describes the force-velocity
%relationship. Note that alpha is a piecewise linear
%function (approximating a hyperbola).
%Inputs
%vc = shortening velocity of contractile element.
%alpham = coefficient for shortening
%alphap = coefficient for lengthening
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global LO P0 Lis mu0 mu1 lambda2 alpham alphap alphamax k1 k2 k30 k40 k5 km1 km2 C S 

%if(calculate_vc(Ca, Caf, P, V, L)<0),
    alpha = 1 + alpham*calculate_vc(Ca, Caf, P, V, L);  %assume length-tension curve has constant slope
%else
%    alpha = 1 + alphap*calculate_vc(Ca, Caf, P, V, L);
%end

end

