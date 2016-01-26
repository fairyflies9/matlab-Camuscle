function [ lc ] = calculate_lc( Ca, Caf, P, V, L )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to determine the current length of the contractile element, lc
%Inputs
%mu = effective stiffness of the series elastic element
%L = sum of contractile and series elements
%P = current total force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global LO P0 Lis mu0 mu1 lambda2 alpham alphap alphamax k1 k2 k30 k40 k5 km1 km2 C S 

lc = L-P/calculate_mu(Ca, Caf, P, V, L);

end

