function [ mu ] = calculate_mu( Ca, Caf, P, V, L )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to calculate the effective stiffness of the series elastic
%element, mu.
%Inputs
%mu0 = Stiffness of SE when Caf=0
%mu1 = Gradient of mu against Caf
%Caf = non-dimensional concentration of Ca-bound filaments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global LO P0 Lis mu0 mu1 lambda2 alpham alphap alphamax k1 k2 k30 k40 k5 km1 km2 C S 

mu = mu0+mu1*Caf;

end

