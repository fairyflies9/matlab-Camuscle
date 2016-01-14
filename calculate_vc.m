function [ vc ] = calculate_vc( Ca, Caf, P, V, L )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to calculate the shortening velocity of the contractile element.
%Inputs
%V = shortening velocity of entire element.
%dP = rate of change of force
%mu = effective stiffness series elastic element
%mu1 = Gradient of mu against Caf
%P = total force
%dCaf = rate of change of non-dimensional concentration of Ca-bound filaments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global LO P0 Lis mu0 mu1 lambda2 alpham alphap alphamax k1 k2 k30 k40 k5 km1 km2 C S 

vc = V - calculate_dP(Ca, Caf, P, V, L)/calculate_mu(Ca, Caf, P, V, L) + mu1*P*calculate_dCaf(Ca, Caf, P, V, L)/((calculate_mu(Ca, Caf, P, V, L))^2);


end

