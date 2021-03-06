function [ Pc ] = calculate_Pc( Ca, Caf, P, V, L )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to calculate normalized force generated by the contractile
%element.
%Inputs
%lambda = function that describes the length-tension relationship
%alpha = function describes the force-velocity relationship.
%Caf = non-dimensional concentration of Ca-bound filaments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global LO P0 Lis mu0 mu1 lambda2 alpham alphap alphamax k1 k2 k30 k40 k5 km1 km2 C S 

Pc = calculate_lambda(Ca, Caf, P, V, L)*(1+calculate_alpha1(Ca, Caf, P, V, L)*(V-calculate_dP(Ca, Caf, P, V, L)/calculate_mu(Ca, Caf, P, V, L)+mu1*P*calculate_dCaf(Ca, Caf, P, V, L)/(calculate_mu(Ca, Caf, P, V, L))^2))*Caf;

end

