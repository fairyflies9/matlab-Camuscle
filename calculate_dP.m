function [ dP ] = calculate_dP( Ca, Caf, P, V, L )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to calculate the rate of change of the force
%Inputs
%Caf = nondimensional concentration of filament bound Ca.
%V = shortening velocity of series + contractile elements.
%P = force
%mu = effective stiffness series elastic element.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global LO P0 Lis mu0 mu1 lambda2 alpham alphap alphamax k1 k2 k30 k40 k5 km1 km2 C S 

dP = (calculate_lambda(Ca, Caf, P, V, L)*Caf*(1+calculate_alpha1(Ca, Caf, P, V, L)*V + calculate_alpha1(Ca, Caf, P, V, L)*mu1*P*(calculate_dCaf(Ca, Caf, P, V, L))/((calculate_mu(Ca, Caf, P, V, L))^2))-P)/((1/k5)+calculate_lambda(Ca, Caf, P, V, L)*calculate_alpha1(Ca, Caf, P, V, L)*Caf/calculate_mu(Ca, Caf, P, V, L));

end

