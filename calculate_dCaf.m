function [ dCaf ] = calculate_dCaf( Ca, Caf, P, V, L )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function that calculates the rate of change of filament bound Ca2+.
%Inputs
%k3 = Rate constant, Ca2+ binding to filaments
%k4 = Rate constant, Ca2+ release from filaments
%Ca = nondimensional concentration of free Ca2+
%Caf = nondimensional concentration of filament bound Ca2+
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global LO P0 Lis mu0 mu1 lambda2 alpham alphap alphamax k1 k2 k30 k40 k5 km1 km2 C S 

dCaf = -(k40*Caf-k30*Ca)*(1-Caf);

end

