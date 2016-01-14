function [ lambda ] = calculate_lambda( Ca, Caf, P, V, L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to calculate lambda, a parabolic function to describe the
%length-tension curve
%Inputs
%lambda2 = constant to describe how force changes as you move away from
%optimal length.
%lc = length of contractile element
%L0 = Optimal length of SE + CE.
%P = force
%mu = effective stiffness series elastic element.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global LO P0 Lis mu0 mu1 lambda2 alpham alphap alphamax k1 k2 k30 k40 k5 km1 km2 C S 

lc0=2.6;    %lc0 is is the length of the CE at the optimum length, L0
lambda = 1+lambda2*(calculate_lc(Ca, Caf, P, V, L)-lc0)^2;

end

