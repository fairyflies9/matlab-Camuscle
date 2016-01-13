
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Inputs for Ca2+ - muscle model for lamprey from Williams et al.%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L0 = 2.94;              % mm, Optimal length 
P0 = 67;                % mN*mm^2, Maximal tetanic force 
Lis = 2.70;             % mm, Length in situ 
mu0 = 1;                % Stiffness of SE when Caf=0
mu1 = 23;               % Gradient of mu against Caf 
lambda2 = -20;          % Coefficient of lambda(lc) 
alpham = 0.80;          % Coefficient of alpha(vc) for vc<0, s 
alphap = 2.90;          % Coefficient of alpha(vc) for vc>0, s
alphamax = 1.8;         % Maximum value for alpha(vc)
k1 = 9;                 % Rate constant, Ca2+ binding in SR, s^-1 
k2 = 50;                % Rate constant, Ca2+ release from SR, s^-1
k30 = 40;               % Rate constant, Ca2+ binding to filaments, s^-1
k40 = 19.4;             % Rate constant, Ca2+ release from filaments, s^-1
k5 = 100;               % Rate constant, transfer of force from CE to SE, s^-1
km1 = 15;               % Rate constant of m increase, s^-1
km2 = 10;               % Rate constant of m decrease, s^-1
C = 2;                  % total dimensionless Ca2+ concentration
S = 6;                  % total dimensionless concentrations of sarcoplasmic-reticular binding sites