%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filename: ca_muscle.m
%units:
% length: mm = millimeters  (milli=10^-3)
% time:   s  = seconds
% force:  mN = milliNewtons (milli=10^-3; Newton=Kg m/s^2)   
% velocity:  mm/s = millimeters per second 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clf
global LO P0 Lis mu0 mu1 lambda2 alpham alphap alphamax k1 k2 k30 k40 k5 km1 km2 C S
LO = 2.94;              % mm, Optimal length 
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

tfinal = 1;                     %final simulation time in seconds
N = 1000000;                    % number of time steps
dt = tfinal/(N-1);
P_0 = 0;                        %Initial force
Ca_0 = 1;                       %Initial concentration free Ca
Caf_0 = 0;                      %Initial concentration filament bound Ca
%init_cond = [Ca_0 Caf_0 P_0];    %Initial conditions
init_cond = [Ca_0 Caf_0];
time = 0:dt:tfinal;             %vector of times to solve ODEs
[T,X] = ode23s(@deriv,time,init_cond);

P = zeros(N,1);
L = LO;              %get length L and velocity V of SE+CE
V = 0;
for(i=1:N-1),
    %P(i+1) = P(i)+dt*(k5*(calculate_Pc(X(i,1), X(i,2), P(i), V, L)-P(i))); 
end
%plot results
figure(1)
plot(T,X)
xlabel('time (sec)');
ylabel('dimensionless concentration');
legend('Ca', 'Caf');

figure(2)
plot(T,P)
xlabel('time (sec)');
ylabel('force (mN)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




