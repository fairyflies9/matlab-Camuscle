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
global LO P0 mu0 mu1 lambda2 alpham alphap k1 k2 k30 k40 k5 C S

LO = 2.94;              % mm, Optimal length 
P0 = 67;                % maximal force mN mm^-2

mu0 = 1;                % Stiffness of SE when Caf=0
mu1 = 23;               % Gradient of mu against Caf 

lambda2 = -20;          % Coefficient of lambda(lc) 
alpham = 0.80;          % Coefficient of alpha(vc) for vc<0, s 
alphap = 2.90;          % Coefficient of alpha(vc) for vc>0, s

k1 = 9;                 % Rate constant, Ca2+ binding in SR, s^-1 
k2 = 50;                % Rate constant, Ca2+ release from SR, s^-1
k30 = 40;               % Rate constant, Ca2+ binding to filaments, s^-1
k40 = 19.4;             % Rate constant, Ca2+ release from filaments, s^-1
k5 = 100;               % Rate constant, transfer of force from CE to SE, s^-1

C = 2;                  % total dimensionless Ca2+ concentration
S = 6;                 % total dimensionless concentrations of sarcoplasmic-reticular binding sites



tfinal = 1;                     %final simulation time in seconds
N = 100000;                    % number of time steps
dt = tfinal/(N-1);
time = 0:dt:tfinal;             %vector of times to solve ODEs

% --- Initial conditions
P_0 = 0;                        %Initial force
Ca_0 = 0.2;                       %Initial concentration free Ca
Caf_0 = 0.9;                      %Initial concentration filament bound Ca
% ---

V = 0;
% L = LO;
lc0=2.6;    % lenght of CE at optimum length L0
L = 2.6;
% --- functions
mu = @(Caf) mu0+mu1*Caf;
alpha1 = alpham;

% alpha = @(vc) 1 + alpham*vc;

P = zeros(N,1);
lc = zeros(N,1);
Ca = zeros(N,1);
Caf = zeros(N,1);

P(1) = P_0;
Ca(1) = Ca_0;
Caf(1) = Caf_0;

lc(1) = L - P(1)./mu(Caf(1));
lambda = P0*(1+lambda2*(lc(1) - lc0).^2);

for i=1:N-1;
    
    if (dt*i)<0.1
        % on
        tmp = k1*(C-Ca(i)-Caf(i));
    else
        tmp = k2*Ca(i)*(C-S-Ca(i)-Caf(i));
    end
    
    RHS_Ca = (k40*Caf(i)-k30*Ca(i))*(1-Caf(i)) + tmp;
    RHS_Caf = -(k40*Caf(i) - k30*Ca(i))*(1-Caf(i));
    
    RHS_P = (lambda*Caf(i)*( 1+alpha1*V+alpha1*mu1*P(i)*RHS_Caf/mu(Caf(i))^2 )-P(i))...
        / (1/k5 + lambda*lc(i)*alpha1*Caf(i)/mu(Caf(i)));
              
    % compute at next time step
    Ca(i+1) = Ca(i) + RHS_Ca*dt;
    Caf(i+1) = Caf(i) + RHS_Caf*dt;
    
    P(i+1) = P(i) + RHS_P*dt;
    
    lc(i+1) = L - P(i+1)./mu(Caf(i+1));
    lambda = P0*(1+lambda2*(lc(i+1) - lc0).^2);
    
end


%plot results
figure(1)
plot(time,Ca,time,Caf)
xlabel('time (sec)');
ylabel('dimensionless concentration');
legend('Ca', 'Caf');

figure(2)
plot(time,P,time,lc)
xlabel('time (sec)');
ylabel('force (mN)');
legend('force','lc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




