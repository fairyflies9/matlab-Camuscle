%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filename: ca_muscle_v3.m
%units:
% length: mm = millimeters  (milli=10^-3)
% time:   s  = seconds
% force:  mN = milliNewtons (milli=10^-3; Newton=Kg m/s^2)
% velocity:  mm/s = millimeters per second
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clf
global L0 P0 mu0 mu1 lambda2 alpham alphap k1 k2 k30 k40 k5 C S

L0 = 2.94;              % mm, Optimal length
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
    
km1 = 15;               % Rate m shortening
km2 = 10;               % Rate m returning to 1

C = 200;                  % total dimensionless Ca2+ concentration
S = 200;                  % total dimensionless concentrations of sarcoplasmic-reticular binding sites

% --- more constants
V = 0;                  % rate of change of L
lc0=2.6/L0;             % length of CE at optimum length L0 / non-dimensionalized by L0
L = 2.94/L0;            % maximal length of muscle / non-dimensionalized by L0

% --- function
mu = @(Caf) mu0+mu1*Caf;

% ---
tfinal = 1;                    %final simulation time in seconds
N = 100000;                    % number of time steps
dt = tfinal/(N-1);
time = 0:dt:tfinal;            %vector of times to solve ODEs

% --- Initial conditions
P = zeros(N,1);
lc = zeros(N,1);
Ca = zeros(N,1);
Caf = zeros(N,1);
m = zeros(N,1);
vc = zeros(N,1);

P(1) = 0;           % inital force
Ca(1) = 0;          % inital free-calcium
Caf(1) = 0;         % initial bound calcium to filaments
m(1) = 1;           % initial condition - m controls binding and release rates

lc(1) = L - P(1)./mu(Caf(1));
lambda = (1+lambda2*(lc(1) - lc0).^2);
alpha1 = alphap;        % assuming vc>0 initially

% generate the spike-train forcing
% nu = 5 frequency
% tau = 20ms decay time
% k1b = kbar(tfinal,dt,k1/5,10,0.2);
k1b = k1*kbar2(tfinal,dt,k1/5,5,0.2);

% --- time stepping with forward Euler
for i=1:N-1;
    
    % calculate all RHS of DEQs at current time step (i)
    tmp = k2*Ca(i)*(C-S-Ca(i)-Caf(i));  % uptake of calcium always happens

    % account for stimulus
    tmp = tmp + k1b(i)*(C-Ca(i)-Caf(i));
    
    k3 = k30/sqrt(m(i));
    k4 = k40/sqrt(m(i));
    
    % removed weird non-biological term
     RHS_Ca = k4*Caf(i) - k3*Ca(i)*(1-Caf(i)) + tmp;
     RHS_Caf = -k4*Caf(i) + k3*Ca(i)*(1-Caf(i));
    
%     RHS_Ca = (k4*Caf(i)-k3*Ca(i))*(1-Caf(i)) + tmp;
%     RHS_Caf = -(k4*Caf(i) - k3*Ca(i))*(1-Caf(i));
    
    RHS_P = (lambda*Caf(i)*( 1+alpha1*V+alpha1*mu1*P(i)*RHS_Caf/mu(Caf(i))^2 )-P(i))...
        / (1/k5 + lambda*lc(i)*alpha1*Caf(i)/mu(Caf(i)));
    
    vc(i) = V - RHS_P/mu(Caf(i));
    
    if vc(i)<0
        % when shortening
        RHS_m = -km1*P(i)*vc(i);
        alpha1 = alpham;
    else
        RHS_m = -km2*(m(i)-1);
        alpha1 = alphap;
    end
    
    % compute at next time step
    Ca(i+1) = Ca(i) + RHS_Ca*dt;
    Caf(i+1) = Caf(i) + RHS_Caf*dt;
    
    P(i+1) = P(i) + RHS_P*dt;
    
    m(i+1) = m(i) + RHS_m*dt;
    
    lc(i+1) = L - P(i+1)./mu(Caf(i+1));
    lambda = (1+lambda2*(lc(i+1) - lc0).^2);
    
end
RHS_P = (lambda*Caf(N)*( 1+alpha1*V+alpha1*mu1*P(N)*RHS_Caf/mu(Caf(N))^2 )-P(N))...
    / (1/k5 + lambda*lc(N)*alpha1*Caf(N)/mu(Caf(N)));

vc(N) = V - RHS_P/mu(Caf(N));

% --- plot results
figure(3);clf;
plot(time,Ca,time,Caf,'linewidth',2);hold on
plot(time,k1b,'k--','linewidth',2)
xlabel('time (sec)');
ylabel('dimensionless concentration');
legend('Ca', 'Caf','k1');
ylim([0 3])
set(gca,'fontsize',18)

figure(4)
plot(time,P,time,lc,time,vc,time,m,'linewidth',2)
xlabel('time (sec)');
ylabel('force (mN)');
legend('force','lc','vc','m')
set(gca,'fontsize',18)

%figure(2)
%plot(time,C-Ca-Caf, time,C-S-Ca-Caf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




