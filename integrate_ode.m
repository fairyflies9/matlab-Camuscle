% script containing integration scheme

% --- function
mu = @(Caf) mu0+mu1*Caf;

% ---
tfinal = 1;                    %final simulation time in seconds
N = 10000;                    % number of time steps
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

% --- time stepping with forward Euler
for i=1:N-1;
    
    % calculate all RHS of DEQs at current time step (i)
    tmp = k2*Ca(i)*(C-S-Ca(i)-Caf(i));  % uptake of calcium always happens
    if (dt*i)<1
        % stimulus on
        tmp = tmp + k1*(C-Ca(i)-Caf(i));
    end
    
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