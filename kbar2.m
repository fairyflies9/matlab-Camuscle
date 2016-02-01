function g = kbar2(tfinal,dt,a,nu,tau)
% dt time step size
% a height to add each spike time
% nu frequency of arriving spikes
% tau decay constant
% f vector value of kbar

nstepmax = ceil(tfinal/dt)+1;
f = zeros(nstepmax,1);
g = zeros(nstepmax,1);

% spike times
ts = [0:1/nu:tfinal];
nt = 2;     % the next spike time to check

% for computing speed
expt = exp(-dt/tau);

time = 0;   % current time
f(1) = a;  % initial condition
g(1) = 0;

for nstep = 2:nstepmax
    if (ts(nt)>time) && (ts(nt)<time+dt)
        % spike to be used
        f(nstep) = f(nstep-1)*expt + a*exp(-(time+dt-ts(nt))/tau);
        nt = nt+1;
    else
        f(nstep) = f(nstep-1)*expt;
    end
    
    g(nstep) = g(nstep-1) + (f(nstep-1)*(1-g(nstep-1)) - g(nstep-1)/tau)*dt;
    
    time = time+dt;
end

end % end function