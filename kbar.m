function f = kbar(tfinal,dt,nu,tau)
global k1
% dt time step size
% nu frequency of arriving spikes
% tau decay constant
% f vector value of kbar

nstepmax = ceil(tfinal/dt)+1;
f = zeros(nstepmax,1);

% spike times
ts = [0:1/nu:tfinal];
nt = 2;     % the next spike time to check

% for computing speed
expt = exp(-dt/tau);

time = 0;   % current time
f(1) = k1;  % initial condition
for nstep = 2:nstepmax
    if (ts(nt)>time) && (ts(nt)<time+dt)
        % spike to be used
        f(nstep) = f(nstep-1)*expt + k1*exp(-(time+dt-ts(nt))/tau);
        nt = nt+1;
    else
        f(nstep) = f(nstep-1)*expt;
    end
    
    time = time+dt;
end

end % end function