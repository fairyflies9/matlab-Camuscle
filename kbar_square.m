function f = kbar_square(tfinal,dt,a,nu,b)
% dt time step size
% a height to add each spike time
% nu frequency of arriving spikes
% b square wave width
% f vector value of kbar

nstepmax = ceil(tfinal/dt)+1;
f = zeros(nstepmax,1);

% spike times
ts = [0:1/nu:tfinal];
ns = round(b/dt); % time stimulus on for

for k=1:length(ts)
    ind1 = floor(ts(k)/dt) + 1;
    ind2 = min(ind1+ns,nstepmax);
    
    f(ind1:ind2) = a;
end