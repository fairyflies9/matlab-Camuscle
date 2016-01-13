function dy=deriv(t,x)
%Function to solve system of differential equations for combined muscle
%model. 
global LO P0 Lis mu0 mu1 lambda2 alpham alphap alphamax k1 k2 k30 k40 k5 km1 km2 C S 

dy = zeros(2,1);                %initialize as column vector
L = LO;              %get length L and velocity V of SE+CE
V = 0;

%Solve system of differential equations
if(stimulation(t)>0),           %check if stimulus on or off
    dy(1) = (k40*x(2)-k30*x(1))*(1-x(2))+k1*(C-x(1)-x(2));
else
    dy(1) = (k40*x(2)-k30*x(1))*(1-x(2))+k2*(x(1)*(C-S-x(1)-x(2)));
end
dy(2) = -(k40*x(2)-k30*x(1))*(1-x(2));
end
