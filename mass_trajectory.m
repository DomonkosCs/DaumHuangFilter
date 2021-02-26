
eta = 1.64e-4;
[t,y] = ode45(@(t,y) odefun_param(t,y,eta),[0,60],[91500;6100;0.06]);

function dxdt = odefun_param(t,x,eta)
    dxdt = zeros(3,1);
    dxdt(3) = 0;
    dxdt(2) = -x(3)^2*x(2)^2*exp(-eta*x(1));
    dxdt(1) = -x(2);
end

function y = measurement(a,b,x1,sigma,freq)
     y = sqrt(b^2+(x1-a)^2) + normrnd(0,sigma)
end