% Euler-Maruyama method for SDE in case of a 1DoF pendulum
tic
t_bounds = [0,5]; % sec
N = 4001;
dt = (t_bounds(2) - t_bounds(1)) / (N-1);
ts = linspace(t_bounds(1),t_bounds(2),N);
pd = makedist('Normal',0,sqrt(dt));

x_init = [deg2rad(5);0];
mc_run = 1;

g = 9.81;
l = 0.2;
A = [0, 1; ...
    -3*g/2/l, 0];
sigma_0 = 0.4;
sigma = [0; sigma_0];

x_mc = zeros(mc_run,numel(ts));
for i = 1:mc_run
    xs = emPendulum(N,x_init,ts,pd,dt,A,sigma);
    x_mc(i,:) = xs(1,:); 
end
[~,sigma_hat] = inspectNoise(ts,x_mc);
%plot(ts,sigma_hat)
toc
function [mu_hat,sigma_hat] = inspectNoise(ts,x)
    mu_hat = zeros(1,numel(ts));
    sigma_hat = zeros(1,numel(ts));
    for i = 1:numel(ts)
        [mu_hat(1,i),sigma_hat(1,i)] = normfit(x(:,i));
    end
end

function xs = emPendulum(N,x_init,ts,pd,dt,A,sigma)
    xs = zeros(2,N);
    xs(:,1) = x_init;
    for i = 2:numel(ts)
        dW = random(pd);
        x = xs(:,i-1);
        xs(:,i) =  x + A*x*dt + sigma*dW;
    end
end