% Euler-Maruyama method for SDE in case of a 1DoF pendulum. Nonlinear.
tic
t_bounds = [0,5]; % sec
N = 15001;
dt = (t_bounds(2) - t_bounds(1)) / (N-1);
ts = linspace(t_bounds(1),t_bounds(2),N);

x_init = [deg2rad(170);0];
mc_run = 100;

g = 9.81;
l = 0.2;
a = -3*g/2/l;
sigma_0 = 0.1;
sigma = [0; sigma_0];
determ_fcn = @(x) ([0,1;0,0]*x+[0,0;a,0]*sin(x));
stoch_fcn = @(x) [0;sigma_0];

x_mc = zeros(mc_run,numel(ts));
parfor i = 1:mc_run
    xs = emPendulum(N,x_init,ts,dt,determ_fcn,stoch_fcn);
    x_mc(i,:) = xs(1,:); 
end
[~,sigma_hat] = inspectNoise(ts,x_mc);
%plot(ts,sigma_hat)
toc
%plot(ts,sigma_hat)
function [mu_hat,sigma_hat] = inspectNoise(ts,x)
    mu_hat = zeros(1,numel(ts));
    sigma_hat = zeros(1,numel(ts));
    for i = 1:numel(ts)
        [mu_hat(1,i),sigma_hat(1,i)] = normfit(x(:,i));
    end
end

function xs = emPendulum(N,x_init,ts,dt,determ_fcn,stoch_fcn)
    xs = zeros(2,N);
    xs(:,1) = x_init;
    for i = 2:numel(ts)
        dW = sqrt(dt) * randn(1);
        x = xs(:,i-1);
        %xs(:,i) =  x + ([0,1;0,0]*x+[0,0;a,0]*sin(x))*dt + sigma * dW;
        xs(:,i) =  x + determ_fcn(x)*dt + stoch_fcn(x) * dW;
    end
end