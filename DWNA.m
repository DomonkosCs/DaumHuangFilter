
% time integral of wiener process

SIMU_TIME = 30;
D_T = 0.05;
MONTE_CARLO_RUNS = 1000;

timespan = 0:D_T:SIMU_TIME;
SIGMA = 10;
Q = [1/3*D_T^3,1/2*D_T^2;1/2*D_T^2,D_T]*SIGMA^2;
Q = 0.5*(Q+Q');


x = zeros(MONTE_CARLO_RUNS,numel(timespan));
x_ode = zeros(MONTE_CARLO_RUNS,numel(timespan));

mu_hat = zeros(2,SIMU_TIME/D_T+1);
sigma_hat = zeros(2,SIMU_TIME/D_T+1);

parfor i = 1:MONTE_CARLO_RUNS
    noise = mvnrnd([0;0],Q,numel(timespan))';
    x(i,:) = recursion([0;0],timespan,noise);
    x_ode_full = ode4_noise(@(t,x) odefun(t,x,SIGMA),timespan,[0;0],noise);
    x_ode(i,:) = x_ode_full(:,1);
end

for i = 1:SIMU_TIME/D_T+1
[mu_hat(1,i),sigma_hat(1,i)] = normfit(x(:,i));
[mu_hat(2,i),sigma_hat(2,i)] = normfit(x_ode(:,i));
end
hold on 
plot(timespan,sigma_hat)

function x_1 = recursion(x_init,timespan,noise)
    dT = timespan(2)-timespan(1);
    x(:,1) = x_init;
    F = [1,dT;0,1];
    
    for i = 1:numel(timespan)-1
        %noise = mvnrnd([0;0],Q)';
        x(:,i+1) = F*x(:,i) + noise(:,i);
    end
    x_1 = x(1,:);
end

function dxdt = odefun(t,x,sigma)
    %Q = [0,0;0,sigma^2];
    %noise = mvnrnd([0;0],Q);
    dxdt = zeros(2,1);
    dxdt(1) = x(2);
    dxdt(2) = 0;
end
