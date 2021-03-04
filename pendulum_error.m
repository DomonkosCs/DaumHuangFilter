
g = 9.81;
l = 0.8;
a = -3*g/2/l;
A = [0,1;a,0];
G = [0;1];

SIMU_TIME = 10;
D_T = 0.04 ;
timespan = 0:D_T:SIMU_TIME;
INIT_DEG = 100;

SIGMA = 1;
MONTE_CARLO_RUNS = 1000;
Q = generateNoiseCovar(A,G,SIGMA,D_T); % linearized!

x_noisy = zeros(MONTE_CARLO_RUNS,SIMU_TIME/D_T+1);
mu_hat = zeros(1,SIMU_TIME/D_T+1);
sigma_hat = zeros(1,SIMU_TIME/D_T+1);

%[t,x] = ode23(@(t,x) odefun_param(t,x,a,A,G,sigma,dT),[0,10],[deg2rad(179);0]) ;
x_noiseless = ode4(@(t,x) odefun(t,x,a,Q,true),timespan,[deg2rad(INIT_DEG);0]);


parfor i = 1:MONTE_CARLO_RUNS
    x_noisy_full = ode4(@(t,x) odefun(t,x,a,Q,false),timespan,[deg2rad(INIT_DEG);0]);
    x_noisy(i,:) = x_noisy_full(:,1);
end

for i = 1:SIMU_TIME/D_T+1
[mu_hat(1,i),sigma_hat(1,i)] = normfit(x_noisy(:,i));
end
hold on
plot(timespan,sigma_hat)

function dxdt = odefun(t,x,a,Q,noiseless,pd)
    if (noiseless)
        noise = [0;0];
    else
        noise = mvnrnd([0;0],Q); % decompose
    end
    dxdt = zeros(2,1);
    dxdt(1) = 0 * x(1) + 1 * x(2) + noise(1);
    dxdt(2) = a * x(1) + 0 * x(2) + noise(2);
end

% a Q mátrix tuti jó. 
% 
