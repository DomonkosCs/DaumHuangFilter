% Euler-Maruyama method for SDE in case of a 1DoF pendulum
function [xs,ts] = stochastic_pendulum_em( ...
        simu_time, ...
        dt, ...
        x_init, ...
        determ_fcn, ...
        stoch_fcn)
    
    dim = numel(x_init);
    N = simu_time/dt+1;
    xs = zeros(dim,N);
    xs(:,1) = x_init;
    ts = 0:dt:simu_time;
    for i = 2:N
        dW = sqrt(dt) * randn(1);
        x = xs(:,i-1);
        xs(:,i) =  x + determ_fcn(x)*dt + stoch_fcn(x) * dW;
    end

end