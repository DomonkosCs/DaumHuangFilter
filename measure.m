
simu_time = 5;
dt = 0.0001;
x_init = [deg2rad(30);0];
determ_fcn = @(x) ([0,1;0,0]*x+[0,0;-70,0]*sin(x));
stoch_fcn = @(x) [0;1];
measure_fcn = @(x) [1,0]*x;
dt_sample = 0.003443667;
[x,tx] = stochastic_pendulum_em(simu_time,dt,x_init,determ_fcn,stoch_fcn);
simu_time = tx(end);

y = measureState(x,tx,dt_sample,measure_fcn,0.01);

plot(tx,x(1,:))
hold on
plot([0:dt_sample:simu_time],y,'ro')

function y = measureState(x,tx,dt_sample,measure_fcn,R)
    simu_time = tx(end);
    t_sample = 0:dt_sample:simu_time;
    x1_sample = interp1(tx,x(1,:),t_sample);
    x2_sample = interp1(tx,x(2,:),t_sample);
    x_sample = [x1_sample; x2_sample];
    y = measure_fcn(x_sample) + sqrt(R)*randn(size(x1_sample));
end