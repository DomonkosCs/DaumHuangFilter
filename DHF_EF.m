%% 1D test case

% TODO: function descriptions

% plant_function = @(x,k) x/2 + 25*x/(1+ x^2) + 8*cos(1.2*k);
% measurement_function = @(x) x^2/20;
% 

%------- 1D
% state_equation = @(xk_prev,k) xk_prev./2 + 25.*xk_prev./(1 + xk_prev.^2) + 8*cos(1.2*k);
% output_equation = @(xk) xk.^2 ./ 20;
%-------

%------- 2D Pendulum parameters
close all

GRAVITY_SI = 9.81;
LENGTH_METER = 0.8;
TIME_STEP_SEC = 0.01; % correlation between the time step and the noise!
SIMU_TIME_SEC = 5;
a = -3*GRAVITY_SI/2/LENGTH_METER;

%-- TODO!
system_matr_A = [1,TIME_STEP_SEC;0,1]; % discretized, nonlinear
A_nonlin = [0;TIME_STEP_SEC*a];
output_matr_C = [1,0];
INIT_STATE = [deg2rad(-110);0]; % , or ; ?

state_fcn = @(xk_prev) system_matr_A*xk_prev+A_nonlin*sin(xk_prev(1));
output_fcn = @(x) output_matr_C*x;
dim_x = size(system_matr_A,1);
dim_y = size(output_matr_C,1);

% TODO: redundance

determ_fcn = @(x) ([0,1;0,0]*x+[0,0;a,0]*sin(x)); % continuous, nonlinear
STD = 1;
stoch_fcn = @(x) STD*[0;1]; % = noise_gain_matrix_G
%-------
%------- Filter parameters
noise_gain_matrix_G = [0;1];
PROCESS_NOISE_COVAR_Q = generateNoiseCovar([0,1;a,0],[0;1],STD,TIME_STEP_SEC); % TODO
MEASUR_NOISE_VAR_R = 1e-2*eye(dim_y);
PARTICLE_COUNT = 30;

%------- Monte Carlo parameters
MONTE_CARLO_RUNS = 100;
filter_rms = zeros(2,MONTE_CARLO_RUNS,numel(PARTICLE_COUNT));

%-------

%------
if (true)
    
equations = systemEquationsHandler( ...
            state_fcn, output_fcn, dim_x);

[state_array,measure_array] = simulateAndMeasure( ...
    SIMU_TIME_SEC, 0.0001, INIT_STATE, determ_fcn, stoch_fcn, ...
    output_fcn, TIME_STEP_SEC,MEASUR_NOISE_VAR_R);

[ekf_state,dhf_state] = exactFlowFilter(...
            equations,INIT_STATE,PARTICLE_COUNT,SIMU_TIME_SEC,TIME_STEP_SEC,PROCESS_NOISE_COVAR_Q, ...
            MEASUR_NOISE_VAR_R,measure_array);  
figure(1)
hold on
plot([0:TIME_STEP_SEC:SIMU_TIME_SEC],state_array(1,:))
plot([0:TIME_STEP_SEC:SIMU_TIME_SEC],dhf_state(1,:))
plot([0:TIME_STEP_SEC:SIMU_TIME_SEC],ekf_state(1,:))
legend('true','dhf','ekf')
figure(2)
hold on
plot([0:TIME_STEP_SEC:SIMU_TIME_SEC],state_array(2,:))
plot([0:TIME_STEP_SEC:SIMU_TIME_SEC],dhf_state(2,:))
plot([0:TIME_STEP_SEC:SIMU_TIME_SEC],ekf_state(2,:))
legend('true','dhf','ekf')
end

if (false)
    parfor i = 1:MONTE_CARLO_RUNS
        equations = systemEquationsHandler( ...
            state_fcn, output_fcn, dim_x);
        [state_array,measure_array] = simulateAndMeasure( ...
            SIMU_TIME_SEC, 0.0001, INIT_STATE, determ_fcn, stoch_fcn, ...
            output_fcn, TIME_STEP_SEC,MEASUR_NOISE_VAR_R);

        [ekf_state,dhf_state] = exactFlowFilter(...
            equations,INIT_STATE,PARTICLE_COUNT,SIMU_TIME_SEC,TIME_STEP_SEC,PROCESS_NOISE_COVAR_Q, ...
            MEASUR_NOISE_VAR_R,measure_array);  

        filter_rms(:,i) = [rms(state_array(1,:)-ekf_state(1,:)); ...
                           rms(state_array(1,:)-dhf_state(1,:))];
    end
    mean(filter_rms,2)
end

function [EKF_mean_m,filtered_state] = exactFlowFilter( ...
        system_equations, ...
        init_state, ...
        particle_count, ...
        end_time, ...
        delta_time, ...
        process_noise_variance_Q, ...
        measurement_noise_variance_R, ...
        measurement_z)
    
    state_dim = numel(init_state);
    LAMBDA_DIV_POINTS = 11;
    d_lambda = 1/(LAMBDA_DIV_POINTS-1);
    time_div_points = end_time/delta_time+1;
    
    %preallocate
    filtered_state = zeros(state_dim,time_div_points); % to store the DHF filtered states
    state_array = zeros( ...
        state_dim,particle_count,time_div_points);
    EKF_mean_m = zeros(state_dim,time_div_points);
    EKF_covar_P = zeros(state_dim,state_dim,time_div_points);

    %init
    EKF_covar_P(:,:,1) = 10*process_noise_variance_Q;
    state_array(:,:,1) = mvnrnd(init_state,EKF_covar_P(:,:,1),particle_count)'; %'!
    filtered_state(:,1) = mean(state_array(:,:,1),2);
    EKF_mean_m(:,1) = filtered_state(:,1);
    
    for i = 2:time_div_points
        state_array(:,:,i) = system_equations.plant_with_noise( ...
            state_array(:,:,i-1),process_noise_variance_Q); % !!!!!
        filtered_state(:,i) = mean(state_array(:,:,i),2);
        [predicted_EKF_mean_m,predicted_EKF_covar_P] = EKFPrediction( ...
            system_equations, ...
            EKF_mean_m(:,i-1), ...
            EKF_covar_P(:,:,i-1), ...
            process_noise_variance_Q);
        
        for j = 1:LAMBDA_DIV_POINTS
            lambda = j * d_lambda;
            
            [B, b] = calculateSlope( ...
                system_equations, lambda, filtered_state(:,i),...
                filtered_state(:,i), measurement_noise_variance_R, ...
                predicted_EKF_covar_P, measurement_z(i));
            
            for k = 1:particle_count
                dx_dlambda = B * state_array(:,k,i) + b;
                state_array(:,k,i) = state_array(:,k,i) + dx_dlambda * d_lambda; 
            end
            filtered_state(:,i) = mean(state_array(:,:,i),2); % re evaluate the avg
        end
        
        [EKF_mean_m(:,i),EKF_covar_P(:,:,i)] = EKFUpdate( ...
            system_equations, ...
            predicted_EKF_mean_m, ...
            predicted_EKF_covar_P, ...
            measurement_noise_variance_R, ...
            measurement_z(i));
    end
end % TODO: get rid of the for loops -> vectorize

function [B, b] = calculateSlope( ...
        system_equations, ...
        lambda, ...
        currentState, ...
        state_avg, ...
        measurement_noise_variance_R, ...
        predicted_EKF_covar_P, ...
        measurement_z)
    
    linear_measur_eval = system_equations.linear_sensor_matr(currentState); 
    B = -1/2*predicted_EKF_covar_P*linear_measur_eval' ...
        /(lambda*linear_measur_eval*predicted_EKF_covar_P ...
        *linear_measur_eval' + measurement_noise_variance_R) * linear_measur_eval;
    b = (eye(size(B)) - 2*lambda*B) ...
        * ((eye(size(B)) + lambda*B)*predicted_EKF_covar_P*linear_measur_eval'...
        / measurement_noise_variance_R*measurement_z + B*state_avg);
end




