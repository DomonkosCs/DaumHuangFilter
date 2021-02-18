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
GRAVITY_SI = 9.81;
LENGTH_METER = 0.8;
TIME_STEP_SEC = 0.001; % az időépés és a zaj között összefüggés van!
SIMU_TIME_SEC = 5;
system_matr_A = [1,TIME_STEP_SEC;0,1];
A_nonlin = [0;-TIME_STEP_SEC*3*GRAVITY_SI/2/LENGTH_METER];
output_matr_C = [1,0];
INIT_STATE = [deg2rad(30),0];

state_equation = @(xk_prev) system_matr_A*xk_prev+A_nonlin*sin(xk_prev(1));
output_equation = @(xk) output_matr_C*xk;
dim_x = size(system_matr_A,1);
dim_y = size(output_matr_C,1);
%-------
%------- Filter parameters
PROCESS_NOISE_VAR_Q = 1e-4*eye(dim_x);
MEASUR_NOISE_VAR_R = 1e-2*eye(dim_y);
NOISE_GAIN_G = [1,0;0.2,0];
PARTICLE_COUNT = 15;

%------- Monte Carlo parameters
MONTE_CARLO_RUNS = 1000;
%-------

true_state = zeros(2,SIMU_TIME_SEC/TIME_STEP_SEC+1);
[time,true_state(1,:),true_state(2,:)] = pendulumSimulation( ...
    LENGTH_METER, TIME_STEP_SEC, SIMU_TIME_SEC, INIT_STATE);

equations = systemEquationsHandler( ...
    state_equation, output_equation, dim_x, dim_y);
[state_array,measur_array] = simulateProcess( ...
    equations, INIT_STATE, TIME_STEP_SEC, SIMU_TIME_SEC, NOISE_GAIN_G, ...
    PROCESS_NOISE_VAR_Q, MEASUR_NOISE_VAR_R);

[dhf_state,ekf_state] = exactFlowFilter(...
    equations,INIT_STATE,PARTICLE_COUNT,SIMU_TIME_SEC,TIME_STEP_SEC,NOISE_GAIN_G,PROCESS_NOISE_VAR_Q, ...
    MEASUR_NOISE_VAR_R,measur_array);

scatter(true_state(1,:),true_state(2,:))
hold on
scatter(state_array(1,:),state_array(2,:))
figure(2)
plot(time, true_state(1,:))
hold on
plot(time,state_array(1,:))
plot(time,measur_array)
plot(time,dhf_state(1,:))
plot(time,ekf_state(1,:))
legend('true state','noisy state','measurement','dhf-ef','ekf');

% [dhf_filtered_state,ekf_filtered_state] = exactFlowFilter(equations,kParticleCount,kMaxProcessTime,kDeltaTime,...
%             kProcessNoiseVarianceQ,kMeasurementNoiseVarianceR,measurement_array);
% plot(dhf_filtered_state)
% hold on
% plot(state_array)
% plot(ekf_filtered_state)
% legend('dhf','true','ekf')

% filter_RMSE = zeros(kMonteCarloRuns,kMaxProcessTime);

% for i=200
%     parfor j = 1:kMonteCarloRuns
%         [state_array,measurement_array] = ...
%             simulateProcess(equations,kInitState,kDeltaTime,i, ...
%             kProcessNoiseVarianceQ,kMeasurementNoiseVarianceR);
%         [dhf_filtered_state,ekf_filtered_state] = exactFlowFilter(equations,kParticleCount,i,kDeltaTime,...
%             kProcessNoiseVarianceQ,kMeasurementNoiseVarianceR,measurement_array);
%         filter_RMSE(j,i) = sqrt(mean((state_array-ekf_filtered_state).^2));
%     end
% end
% filter_RMSE_mean = mean(filter_RMSE);



function systemEquations = systemEquationsHandler( ...
        state_equation, ...
        output_equation, ...
        dim_x, ...
        dim_y)

    % anonym. -> symbolic for Jacobian
    xk_prev_array = sym('xk_prev',[dim_x,1]);
    xk_array = sym('xk',[dim_x,1]);
    f = state_equation(xk_prev_array);
    h = output_equation(xk_array);
    
    white_noise = @(covar,dim,num_of_states) mvnrnd(zeros(dim,1), ...
        covar,num_of_states)';
    plant_with_noise = @(xk_prev,G,Q) state_equation(xk_prev) ...
        + G*white_noise(Q,dim_x,size(xk_prev,2));
    sensor_with_noise = @(xk,R) output_equation(xk) ...
        + white_noise(R,dim_y,size(xk,2));
    
    linear_system_matr = matlabFunction( ...
        jacobian(f,xk_prev_array) ...
        + 1e-100*xk_prev_array,'vars',{xk_prev_array}); % !!!
    linear_sensor_matr = matlabFunction( ...
        jacobian(h,xk_array) ...
        + 1e-100*sum(xk_array),'vars',{xk_array}); % !!!
    systemEquations.plant_with_noise = plant_with_noise;
    systemEquations.plant_noiseless = state_equation;
    systemEquations.sensor_with_noise = sensor_with_noise;
    systemEquations.sensor_noiseless = output_equation;
    systemEquations.linear_system_matr = linear_system_matr;
    systemEquations.linear_sensor_matr = linear_sensor_matr;
end

function [state_array,measurement_array] = simulateProcess( ...
        system_equations, ...
        init_state, ...
        delta_time, ...
        end_time, ...
        noise_gain_G, ...
        process_noise_variance_Q, ...
        measurement_noise_variance_R)
 
    state_array = generateState(system_equations,init_state, ...
        delta_time,end_time,noise_gain_G,process_noise_variance_Q);
    measurement_array = system_equations.sensor_with_noise(state_array, ...
        measurement_noise_variance_R);

end

function state_x = generateState( ...
        system_equations, ...
        init_state, ...
        delta_time, ...
        end_time, ...
        noise_gain_G, ...
        process_noise_variance_Q)
    
    kStartTime = 0;
    time_eval_points = kStartTime:delta_time:end_time;
    state_x = zeros(numel(init_state),numel(time_eval_points));
    
    state_x(:,1) = init_state;
    for i = 1:end_time/delta_time
        current_state = state_x(:,i);
        state_x(:,i+1) = system_equations.plant_with_noise( ...
            current_state,noise_gain_G,process_noise_variance_Q);
    end
end

function [EKF_mean_m,filtered_state] = exactFlowFilter( ...
        system_equations, ...
        init_state, ...
        particle_count, ...
        end_time, ...
        delta_time, ...
        noise_gain_G, ...
        process_noise_variance_Q, ...
        measurement_noise_variance_R, ...
        measurement_z)
    
    state_dim = numel(init_state);
    LAMBDA_DIV_POINTS = 11;
    INIT_DISTR_COVAR = 1*eye(state_dim);
    
    d_lambda = 1/(LAMBDA_DIV_POINTS-1);
    time_div_points = end_time/delta_time+1;

    filtered_state = zeros(state_dim,time_div_points);
    
    % 1st dim: state vector; 2nd dim: particles
    % 3rd dim: div in lambda; 4th dim: division in time 
    % TODO: not efficient to store all the variables in a giant matrix
    state_array = zeros( ...
        state_dim,particle_count,LAMBDA_DIV_POINTS,time_div_points);
    
    predicted_EKF_mean_m = zeros(state_dim,time_div_points);
    predicted_EKF_covar_P = zeros(state_dim,state_dim,time_div_points);
    EKF_mean_m = zeros(state_dim,time_div_points);
    EKF_covar_P = zeros(state_dim,state_dim,time_div_points);
    
    % initialization from a gaussian, centered around the initial cond (if known)
    EKF_covar_P(:,:,1) = INIT_DISTR_COVAR;
    state_array(:,:,1,1) = mvnrnd(init_state,EKF_covar_P(:,:,1),particle_count)'; %'!
    
    filtered_state(:,1) = mean(state_array(:,:,1,1),2);
    EKF_mean_m(:,1) = filtered_state(:,1);
    
    % iterate through time
    for i = 2:time_div_points
        state_array(:,:,1,i) = system_equations.plant_with_noise( ...
            state_array(:,:,1,i-1),noise_gain_G,process_noise_variance_Q);
        filtered_state(:,i) = mean(state_array(:,:,1,i),2);
        [predicted_EKF_mean_m(:,i),predicted_EKF_covar_P(:,:,i)] = EKFPrediction( ...
            system_equations, ...
            EKF_mean_m(:,i-1), ...
            EKF_covar_P(:,:,i-1), ...
            noise_gain_G, ...
            process_noise_variance_Q);
        % iterate through the pseudo time lambda
        for j = 1:LAMBDA_DIV_POINTS
            lambda = j * d_lambda;
            % iterate through the particles
            for k = 1:particle_count
                dx_dlambda = calculateSlope(system_equations,lambda,state_array(:,k,j,i),filtered_state(:,i), ...
                    measurement_noise_variance_R,predicted_EKF_covar_P(:,:,i),measurement_z(i));
                state_array(:,k,j+1,i) = state_array(:,k,j,i) + dx_dlambda * d_lambda; 
            end
            filtered_state(:,i) = mean(state_array(:,:,j+1,i),2); % re evaluate the avg
        end
        [EKF_mean_m(:,i),EKF_covar_P(:,:,i)] = EKFUpdate(system_equations,predicted_EKF_mean_m(:,i), ...
            predicted_EKF_covar_P(:,:,i),measurement_noise_variance_R,measurement_z(i));
    end
end % TODO: get rid of the for loops -> vectorize

function dx_dlambda = calculateSlope( ...
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
    b = (eye(size(B)) + 2*lambda*B) ...
        * ((eye(size(B)) + lambda*B)*predicted_EKF_covar_P*linear_measur_eval'...
        / measurement_noise_variance_R*measurement_z + B*state_avg);
    
    dx_dlambda = B*currentState + b;
end

function [predicted_EKF_mean_m,predicted_EKF_covar_P] = EKFPrediction( ...
        system_equations, ...
        EKF_mean_m, ...
        EKF_covar_P, ...
        noise_gain_G, ...
        process_noise_variance_Q)
    
    predicted_EKF_mean_m = system_equations.plant_noiseless(EKF_mean_m);
    linear_plant_eval = system_equations.linear_system_matr(EKF_mean_m);
    predicted_EKF_covar_P = linear_plant_eval*EKF_covar_P*linear_plant_eval' ...
        +  noise_gain_G*process_noise_variance_Q*noise_gain_G'; 
    % TODO: make sure that the linearized system matrix is evaluated at the
    % EKF_mean_m point.
end

function [EKF_mean_m,EKF_covar_P] = EKFUpdate( ...
        system_equations, ...
        predicted_EKF_mean_m, ... 
        predicted_EKF_covar_P, ...
        measurement_noise_variance_R, ...
        measurement_z)
    
    linear_measur_eval = system_equations.linear_sensor_matr(predicted_EKF_mean_m);
    kalman_gain_K = predicted_EKF_covar_P * linear_measur_eval' ...
        /(linear_measur_eval*predicted_EKF_covar_P*linear_measur_eval' ...
        + measurement_noise_variance_R);
    EKF_mean_m = predicted_EKF_mean_m + kalman_gain_K ...
        * (measurement_z - system_equations.sensor_noiseless(predicted_EKF_mean_m));
    EKF_covar_P = ...
        (eye(size(predicted_EKF_covar_P)) - kalman_gain_K*linear_measur_eval) ...
        * predicted_EKF_covar_P;
    
    % z is only a scalar. When the measurement error is weighted according
    % to the Kalman gain, the measurement also modifies the non measured
    % state element. Problematic?
end

